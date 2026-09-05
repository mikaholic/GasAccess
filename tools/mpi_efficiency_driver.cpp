#include "efficiency_fixture.hpp"
#include "gasaccess/distributed_deposition_updater.hpp"
#include "gasaccess/distributed_exterior_classifier.hpp"
#include "gasaccess/mpi_gas_grid.hpp"
#include "gasaccess/spparks_adapter.hpp"
#include "mock_spparks.hpp"
#include "mock_spparks_fixture.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

namespace {

using gasaccess::Atom;
using gasaccess::AtomView;
using gasaccess::DistributedDepositionUpdateResult;
using gasaccess::DistributedDepositionUpdater;
using gasaccess::DistributedExteriorClassifier;
using gasaccess::DistributedGasAccessibilityQuery;
using gasaccess::DistributedGasGrid;
using gasaccess::GasState;
using gasaccess::GridDimensions;
using gasaccess::MpiDecompositionSpec;
using gasaccess::Point3;
using gasaccess::SpparksAdapterConfig;
using gasaccess::SpparksAtomBuffer;
using gasaccess::VoxelCoord;
using gasaccess::testing::EfficiencyRepairCase;
using gasaccess::testing::EfficiencyRepairFixture;
using gasaccess::testing::MockScenario;
using gasaccess::testing::MockSpparksApp;
using gasaccess::testing::MockSpparksDomain;

enum class Operation {
    Initialization,
    Query,
    Repair
};

struct Options {
    Operation operation = Operation::Repair;
    MockScenario scenario = MockScenario::MillionSlab;
    EfficiencyRepairCase repair_case = EfficiencyRepairCase::Worst;
    GridDimensions dimensions{128, 128, 128};
    double minimum_measured_seconds = 1.0;
    std::uint64_t minimum_repetitions = 10;
    std::uint64_t maximum_repetitions = 1000000;
    std::uint64_t minimum_queries_per_rank = 1000000;
    std::uint64_t warmup_repetitions = 3;
    bool show_help = false;
};

struct Statistics {
    std::uint64_t count = 0;
    double sum = 0.0;
    double sum_squares = 0.0;
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = 0.0;

    void add(double value)
    {
        ++count;
        sum += value;
        sum_squares += value * value;
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }

    double average() const noexcept
    {
        return count == 0 ? 0.0 : sum / static_cast<double>(count);
    }

    double standard_deviation() const noexcept
    {
        if (count == 0) {
            return 0.0;
        }
        const auto mean = average();
        const auto variance = std::max(
            0.0,
            sum_squares / static_cast<double>(count) - mean * mean);
        return std::sqrt(variance);
    }
};

struct RepairMetrics {
    std::uint64_t global_visited = 0;
    std::uint64_t global_closed = 0;
    std::uint64_t global_changed = 0;
    std::uint64_t participating_ranks = 0;
    std::uint64_t closing_ranks = 0;
    std::uint64_t communication_rounds = 0;
    std::uint64_t sent_frontier_entries = 0;
    std::uint64_t received_frontier_entries = 0;
    std::uint64_t minimum_rank_visited = 0;
    std::uint64_t maximum_rank_visited = 0;
    std::uint64_t minimum_rank_closed = 0;
    std::uint64_t maximum_rank_closed = 0;
};

struct InitializationTimes {
    double total = 0.0;
    double grid_construction = 0.0;
    double atom_adaptation_and_voxelization = 0.0;
    double classification = 0.0;
};

void check_mpi(int error_code, const char* operation)
{
    if (error_code != MPI_SUCCESS) {
        throw std::runtime_error(std::string(operation) + " failed");
    }
}

std::uint64_t parse_uint64(const std::string& value, const char* option)
{
    std::size_t parsed_count = 0;
    const auto parsed = std::stoull(value, &parsed_count, 10);
    if (parsed_count != value.size()) {
        throw std::invalid_argument(std::string(option) + " requires an integer");
    }
    return static_cast<std::uint64_t>(parsed);
}

double parse_nonnegative_double(const std::string& value, const char* option)
{
    std::size_t parsed_count = 0;
    const auto parsed = std::stod(value, &parsed_count);
    if (parsed_count != value.size() || parsed < 0.0 || !std::isfinite(parsed)) {
        throw std::invalid_argument(
            std::string(option) + " requires a finite nonnegative value");
    }
    return parsed;
}

Operation parse_operation(const std::string& value)
{
    if (value == "initialization") {
        return Operation::Initialization;
    }
    if (value == "query") {
        return Operation::Query;
    }
    if (value == "repair") {
        return Operation::Repair;
    }
    throw std::invalid_argument(
        "operation must be initialization, query, or repair");
}

const char* operation_name(Operation operation) noexcept
{
    switch (operation) {
    case Operation::Initialization:
        return "initialization";
    case Operation::Query:
        return "query";
    case Operation::Repair:
        return "repair";
    }
    return "unknown";
}

Options parse_options(int argument_count, char** arguments)
{
    Options options{};
    for (int index = 1; index < argument_count; ++index) {
        const std::string argument(arguments[index]);
        const auto require_value = [&]() -> std::string {
            if (index + 1 >= argument_count) {
                throw std::invalid_argument(argument + " requires a value");
            }
            return arguments[++index];
        };
        if (argument == "--help") {
            options.show_help = true;
        } else if (argument == "--operation") {
            options.operation = parse_operation(require_value());
        } else if (argument == "--scenario") {
            options.scenario = gasaccess::testing::parse_mock_scenario(
                require_value());
        } else if (argument == "--case") {
            options.repair_case =
                gasaccess::testing::parse_efficiency_repair_case(require_value());
        } else if (argument == "--nx") {
            options.dimensions.x = parse_uint64(require_value(), "--nx");
        } else if (argument == "--ny") {
            options.dimensions.y = parse_uint64(require_value(), "--ny");
        } else if (argument == "--nz") {
            options.dimensions.z = parse_uint64(require_value(), "--nz");
        } else if (argument == "--min-measured-seconds") {
            options.minimum_measured_seconds = parse_nonnegative_double(
                require_value(), "--min-measured-seconds");
        } else if (argument == "--min-repetitions") {
            options.minimum_repetitions = parse_uint64(
                require_value(), "--min-repetitions");
        } else if (argument == "--max-repetitions") {
            options.maximum_repetitions = parse_uint64(
                require_value(), "--max-repetitions");
        } else if (argument == "--query-count") {
            options.minimum_queries_per_rank = parse_uint64(
                require_value(), "--query-count");
        } else if (argument == "--warmup-repetitions") {
            options.warmup_repetitions = parse_uint64(
                require_value(), "--warmup-repetitions");
        } else {
            throw std::invalid_argument("unknown argument: " + argument);
        }
    }
    if (options.minimum_repetitions == 0
        || options.maximum_repetitions < options.minimum_repetitions) {
        throw std::invalid_argument(
            "repetition bounds must be positive and ordered");
    }
    if (options.minimum_queries_per_rank == 0) {
        throw std::invalid_argument("query count must be positive");
    }
    return options;
}

void print_help()
{
    std::cout
        << "Usage: gasaccess_mpi_efficiency_driver [options]\n"
        << "  --operation initialization|query|repair\n"
        << "  --scenario open-trench|sealed-trench|million-slab\n"
        << "  --case baseline|best|medium|worst\n"
        << "  --nx N --ny N --nz N\n"
        << "  --min-measured-seconds S\n"
        << "  --min-repetitions N --max-repetitions N\n"
        << "  --query-count N\n"
        << "  --warmup-repetitions N\n";
}

double allreduce_max(double local_value)
{
    double global_value = 0.0;
    check_mpi(MPI_Allreduce(
        &local_value,
        &global_value,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD), "MPI_Allreduce(max double)");
    return global_value;
}

std::uint64_t allreduce_sum(std::uint64_t local_value)
{
    std::uint64_t global_value = 0;
    check_mpi(MPI_Allreduce(
        &local_value,
        &global_value,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD), "MPI_Allreduce(sum uint64)");
    return global_value;
}

std::uint64_t allreduce_min(std::uint64_t local_value)
{
    std::uint64_t global_value = 0;
    check_mpi(MPI_Allreduce(
        &local_value,
        &global_value,
        1,
        MPI_UINT64_T,
        MPI_MIN,
        MPI_COMM_WORLD), "MPI_Allreduce(min uint64)");
    return global_value;
}

std::uint64_t allreduce_max(std::uint64_t local_value)
{
    std::uint64_t global_value = 0;
    check_mpi(MPI_Allreduce(
        &local_value,
        &global_value,
        1,
        MPI_UINT64_T,
        MPI_MAX,
        MPI_COMM_WORLD), "MPI_Allreduce(max uint64)");
    return global_value;
}

std::uint64_t peak_rss_bytes() noexcept
{
#if defined(__unix__) || defined(__APPLE__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0 || usage.ru_maxrss < 0) {
        return 0;
    }
#if defined(__APPLE__)
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024U;
#endif
#else
    return 0;
#endif
}

bool global_query(
    const DistributedGasGrid& grid,
    const VoxelCoord& coordinate)
{
    std::uint64_t local_owner = 0;
    std::uint64_t local_accessible = 0;
    if (grid.owns(coordinate)) {
        local_owner = 1;
        local_accessible = DistributedGasAccessibilityQuery(grid)
                .is_site_accessible(grid.voxel_center(coordinate))
            ? 1U
            : 0U;
    }
    if (allreduce_sum(local_owner) != 1) {
        throw std::logic_error("query coordinate does not have one owner");
    }
    return allreduce_sum(local_accessible) != 0;
}

MpiDecompositionSpec make_repair_decomposition(
    const EfficiencyRepairFixture& fixture,
    int rank,
    int process_count)
{
    const auto rank_count = static_cast<std::uint64_t>(process_count);
    const auto local_x = fixture.grid_spec.dimensions.x / rank_count;
    const auto unsigned_rank = static_cast<unsigned int>(rank);
    const auto x_begin = local_x * static_cast<std::uint64_t>(unsigned_rank);
    MpiDecompositionSpec decomposition{};
    decomposition.communicator = MPI_COMM_WORLD;
    decomposition.global_lower = {0.0, 0.0, 0.0};
    decomposition.global_upper = {
        static_cast<double>(fixture.grid_spec.dimensions.x),
        static_cast<double>(fixture.grid_spec.dimensions.y),
        static_cast<double>(fixture.grid_spec.dimensions.z)};
    decomposition.local_lower = {
        static_cast<double>(x_begin), 0.0, 0.0};
    decomposition.local_upper = {
        static_cast<double>(x_begin + local_x),
        decomposition.global_upper.y,
        decomposition.global_upper.z};
    decomposition.process_grid = fixture.process_grid;
    decomposition.process_location = {
        static_cast<std::int64_t>(unsigned_rank), 0, 0};
    decomposition.atom_ghost_distance = 0.0;
    decomposition.maximum_excluded_radius = 0.0;
    return decomposition;
}

std::unique_ptr<DistributedGasGrid> prepare_repair_grid(
    EfficiencyRepairCase repair_case,
    const EfficiencyRepairFixture& fixture,
    int rank,
    int process_count)
{
    auto grid = std::make_unique<DistributedGasGrid>(
        fixture.grid_spec,
        make_repair_decomposition(fixture, rank, process_count));
    const auto& range = grid->owned_range();
    for (auto z = range.begin.z; z < range.end.z; ++z) {
        for (auto y = range.begin.y; y < range.end.y; ++y) {
            for (auto x = range.begin.x; x < range.end.x; ++x) {
                const VoxelCoord coordinate{x, y, z};
                if (gasaccess::testing::is_efficiency_fixture_solid(
                        repair_case, fixture, coordinate)) {
                    grid->set_owned_gas_state(coordinate, GasState::Solid);
                }
            }
        }
    }
    DistributedExteriorClassifier().classify(*grid);
    const auto global_closed = allreduce_sum(
        grid->owned_gas_state_count(GasState::ClosedVoid));
    if (global_closed != fixture.expected_initial_closed_count) {
        throw std::logic_error("repair fixture initial closed count is incorrect");
    }
    return grid;
}

RepairMetrics validate_repair_result(
    EfficiencyRepairCase repair_case,
    const EfficiencyRepairFixture& fixture,
    const DistributedGasGrid& grid,
    const DistributedDepositionUpdateResult& result,
    int process_count)
{
    const bool baseline = repair_case == EfficiencyRepairCase::DetectionBaseline;
    if (baseline) {
        if (result.geometry_changed() || result.used_distributed_repair()) {
            throw std::logic_error("detection baseline unexpectedly changed geometry");
        }
    } else if (!result.geometry_changed()
               || !result.used_distributed_repair()
               || result.used_full_reclassification()) {
        throw std::logic_error("repair benchmark did not use incremental repair");
    }

    RepairMetrics metrics{};
    metrics.global_visited = allreduce_sum(
        result.local_repair_visited_voxel_count);
    metrics.global_closed = allreduce_sum(
        result.local_repair_closed_voxel_count);
    metrics.global_changed = allreduce_sum(
        static_cast<std::uint64_t>(result.changed_owned_voxel_coords.size()));
    metrics.participating_ranks = allreduce_sum(
        result.local_repair_visited_voxel_count == 0 ? 0U : 1U);
    metrics.closing_ranks = allreduce_sum(
        result.local_repair_closed_voxel_count == 0 ? 0U : 1U);
    metrics.communication_rounds = allreduce_max(
        result.repair_communication_round_count);
    metrics.sent_frontier_entries = allreduce_sum(
        result.sent_frontier_entry_count);
    metrics.received_frontier_entries = allreduce_sum(
        result.received_frontier_entry_count);
    metrics.minimum_rank_visited = allreduce_min(
        result.local_repair_visited_voxel_count);
    metrics.maximum_rank_visited = allreduce_max(
        result.local_repair_visited_voxel_count);
    metrics.minimum_rank_closed = allreduce_min(
        result.local_repair_closed_voxel_count);
    metrics.maximum_rank_closed = allreduce_max(
        result.local_repair_closed_voxel_count);

    const auto expected_final_closed = fixture.expected_initial_closed_count
        + fixture.expected_newly_closed_count;
    const auto global_final_closed = allreduce_sum(
        grid.owned_gas_state_count(GasState::ClosedVoid));
    if (global_final_closed != expected_final_closed
        || metrics.global_closed != fixture.expected_newly_closed_count
        || metrics.sent_frontier_entries != metrics.received_frontier_entries) {
        throw std::logic_error("repair benchmark state or communication mismatch");
    }
    const auto expected_changed = baseline
        ? std::uint64_t{0}
        : fixture.expected_newly_closed_count + 1U;
    if (metrics.global_changed != expected_changed) {
        throw std::logic_error("repair benchmark changed-coordinate count is wrong");
    }
    if (!global_query(grid, fixture.outside_probe)) {
        throw std::logic_error("repair benchmark closed the outside control probe");
    }
    if (!baseline && global_query(grid, fixture.cavity_probe)) {
        throw std::logic_error("repair benchmark cavity remained accessible");
    }
    if (repair_case == EfficiencyRepairCase::Best
        && (metrics.participating_ranks != 1 || metrics.closing_ranks != 1)) {
        throw std::logic_error("best repair escaped its owning rank");
    }
    if (repair_case == EfficiencyRepairCase::Worst
        && (metrics.participating_ranks
                != static_cast<std::uint64_t>(process_count)
            || metrics.closing_ranks
                != static_cast<std::uint64_t>(process_count)
            || metrics.minimum_rank_visited == 0
            || metrics.minimum_rank_closed == 0)) {
        throw std::logic_error(
            "worst repair did not visit and close voxels on every MPI rank");
    }
    return metrics;
}

double run_repair_once(
    const Options& options,
    int rank,
    int process_count,
    RepairMetrics* metrics)
{
    const auto fixture = gasaccess::testing::make_efficiency_repair_fixture(
        options.repair_case, options.dimensions, process_count);
    auto grid = prepare_repair_grid(
        options.repair_case, fixture, rank, process_count);
    const auto expected_before =
        options.repair_case != EfficiencyRepairCase::DetectionBaseline;
    if (global_query(*grid, fixture.cavity_probe) != expected_before) {
        throw std::logic_error("repair fixture has wrong initial accessibility");
    }

    const Atom deposited_atom{grid->voxel_center(fixture.opening), 0.0};
    DistributedDepositionUpdater updater(0.0);
    check_mpi(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(repair timing)");
    const double start = MPI_Wtime();
    const auto result = updater.apply_deposition(
        *grid, AtomView{&deposited_atom, 1});
    const double local_elapsed = MPI_Wtime() - start;
    const double elapsed = allreduce_max(local_elapsed);
    const auto actual_metrics = validate_repair_result(
        options.repair_case, fixture, *grid, result, process_count);
    if (metrics != nullptr) {
        *metrics = actual_metrics;
    }
    return elapsed;
}

double run_prepared_detection_baseline_once(
    const EfficiencyRepairFixture& fixture,
    DistributedGasGrid& grid,
    DistributedDepositionUpdater& updater,
    const Atom& deposited_atom,
    int process_count,
    RepairMetrics* metrics)
{
    check_mpi(
        MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(baseline timing)");
    const double start = MPI_Wtime();
    const auto result = updater.apply_deposition(
        grid, AtomView{&deposited_atom, 1});
    const double elapsed = allreduce_max(MPI_Wtime() - start);
    const auto actual_metrics = validate_repair_result(
        EfficiencyRepairCase::DetectionBaseline,
        fixture,
        grid,
        result,
        process_count);
    if (metrics != nullptr) {
        *metrics = actual_metrics;
    }
    return elapsed;
}

void print_statistics(const Statistics& statistics)
{
    std::cout << "repetition_count=" << statistics.count << '\n'
              << "total_measured_seconds=" << statistics.sum << '\n'
              << "average_seconds=" << statistics.average() << '\n'
              << "minimum_seconds="
              << (statistics.count == 0 ? 0.0 : statistics.minimum) << '\n'
              << "maximum_seconds=" << statistics.maximum << '\n'
              << "standard_deviation_seconds="
              << statistics.standard_deviation() << '\n';
}

int run_repair(const Options& options, int rank, int process_count)
{
    Statistics statistics{};
    RepairMetrics metrics{};
    if (options.repair_case == EfficiencyRepairCase::DetectionBaseline) {
        const auto fixture = gasaccess::testing::make_efficiency_repair_fixture(
            options.repair_case, options.dimensions, process_count);
        auto grid = prepare_repair_grid(
            options.repair_case, fixture, rank, process_count);
        const Atom deposited_atom{grid->voxel_center(fixture.opening), 0.0};
        DistributedDepositionUpdater updater(0.0);
        for (std::uint64_t warmup = 0;
             warmup < options.warmup_repetitions;
             ++warmup) {
            (void)run_prepared_detection_baseline_once(
                fixture,
                *grid,
                updater,
                deposited_atom,
                process_count,
                nullptr);
        }
        while ((statistics.count < options.minimum_repetitions
                || statistics.sum < options.minimum_measured_seconds)
               && statistics.count < options.maximum_repetitions) {
            statistics.add(run_prepared_detection_baseline_once(
                fixture,
                *grid,
                updater,
                deposited_atom,
                process_count,
                &metrics));
        }
    } else {
        for (std::uint64_t warmup = 0;
             warmup < options.warmup_repetitions;
             ++warmup) {
            (void)run_repair_once(options, rank, process_count, nullptr);
        }
        while ((statistics.count < options.minimum_repetitions
                || statistics.sum < options.minimum_measured_seconds)
               && statistics.count < options.maximum_repetitions) {
            statistics.add(run_repair_once(
                options, rank, process_count, &metrics));
        }
    }
    const auto global_peak_rss = allreduce_sum(peak_rss_bytes());

    if (rank == 0) {
        const auto fixture = gasaccess::testing::make_efficiency_repair_fixture(
            options.repair_case, options.dimensions, process_count);
        std::cout << "repair_case="
                  << gasaccess::testing::efficiency_repair_case_name(
                         options.repair_case)
                  << '\n'
                  << "expected_newly_closed_voxel_count="
                  << fixture.expected_newly_closed_count << '\n'
                  << "repair_visited_voxel_count=" << metrics.global_visited
                  << '\n'
                  << "repair_closed_voxel_count=" << metrics.global_closed
                  << '\n'
                  << "changed_voxel_count=" << metrics.global_changed << '\n'
                  << "participating_rank_count=" << metrics.participating_ranks
                  << '\n'
                  << "closing_rank_count=" << metrics.closing_ranks << '\n'
                  << "repair_communication_round_count="
                  << metrics.communication_rounds << '\n'
                  << "sent_frontier_entry_count="
                  << metrics.sent_frontier_entries << '\n'
                  << "received_frontier_entry_count="
                  << metrics.received_frontier_entries << '\n'
                  << "minimum_rank_visited_voxel_count="
                  << metrics.minimum_rank_visited << '\n'
                  << "maximum_rank_visited_voxel_count="
                  << metrics.maximum_rank_visited << '\n'
                  << "minimum_rank_closed_voxel_count="
                  << metrics.minimum_rank_closed << '\n'
                  << "maximum_rank_closed_voxel_count="
                  << metrics.maximum_rank_closed << '\n'
                  << "peak_rss_bytes_sum=" << global_peak_rss << '\n';
        print_statistics(statistics);
        std::cout << "average_seconds_per_visited_voxel="
                  << (metrics.global_visited == 0
                          ? 0.0
                          : statistics.average()
                              / static_cast<double>(metrics.global_visited))
                  << '\n';
    }
    return 0;
}

struct MockWorkflow {
    gasaccess::testing::MockFixtureSpec fixture{};
    GridDimensions process_grid{};
    MockSpparksDomain domain;
    MockSpparksApp app;

    MockWorkflow(MockScenario scenario, int process_count)
        : fixture(gasaccess::testing::make_mock_fixture_spec(scenario)),
          process_grid{static_cast<std::uint64_t>(process_count), 1, 1},
          domain(
              MPI_COMM_WORLD,
              fixture.global_lower,
              fixture.global_upper,
              process_grid,
              fixture.spparks_periodic),
          app(domain)
    {
        if (fixture.dimensions.x % process_grid.x != 0) {
            throw std::invalid_argument(
                "mock fixture x dimension must be divisible by rank count");
        }
        app.set_owned_atoms(gasaccess::testing::make_mock_atoms(
            scenario, fixture, &domain));
        app.synchronize_ghost_atoms(fixture.ghost_distance);
    }
};

SpparksAdapterConfig make_adapter_config(const MockWorkflow& workflow)
{
    SpparksAdapterConfig config{};
    config.requested_spacing = workflow.fixture.requested_spacing;
    config.gas_periodic = workflow.fixture.gas_periodic;
    config.reservoir_faces = workflow.fixture.reservoir_faces;
    config.atom_ghost_distance = workflow.fixture.ghost_distance;
    config.maximum_excluded_radius = workflow.fixture.atom_radius
        + workflow.fixture.precursor_radius;
    return config;
}

InitializationTimes run_initialization_once(
    const MockWorkflow& workflow,
    std::uint64_t& global_visited,
    std::uint64_t& global_sent,
    std::uint64_t& global_received)
{
    check_mpi(
        MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(initialization timing)");
    const double start = MPI_Wtime();
    const auto adapter_config = make_adapter_config(workflow);
    const auto grid_spec = gasaccess::make_spparks_grid_spec(
        workflow.domain, adapter_config);
    const auto decomposition_spec = gasaccess::make_spparks_decomposition_spec(
        workflow.domain, workflow.domain.world, adapter_config);
    DistributedGasGrid grid(grid_spec, decomposition_spec);
    const double grid_end = MPI_Wtime();
    SpparksAtomBuffer atom_buffer;
    atom_buffer.assign(
        workflow.app,
        [&workflow](std::size_t atom_index) {
            return workflow.app.radius[atom_index];
        });
    grid.voxelize_owned_atoms(
        atom_buffer.atom_view(), workflow.fixture.precursor_radius);
    const double voxelization_end = MPI_Wtime();
    const auto summary = DistributedExteriorClassifier().classify(grid);
    const double classification_end = MPI_Wtime();
    InitializationTimes times{};
    times.total = allreduce_max(classification_end - start);
    times.grid_construction = allreduce_max(grid_end - start);
    times.atom_adaptation_and_voxelization = allreduce_max(
        voxelization_end - grid_end);
    times.classification = allreduce_max(
        classification_end - voxelization_end);

    const auto global_unclassified = allreduce_sum(
        grid.owned_gas_state_count(GasState::Unclassified));
    const auto global_solid = allreduce_sum(
        grid.owned_gas_state_count(GasState::Solid));
    const auto global_owned_atoms = allreduce_sum(static_cast<std::uint64_t>(
        static_cast<unsigned int>(workflow.app.nlocal)));
    const auto global_outside = allreduce_sum(
        grid.owned_gas_state_count(GasState::OutsideAccessible));
    const auto global_closed = allreduce_sum(
        grid.owned_gas_state_count(GasState::ClosedVoid));
    const auto expected_voxels = grid_spec.dimensions.x
        * grid_spec.dimensions.y * grid_spec.dimensions.z;
    if (global_unclassified != 0 || global_solid != global_owned_atoms
        || global_solid + global_outside + global_closed != expected_voxels) {
        throw std::logic_error("initialization left inconsistent gas states");
    }
    global_visited = allreduce_sum(summary.local_visited_voxel_count);
    global_sent = allreduce_sum(summary.sent_frontier_entry_count);
    global_received = allreduce_sum(summary.received_frontier_entry_count);
    if (global_sent != global_received) {
        throw std::logic_error("initialization frontier counts do not balance");
    }
    return times;
}

int run_initialization(const Options& options, int rank, int process_count)
{
    const MockWorkflow workflow(options.scenario, process_count);
    std::uint64_t visited = 0;
    std::uint64_t sent = 0;
    std::uint64_t received = 0;
    for (std::uint64_t warmup = 0;
         warmup < options.warmup_repetitions;
         ++warmup) {
        (void)run_initialization_once(workflow, visited, sent, received);
    }
    Statistics statistics{};
    Statistics grid_statistics{};
    Statistics voxelization_statistics{};
    Statistics classification_statistics{};
    while ((statistics.count < options.minimum_repetitions
            || statistics.sum < options.minimum_measured_seconds)
           && statistics.count < options.maximum_repetitions) {
        const auto times = run_initialization_once(
            workflow, visited, sent, received);
        statistics.add(times.total);
        grid_statistics.add(times.grid_construction);
        voxelization_statistics.add(times.atom_adaptation_and_voxelization);
        classification_statistics.add(times.classification);
    }
    const auto global_peak_rss = allreduce_sum(peak_rss_bytes());
    if (rank == 0) {
        std::cout << "scenario="
                  << gasaccess::testing::mock_scenario_name(options.scenario)
                  << '\n'
                  << "nx=" << workflow.fixture.dimensions.x << '\n'
                  << "ny=" << workflow.fixture.dimensions.y << '\n'
                  << "nz=" << workflow.fixture.dimensions.z << '\n'
                  << "classification_visited_voxel_count=" << visited << '\n'
                  << "classification_sent_frontier_count=" << sent << '\n'
                  << "classification_received_frontier_count=" << received
                  << '\n'
                  << "peak_rss_bytes_sum=" << global_peak_rss << '\n';
        print_statistics(statistics);
        std::cout << "grid_construction_average_seconds="
                  << grid_statistics.average() << '\n'
                  << "atom_adaptation_and_voxelization_average_seconds="
                  << voxelization_statistics.average() << '\n'
                  << "classification_average_seconds="
                  << classification_statistics.average() << '\n';
    }
    return 0;
}

std::unique_ptr<DistributedGasGrid> initialize_query_grid(
    const MockWorkflow& workflow)
{
    const auto adapter_config = make_adapter_config(workflow);
    auto grid = std::make_unique<DistributedGasGrid>(
        gasaccess::make_spparks_grid_spec(workflow.domain, adapter_config),
        gasaccess::make_spparks_decomposition_spec(
            workflow.domain, workflow.domain.world, adapter_config));
    SpparksAtomBuffer atom_buffer;
    atom_buffer.assign(
        workflow.app,
        [&workflow](std::size_t atom_index) {
            return workflow.app.radius[atom_index];
        });
    grid->voxelize_owned_atoms(
        atom_buffer.atom_view(), workflow.fixture.precursor_radius);
    DistributedExteriorClassifier().classify(*grid);
    return grid;
}

int run_query(const Options& options, int rank, int process_count)
{
    const MockWorkflow workflow(options.scenario, process_count);
    auto grid = initialize_query_grid(workflow);
    std::vector<Point3> coordinates;
    coordinates.reserve(static_cast<std::size_t>(workflow.app.nlocal));
    for (int atom_index = 0; atom_index < workflow.app.nlocal; ++atom_index) {
        coordinates.push_back({
            workflow.app.xyz[atom_index][0],
            workflow.app.xyz[atom_index][1],
            workflow.app.xyz[atom_index][2]});
    }
    if (coordinates.empty()) {
        throw std::logic_error("query benchmark rank has no owned atoms");
    }

    const DistributedGasAccessibilityQuery query(*grid);
    std::uint64_t warmup_accessible = 0;
    for (std::uint64_t warmup = 0;
         warmup < options.warmup_repetitions;
         ++warmup) {
        for (const auto& coordinate : coordinates) {
            if (query.is_site_accessible(coordinate)) {
                ++warmup_accessible;
            }
        }
    }

    check_mpi(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(query timing)");
    std::uint64_t local_query_count = 0;
    std::uint64_t local_accessible_count = 0;
    const double start = MPI_Wtime();
    double local_elapsed = 0.0;
    do {
        for (const auto& coordinate : coordinates) {
            if (query.is_site_accessible(coordinate)) {
                ++local_accessible_count;
            }
        }
        local_query_count += static_cast<std::uint64_t>(coordinates.size());
        local_elapsed = MPI_Wtime() - start;
    } while (local_query_count < options.minimum_queries_per_rank
             || local_elapsed < options.minimum_measured_seconds);

    const auto global_query_count = allreduce_sum(local_query_count);
    const auto global_accessible_count = allreduce_sum(local_accessible_count);
    const auto maximum_elapsed = allreduce_max(local_elapsed);
    const auto local_average = local_elapsed
        / static_cast<double>(local_query_count);
    const auto maximum_rank_average = allreduce_max(local_average);
    const auto global_warmup_accessible = allreduce_sum(warmup_accessible);
    const auto global_peak_rss = allreduce_sum(peak_rss_bytes());
    if (global_accessible_count == 0 && global_warmup_accessible == 0) {
        throw std::logic_error("query benchmark produced no accessible result");
    }
    if (rank == 0) {
        std::cout << "scenario="
                  << gasaccess::testing::mock_scenario_name(options.scenario)
                  << '\n'
                  << "nx=" << workflow.fixture.dimensions.x << '\n'
                  << "ny=" << workflow.fixture.dimensions.y << '\n'
                  << "nz=" << workflow.fixture.dimensions.z << '\n'
                  << "query_count=" << global_query_count << '\n'
                  << "accessible_query_count=" << global_accessible_count
                  << '\n'
                  << "total_measured_seconds=" << maximum_elapsed << '\n'
                  << "average_query_seconds_max_rank="
                  << maximum_rank_average << '\n'
                  << "average_query_nanoseconds_max_rank="
                  << maximum_rank_average * 1.0e9 << '\n'
                  << "global_query_rate_per_second="
                  << static_cast<double>(global_query_count) / maximum_elapsed
                  << '\n'
                  << "peak_rss_bytes_sum=" << global_peak_rss << '\n';
    }
    return 0;
}

int run(const Options& options, int rank, int process_count)
{
    if (rank == 0) {
        std::cout << std::boolalpha << std::fixed << std::setprecision(9)
                  << "driver=gasaccess_mpi_efficiency\n"
                  << "operation=" << operation_name(options.operation) << '\n'
                  << "mpi_ranks=" << process_count << '\n';
        if (options.operation == Operation::Repair) {
            std::cout << "nx=" << options.dimensions.x << '\n'
                      << "ny=" << options.dimensions.y << '\n'
                      << "nz=" << options.dimensions.z << '\n';
        }
        std::cout
                  << "minimum_measured_seconds="
                  << options.minimum_measured_seconds << '\n'
                  << "minimum_repetitions=" << options.minimum_repetitions
                  << '\n'
                  << "warmup_repetitions=" << options.warmup_repetitions
                  << '\n';
    }
    switch (options.operation) {
    case Operation::Initialization:
        return run_initialization(options, rank, process_count);
    case Operation::Query:
        return run_query(options, rank, process_count);
    case Operation::Repair:
        return run_repair(options, rank, process_count);
    }
    throw std::logic_error("unknown efficiency operation");
}

}  // namespace

int main(int argument_count, char** arguments)
{
    if (MPI_Init(&argument_count, &arguments) != MPI_SUCCESS) {
        std::cerr << "MPI_Init failed\n";
        return 1;
    }
    int rank = 0;
    int process_count = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    int status = 0;
    try {
        const auto options = parse_options(argument_count, arguments);
        if (options.show_help) {
            if (rank == 0) {
                print_help();
            }
        } else {
            status = run(options, rank, process_count);
            if (rank == 0) {
                std::cout << "acceptance=pass\n";
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "rank " << rank << ": " << error.what() << '\n';
        status = 1;
    }
    int global_status = 0;
    MPI_Allreduce(
        &status,
        &global_status,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD);
    MPI_Finalize();
    return global_status;
}
