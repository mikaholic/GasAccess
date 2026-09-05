#include "efficiency_fixture.hpp"
#include "gasaccess/distributed_atom_change_updater.hpp"
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
using gasaccess::AtomChangeBatch;
using gasaccess::AtomChangeRepairMode;
using gasaccess::AtomView;
using gasaccess::AccessibilityRepairKind;
using gasaccess::DistributedAtomChangeUpdateResult;
using gasaccess::DistributedAtomChangeUpdater;
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
using gasaccess::testing::EfficiencyChangeKind;
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
    EfficiencyChangeKind change_kind = EfficiencyChangeKind::Deposition;
    EfficiencyRepairCase repair_case = EfficiencyRepairCase::Worst;
    GridDimensions dimensions{128, 128, 128};
    double minimum_measured_seconds = 1.0;
    std::uint64_t minimum_repetitions = 10;
    std::uint64_t maximum_repetitions = 100000000;
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
    std::uint64_t global_blocker_count_changed = 0;
    std::uint64_t global_newly_solid = 0;
    std::uint64_t global_newly_gas = 0;
    std::uint64_t global_changed = 0;
    std::uint64_t closing_visited = 0;
    std::uint64_t opening_visited = 0;
    std::uint64_t closed = 0;
    std::uint64_t opened = 0;
    std::uint64_t closing_participating_ranks = 0;
    std::uint64_t opening_participating_ranks = 0;
    std::uint64_t closing_communication_rounds = 0;
    std::uint64_t opening_communication_rounds = 0;
    std::uint64_t closing_sent_frontier_entries = 0;
    std::uint64_t closing_received_frontier_entries = 0;
    std::uint64_t opening_sent_frontier_entries = 0;
    std::uint64_t opening_received_frontier_entries = 0;
    std::uint64_t minimum_rank_closing_visited = 0;
    std::uint64_t maximum_rank_closing_visited = 0;
    std::uint64_t minimum_rank_opening_visited = 0;
    std::uint64_t maximum_rank_opening_visited = 0;
    std::uint64_t minimum_rank_closed = 0;
    std::uint64_t maximum_rank_closed = 0;
    std::uint64_t minimum_rank_opened = 0;
    std::uint64_t maximum_rank_opened = 0;

    std::uint64_t visited() const noexcept
    {
        return closing_visited + opening_visited;
    }
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
        } else if (argument == "--change-kind") {
            options.change_kind =
                gasaccess::testing::parse_efficiency_change_kind(
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
        << "  --change-kind deposition|desorption|mixed\n"
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

bool voxel_coord_less(
    const VoxelCoord& lhs,
    const VoxelCoord& rhs) noexcept
{
    if (lhs.z != rhs.z) {
        return lhs.z < rhs.z;
    }
    if (lhs.y != rhs.y) {
        return lhs.y < rhs.y;
    }
    return lhs.x < rhs.x;
}

bool in_region(
    const VoxelCoord& coordinate,
    const VoxelCoord& begin,
    const VoxelCoord& end) noexcept
{
    return coordinate.x >= begin.x && coordinate.x < end.x
        && coordinate.y >= begin.y && coordinate.y < end.y
        && coordinate.z >= begin.z && coordinate.z < end.z;
}

GasState initial_repair_state(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case,
    const EfficiencyRepairFixture& fixture,
    const VoxelCoord& coordinate)
{
    if (gasaccess::testing::is_efficiency_fixture_solid(
            change_kind, repair_case, fixture, coordinate)) {
        return GasState::Solid;
    }
    const bool single_closed_cavity =
        change_kind == EfficiencyChangeKind::Desorption
        || (change_kind == EfficiencyChangeKind::Deposition
            && repair_case == EfficiencyRepairCase::DetectionBaseline);
    if (single_closed_cavity
        && in_region(
            coordinate,
            fixture.opening_cavity_begin,
            fixture.opening_cavity_end)) {
        return GasState::ClosedVoid;
    }
    if (change_kind == EfficiencyChangeKind::Mixed
        && repair_case != EfficiencyRepairCase::DetectionBaseline
        && in_region(
            coordinate,
            fixture.opening_cavity_begin,
            fixture.opening_cavity_end)) {
        return GasState::ClosedVoid;
    }
    return GasState::OutsideAccessible;
}

gasaccess::VoxelBlockerCount initial_blocker_count(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case,
    const EfficiencyRepairFixture& fixture,
    const VoxelCoord& coordinate)
{
    if (fixture.has_removed_atom && coordinate == fixture.removed_voxel) {
        return fixture.initial_removed_blocker_count;
    }
    if (fixture.has_added_atom && coordinate == fixture.added_voxel) {
        return fixture.initial_added_blocker_count;
    }
    return gasaccess::testing::is_efficiency_fixture_solid(
               change_kind, repair_case, fixture, coordinate)
        ? gasaccess::VoxelBlockerCount{1}
        : gasaccess::VoxelBlockerCount{0};
}

std::unique_ptr<DistributedGasGrid> prepare_repair_grid(
    EfficiencyChangeKind change_kind,
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
                        change_kind, repair_case, fixture, coordinate)) {
                    grid->set_owned_gas_state(coordinate, GasState::Solid);
                }
            }
        }
    }
    if (fixture.has_added_atom && grid->owns(fixture.added_voxel)) {
        grid->set_owned_blocker_count(
            fixture.added_voxel,
            fixture.initial_added_blocker_count);
    }
    if (fixture.has_removed_atom && grid->owns(fixture.removed_voxel)) {
        grid->set_owned_blocker_count(
            fixture.removed_voxel,
            fixture.initial_removed_blocker_count);
    }
    DistributedExteriorClassifier().classify(*grid);
    const auto global_closed = allreduce_sum(
        grid->owned_gas_state_count(GasState::ClosedVoid));
    if (global_closed != fixture.expected_initial_closed_count) {
        throw std::logic_error("repair fixture initial closed count is incorrect");
    }
    return grid;
}

void require_collectively(bool local_condition, const char* message)
{
    if (allreduce_sum(local_condition ? 0U : 1U) != 0) {
        throw std::logic_error(message);
    }
}

bool is_sorted_unique(const std::vector<VoxelCoord>& coordinates)
{
    return std::is_sorted(
               coordinates.begin(), coordinates.end(), voxel_coord_less)
        && std::adjacent_find(coordinates.begin(), coordinates.end())
            == coordinates.end();
}

std::uint64_t expected_blocker_changes(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case) noexcept
{
    if (change_kind == EfficiencyChangeKind::Mixed) {
        return repair_case == EfficiencyRepairCase::DetectionBaseline ? 0U : 2U;
    }
    return 1U;
}

AccessibilityRepairKind expected_incremental_repair_kind(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case) noexcept
{
    if (repair_case == EfficiencyRepairCase::DetectionBaseline) {
        return AccessibilityRepairKind::None;
    }
    switch (change_kind) {
    case EfficiencyChangeKind::Deposition:
        return AccessibilityRepairKind::Closing;
    case EfficiencyChangeKind::Desorption:
        return AccessibilityRepairKind::Opening;
    case EfficiencyChangeKind::Mixed:
        return AccessibilityRepairKind::Mixed;
    }
    return AccessibilityRepairKind::None;
}

RepairMetrics validate_repair_result(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case,
    const EfficiencyRepairFixture& fixture,
    const DistributedGasGrid& grid,
    const DistributedAtomChangeUpdateResult& result,
    int process_count,
    bool incremental)
{
    const bool baseline = repair_case == EfficiencyRepairCase::DetectionBaseline;
    if (baseline) {
        if (result.geometry_changed()
            || result.used_full_reclassification()
            || result.used_distributed_closing_repair()
            || result.used_distributed_opening_repair()
            || result.repair_kind != AccessibilityRepairKind::None) {
            throw std::logic_error("detection baseline unexpectedly changed geometry");
        }
    } else if (incremental) {
        if (!result.geometry_changed()
            || result.used_full_reclassification()
            || result.repair_kind
                != expected_incremental_repair_kind(change_kind, repair_case)
            || result.used_distributed_closing_repair()
                != (fixture.expected_newly_solid_count != 0)
            || result.used_distributed_opening_repair()
                != (fixture.expected_newly_gas_count != 0)) {
            throw std::logic_error(
                "repair benchmark did not use the expected incremental passes");
        }
    } else if (!result.geometry_changed()
               || !result.used_full_reclassification()
               || result.repair_kind
                   != AccessibilityRepairKind::FullReclassification
               || result.used_distributed_closing_repair()
               || result.used_distributed_opening_repair()) {
        throw std::logic_error(
            "repair benchmark did not use forced full reclassification");
    }

    if (result.global_blocker_count_changed_voxel_count
            != expected_blocker_changes(change_kind, repair_case)
        || result.global_newly_solid_count
            != fixture.expected_newly_solid_count
        || result.global_newly_gas_count != fixture.expected_newly_gas_count) {
        throw std::logic_error("repair benchmark occupancy counts are incorrect");
    }

    RepairMetrics metrics{};
    metrics.global_blocker_count_changed =
        result.global_blocker_count_changed_voxel_count;
    metrics.global_newly_solid = result.global_newly_solid_count;
    metrics.global_newly_gas = result.global_newly_gas_count;
    metrics.global_changed = allreduce_sum(
        static_cast<std::uint64_t>(result.changed_owned_voxel_coords.size()));
    metrics.closing_visited = allreduce_sum(
        result.local_closing_visited_voxel_count);
    metrics.opening_visited = allreduce_sum(
        result.local_opening_visited_voxel_count);
    metrics.closed = allreduce_sum(
        result.local_repair_closed_voxel_count);
    metrics.opened = allreduce_sum(
        result.local_repair_opened_voxel_count);
    metrics.closing_participating_ranks = allreduce_max(
        result.closing_participating_rank_count);
    metrics.opening_participating_ranks = allreduce_max(
        result.opening_participating_rank_count);
    metrics.closing_communication_rounds = allreduce_max(
        result.closing_communication_round_count);
    metrics.opening_communication_rounds = allreduce_max(
        result.opening_communication_round_count);
    metrics.closing_sent_frontier_entries = allreduce_sum(
        result.closing_sent_frontier_entry_count);
    metrics.closing_received_frontier_entries = allreduce_sum(
        result.closing_received_frontier_entry_count);
    metrics.opening_sent_frontier_entries = allreduce_sum(
        result.opening_sent_frontier_entry_count);
    metrics.opening_received_frontier_entries = allreduce_sum(
        result.opening_received_frontier_entry_count);
    metrics.minimum_rank_closing_visited = allreduce_min(
        result.local_closing_visited_voxel_count);
    metrics.maximum_rank_closing_visited = allreduce_max(
        result.local_closing_visited_voxel_count);
    metrics.minimum_rank_opening_visited = allreduce_min(
        result.local_opening_visited_voxel_count);
    metrics.maximum_rank_opening_visited = allreduce_max(
        result.local_opening_visited_voxel_count);
    metrics.minimum_rank_closed = allreduce_min(
        result.local_repair_closed_voxel_count);
    metrics.maximum_rank_closed = allreduce_max(
        result.local_repair_closed_voxel_count);
    metrics.minimum_rank_opened = allreduce_min(
        result.local_repair_opened_voxel_count);
    metrics.maximum_rank_opened = allreduce_max(
        result.local_repair_opened_voxel_count);

    const auto global_final_closed = allreduce_sum(
        grid.owned_gas_state_count(GasState::ClosedVoid));
    if (global_final_closed != fixture.expected_final_closed_count
        || metrics.global_changed != fixture.expected_changed_count) {
        throw std::logic_error("repair benchmark state or communication mismatch");
    }
    require_collectively(
        is_sorted_unique(result.changed_owned_voxel_coords),
        "repair benchmark changed coordinates are not sorted and unique");

    if (incremental
        && (metrics.closed != fixture.expected_newly_closed_count
            || metrics.opened != fixture.expected_newly_opened_count
            || metrics.closing_sent_frontier_entries
                != metrics.closing_received_frontier_entries
            || metrics.opening_sent_frontier_entries
                != metrics.opening_received_frontier_entries)) {
        throw std::logic_error("repair benchmark traversal metrics are incorrect");
    }
    if (!global_query(grid, fixture.outside_probe)) {
        throw std::logic_error("repair benchmark closed the outside control probe");
    }
    if (!baseline && fixture.expected_newly_closed_count != 0
        && global_query(grid, fixture.closing_probe)) {
        throw std::logic_error("repair benchmark closing cavity remained accessible");
    }
    if (!baseline && fixture.expected_newly_opened_count != 0
        && !global_query(grid, fixture.opening_probe)) {
        throw std::logic_error("repair benchmark opening cavity remained closed");
    }

    if (!incremental) {
        return metrics;
    }

    if (repair_case == EfficiencyRepairCase::Best) {
        if (fixture.expected_newly_closed_count != 0
            && metrics.closing_participating_ranks != 1) {
            throw std::logic_error("best closing repair escaped its owning rank");
        }
        if (fixture.expected_newly_opened_count != 0
            && metrics.opening_participating_ranks != 1) {
            throw std::logic_error("best opening repair escaped its owning rank");
        }
    }
    if (repair_case == EfficiencyRepairCase::Worst) {
        const auto rank_count = static_cast<std::uint64_t>(process_count);
        if (fixture.expected_newly_closed_count != 0
            && (metrics.closing_participating_ranks != rank_count
                || metrics.minimum_rank_closing_visited == 0
                || metrics.minimum_rank_closed == 0)) {
            throw std::logic_error(
                "worst closing repair did not traverse every MPI rank");
        }
        if (fixture.expected_newly_opened_count != 0
            && (metrics.opening_participating_ranks != rank_count
                || metrics.minimum_rank_opening_visited == 0
                || metrics.minimum_rank_opened == 0)) {
            throw std::logic_error(
                "worst opening repair did not traverse every MPI rank");
        }
        const auto& dimensions = fixture.grid_spec.dimensions;
        const auto total_voxels = dimensions.x * dimensions.y * dimensions.z;
        if (change_kind == EfficiencyChangeKind::Desorption
            && fixture.expected_newly_opened_count <= total_voxels / 2U) {
            throw std::logic_error(
                "worst desorption fixture does not open more than half the grid");
        }
        if (change_kind == EfficiencyChangeKind::Mixed
            && (fixture.expected_newly_closed_count
                    + fixture.expected_newly_opened_count
                <= total_voxels / 2U)) {
            throw std::logic_error(
                "worst mixed fixture does not affect more than half the grid");
        }
    }
    return metrics;
}

void compare_repair_results(
    const DistributedGasGrid& incremental_grid,
    const DistributedGasGrid& full_grid,
    const DistributedAtomChangeUpdateResult& incremental_result,
    const DistributedAtomChangeUpdateResult& full_result)
{
    require_collectively(
        incremental_result.changed_owned_voxel_coords
            == full_result.changed_owned_voxel_coords,
        "incremental and full changed-coordinate lists differ");
    std::uint64_t local_mismatch_count = 0;
    const auto& range = incremental_grid.owned_range();
    for (auto z = range.begin.z; z < range.end.z; ++z) {
        for (auto y = range.begin.y; y < range.end.y; ++y) {
            for (auto x = range.begin.x; x < range.end.x; ++x) {
                const VoxelCoord coordinate{x, y, z};
                if (incremental_grid.gas_state(coordinate)
                        != full_grid.gas_state(coordinate)
                    || incremental_grid.owned_blocker_count(coordinate)
                        != full_grid.owned_blocker_count(coordinate)) {
                    ++local_mismatch_count;
                }
            }
        }
    }
    if (allreduce_sum(local_mismatch_count) != 0) {
        throw std::logic_error(
            "incremental and full repair grids do not match");
    }
}

void restore_repair_grid(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case,
    const EfficiencyRepairFixture& fixture,
    DistributedGasGrid& grid,
    const DistributedAtomChangeUpdateResult& result)
{
    const auto restore_coordinate = [&](const VoxelCoord& coordinate) {
        if (!grid.owns(coordinate)) {
            return;
        }
        grid.set_owned_blocker_count(
            coordinate,
            initial_blocker_count(
                change_kind, repair_case, fixture, coordinate));
        grid.set_owned_gas_state(
            coordinate,
            initial_repair_state(
                change_kind, repair_case, fixture, coordinate));
    };
    for (const auto& coordinate : result.changed_owned_voxel_coords) {
        restore_coordinate(coordinate);
    }
    if (fixture.has_added_atom) {
        restore_coordinate(fixture.added_voxel);
    }
    if (fixture.has_removed_atom) {
        restore_coordinate(fixture.removed_voxel);
    }
    if (result.geometry_changed()) {
        grid.exchange_ghost_states();
    }
}

struct TimedRepairResult {
    double elapsed_seconds = 0.0;
    DistributedAtomChangeUpdateResult update{};
};

TimedRepairResult time_atom_changes(
    DistributedGasGrid& grid,
    DistributedAtomChangeUpdater& updater,
    const AtomChangeBatch& changes)
{
    check_mpi(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(repair timing)");
    const double start = MPI_Wtime();
    auto result = updater.apply_atom_changes(grid, changes);
    const auto elapsed = allreduce_max(MPI_Wtime() - start);
    return {elapsed, std::move(result)};
}

struct RepairPairSample {
    double incremental_seconds = 0.0;
    double full_seconds = 0.0;
};

RepairPairSample run_repair_pair(
    const Options& options,
    const EfficiencyRepairFixture& fixture,
    DistributedGasGrid& incremental_grid,
    DistributedGasGrid& full_grid,
    DistributedAtomChangeUpdater& incremental_updater,
    DistributedAtomChangeUpdater& full_updater,
    const AtomChangeBatch& changes,
    int process_count,
    bool full_first,
    bool validate,
    RepairMetrics* metrics)
{
    TimedRepairResult incremental{};
    TimedRepairResult full{};
    if (full_first) {
        full = time_atom_changes(full_grid, full_updater, changes);
        incremental = time_atom_changes(
            incremental_grid, incremental_updater, changes);
    } else {
        incremental = time_atom_changes(
            incremental_grid, incremental_updater, changes);
        full = time_atom_changes(full_grid, full_updater, changes);
    }

    if (validate) {
        const auto actual_metrics = validate_repair_result(
            options.change_kind,
            options.repair_case,
            fixture,
            incremental_grid,
            incremental.update,
            process_count,
            true);
        (void)validate_repair_result(
            options.change_kind,
            options.repair_case,
            fixture,
            full_grid,
            full.update,
            process_count,
            false);
        compare_repair_results(
            incremental_grid,
            full_grid,
            incremental.update,
            full.update);
        if (metrics != nullptr) {
            *metrics = actual_metrics;
        }
    }

    restore_repair_grid(
        options.change_kind,
        options.repair_case,
        fixture,
        incremental_grid,
        incremental.update);
    restore_repair_grid(
        options.change_kind,
        options.repair_case,
        fixture,
        full_grid,
        full.update);
    return {incremental.elapsed_seconds, full.elapsed_seconds};
}

double run_single_repair_measurement(
    const Options& options,
    const EfficiencyRepairFixture& fixture,
    DistributedGasGrid& grid,
    DistributedAtomChangeUpdater& updater,
    const AtomChangeBatch& changes)
{
    auto timed = time_atom_changes(grid, updater, changes);
    restore_repair_grid(
        options.change_kind,
        options.repair_case,
        fixture,
        grid,
        timed.update);
    return timed.elapsed_seconds;
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

void print_prefixed_statistics(
    const char* prefix,
    const Statistics& statistics)
{
    std::cout << prefix << "repetition_count=" << statistics.count << '\n'
              << prefix << "total_measured_seconds=" << statistics.sum << '\n'
              << prefix << "average_seconds=" << statistics.average() << '\n'
              << prefix << "minimum_seconds="
              << (statistics.count == 0 ? 0.0 : statistics.minimum) << '\n'
              << prefix << "maximum_seconds=" << statistics.maximum << '\n'
              << prefix << "standard_deviation_seconds="
              << statistics.standard_deviation() << '\n';
}

int run_repair(const Options& options, int rank, int process_count)
{
    const auto fixture = gasaccess::testing::make_efficiency_repair_fixture(
        options.change_kind,
        options.repair_case,
        options.dimensions,
        process_count);
    auto incremental_grid = prepare_repair_grid(
        options.change_kind,
        options.repair_case,
        fixture,
        rank,
        process_count);
    auto full_grid = prepare_repair_grid(
        options.change_kind,
        options.repair_case,
        fixture,
        rank,
        process_count);

    if (!global_query(*incremental_grid, fixture.outside_probe)) {
        throw std::logic_error("repair fixture outside probe is not accessible");
    }
    if (options.repair_case != EfficiencyRepairCase::DetectionBaseline
        && fixture.expected_newly_closed_count != 0
        && !global_query(*incremental_grid, fixture.closing_probe)) {
        throw std::logic_error("repair fixture closing probe is not accessible");
    }
    if (options.repair_case != EfficiencyRepairCase::DetectionBaseline
        && fixture.expected_newly_opened_count != 0
        && global_query(*incremental_grid, fixture.opening_probe)) {
        throw std::logic_error("repair fixture opening probe is not closed");
    }

    std::vector<Atom> added_atoms;
    std::vector<Atom> removed_atoms;
    if (fixture.has_added_atom) {
        added_atoms.push_back({
            incremental_grid->voxel_center(fixture.added_voxel), 0.0});
    }
    if (fixture.has_removed_atom) {
        removed_atoms.push_back({
            incremental_grid->voxel_center(fixture.removed_voxel), 0.0});
    }
    const AtomChangeBatch changes{
        {added_atoms.data(), added_atoms.size()},
        {removed_atoms.data(), removed_atoms.size()}};
    DistributedAtomChangeUpdater incremental_updater(0.0);
    DistributedAtomChangeUpdater full_updater(
        0.0, AtomChangeRepairMode::FullReclassification);

    std::uint64_t sequence = 0;
    RepairMetrics metrics{};
    bool validated = false;
    for (std::uint64_t warmup = 0;
         warmup < options.warmup_repetitions;
         ++warmup, ++sequence) {
        (void)run_repair_pair(
            options,
            fixture,
            *incremental_grid,
            *full_grid,
            incremental_updater,
            full_updater,
            changes,
            process_count,
            sequence % 2U != 0,
            !validated,
            &metrics);
        validated = true;
    }

    Statistics incremental_statistics{};
    Statistics full_statistics{};
    if (!validated) {
        const auto sample = run_repair_pair(
            options,
            fixture,
            *incremental_grid,
            *full_grid,
            incremental_updater,
            full_updater,
            changes,
            process_count,
            sequence % 2U != 0,
            true,
            &metrics);
        ++sequence;
        incremental_statistics.add(sample.incremental_seconds);
        full_statistics.add(sample.full_seconds);
        validated = true;
    }

    const auto needs_more = [&](const Statistics& statistics) {
        return statistics.count < options.minimum_repetitions
            || statistics.sum < options.minimum_measured_seconds;
    };
    while ((needs_more(incremental_statistics)
            || needs_more(full_statistics))
           && (incremental_statistics.count < options.maximum_repetitions
               || full_statistics.count < options.maximum_repetitions)) {
        const bool run_incremental = needs_more(incremental_statistics)
            && incremental_statistics.count < options.maximum_repetitions;
        const bool run_full = needs_more(full_statistics)
            && full_statistics.count < options.maximum_repetitions;
        if (!run_incremental && !run_full) {
            break;
        }
        const auto measure_incremental = [&]() {
            incremental_statistics.add(run_single_repair_measurement(
                options,
                fixture,
                *incremental_grid,
                incremental_updater,
                changes));
        };
        const auto measure_full = [&]() {
            full_statistics.add(run_single_repair_measurement(
                options,
                fixture,
                *full_grid,
                full_updater,
                changes));
        };
        if (run_incremental && run_full && sequence % 2U != 0) {
            measure_full();
            measure_incremental();
        } else {
            if (run_incremental) {
                measure_incremental();
            }
            if (run_full) {
                measure_full();
            }
        }
        ++sequence;
    }
    if (incremental_statistics.count < options.minimum_repetitions
        || incremental_statistics.sum < options.minimum_measured_seconds
        || full_statistics.sum < options.minimum_measured_seconds) {
        throw std::runtime_error(
            "maximum repetitions reached before timing requirements");
    }
    const auto global_peak_rss = allreduce_sum(peak_rss_bytes());

    if (rank == 0) {
        std::cout << "change_kind="
                  << gasaccess::testing::efficiency_change_kind_name(
                         options.change_kind)
                  << '\n'
                  << "repair_case="
                  << gasaccess::testing::efficiency_repair_case_name(
                         options.repair_case)
                  << '\n'
                  << "expected_blocker_count_changed_voxel_count="
                  << expected_blocker_changes(
                         options.change_kind, options.repair_case)
                  << '\n'
                  << "expected_newly_solid_voxel_count="
                  << fixture.expected_newly_solid_count << '\n'
                  << "expected_newly_gas_voxel_count="
                  << fixture.expected_newly_gas_count << '\n'
                  << "expected_newly_closed_voxel_count="
                  << fixture.expected_newly_closed_count << '\n'
                  << "expected_newly_opened_voxel_count="
                  << fixture.expected_newly_opened_count << '\n'
                  << "blocker_count_changed_voxel_count="
                  << metrics.global_blocker_count_changed << '\n'
                  << "newly_solid_voxel_count=" << metrics.global_newly_solid
                  << '\n'
                  << "newly_gas_voxel_count=" << metrics.global_newly_gas
                  << '\n'
                  << "repair_visited_voxel_count=" << metrics.visited()
                  << '\n'
                  << "closing_visited_voxel_count=" << metrics.closing_visited
                  << '\n'
                  << "opening_visited_voxel_count=" << metrics.opening_visited
                  << '\n'
                  << "repair_closed_voxel_count=" << metrics.closed << '\n'
                  << "repair_opened_voxel_count=" << metrics.opened << '\n'
                  << "changed_voxel_count=" << metrics.global_changed << '\n'
                  << "closing_participating_rank_count="
                  << metrics.closing_participating_ranks
                  << '\n'
                  << "opening_participating_rank_count="
                  << metrics.opening_participating_ranks << '\n'
                  << "closing_communication_round_count="
                  << metrics.closing_communication_rounds << '\n'
                  << "opening_communication_round_count="
                  << metrics.opening_communication_rounds << '\n'
                  << "closing_sent_frontier_entry_count="
                  << metrics.closing_sent_frontier_entries << '\n'
                  << "closing_received_frontier_entry_count="
                  << metrics.closing_received_frontier_entries << '\n'
                  << "opening_sent_frontier_entry_count="
                  << metrics.opening_sent_frontier_entries << '\n'
                  << "opening_received_frontier_entry_count="
                  << metrics.opening_received_frontier_entries << '\n'
                  << "minimum_rank_closing_visited_voxel_count="
                  << metrics.minimum_rank_closing_visited << '\n'
                  << "maximum_rank_closing_visited_voxel_count="
                  << metrics.maximum_rank_closing_visited << '\n'
                  << "minimum_rank_opening_visited_voxel_count="
                  << metrics.minimum_rank_opening_visited << '\n'
                  << "maximum_rank_opening_visited_voxel_count="
                  << metrics.maximum_rank_opening_visited << '\n'
                  << "minimum_rank_closed_voxel_count="
                  << metrics.minimum_rank_closed << '\n'
                  << "maximum_rank_closed_voxel_count="
                  << metrics.maximum_rank_closed << '\n'
                  << "minimum_rank_opened_voxel_count="
                  << metrics.minimum_rank_opened << '\n'
                  << "maximum_rank_opened_voxel_count="
                  << metrics.maximum_rank_opened << '\n'
                  << "phase_timing_enabled=false\n"
                  << "peak_rss_bytes_sum=" << global_peak_rss << '\n';
        print_statistics(incremental_statistics);
        print_prefixed_statistics("incremental_", incremental_statistics);
        print_prefixed_statistics("full_reclassification_", full_statistics);
        std::cout << "full_to_incremental_speedup="
                  << (incremental_statistics.average() == 0.0
                          ? 0.0
                          : full_statistics.average()
                              / incremental_statistics.average())
                  << '\n';
        std::cout << "average_seconds_per_visited_voxel="
                  << (metrics.visited() == 0
                          ? 0.0
                          : incremental_statistics.average()
                              / static_cast<double>(metrics.visited()))
                  << '\n'
                  << "full_reclassification_average_seconds_per_grid_voxel="
                  << full_statistics.average()
                      / static_cast<double>(
                          options.dimensions.x
                          * options.dimensions.y
                          * options.dimensions.z)
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
