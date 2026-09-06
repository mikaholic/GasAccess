#include "gasaccess/distributed_exterior_classifier.hpp"
#include "gasaccess/spparks_adapter.hpp"
#include "mock_spparks.hpp"
#include "mock_spparks_fixture.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

namespace {

using gasaccess::DistributedExteriorClassifier;
using gasaccess::DistributedGasAccessibilityQuery;
using gasaccess::DistributedGasGrid;
using gasaccess::GasState;
using gasaccess::Point3;
using gasaccess::SpparksAdapterConfig;
using gasaccess::SpparksAtomBuffer;
using gasaccess::testing::MockScenario;
using gasaccess::testing::MockSpparksApp;
using gasaccess::testing::MockSpparksDomain;

struct Options {
    MockScenario scenario = MockScenario::OpenTrench;
    bool show_help = false;
};

Options parse_options(int argument_count, char** arguments)
{
    Options options{};
    for (int index = 1; index < argument_count; ++index) {
        const std::string argument(arguments[index]);
        if (argument == "--help") {
            options.show_help = true;
        } else if (argument == "--scenario") {
            if (index + 1 >= argument_count) {
                throw std::invalid_argument("--scenario requires a value");
            }
            options.scenario = gasaccess::testing::parse_mock_scenario(
                arguments[++index]);
        } else {
            throw std::invalid_argument("unknown argument: " + argument);
        }
    }
    return options;
}

void print_help()
{
    std::cout
        << "Usage: gasaccess_spparks_mock_driver [options]\n"
        << "  --scenario open-trench|sealed-trench|million-slab\n"
        << "  --help\n";
}

double reduce_max(double local_value)
{
    double global_value = 0.0;
    if (MPI_Reduce(
            &local_value,
            &global_value,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0,
            MPI_COMM_WORLD)
        != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Reduce failed");
    }
    return global_value;
}

std::uint64_t reduce_sum(std::uint64_t local_value)
{
    std::uint64_t global_value = 0;
    if (MPI_Reduce(
            &local_value,
            &global_value,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            0,
            MPI_COMM_WORLD)
        != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Reduce failed");
    }
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

int run(const Options& options, int rank, int process_count)
{
    const auto fixture_spec =
        gasaccess::testing::make_mock_fixture_spec(options.scenario);
    const auto process_grid =
        gasaccess::testing::make_mock_process_grid(process_count);
    const MockSpparksDomain domain(
        MPI_COMM_WORLD,
        fixture_spec.global_lower,
        fixture_spec.global_upper,
        process_grid,
        fixture_spec.spparks_periodic);

    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Barrier failed");
    }
    const double total_start = MPI_Wtime();

    const double generation_start = MPI_Wtime();
    auto owned_atoms = gasaccess::testing::make_mock_atoms(
        options.scenario, fixture_spec, &domain);
    MockSpparksApp app(domain);
    app.set_owned_atoms(std::move(owned_atoms));
    const double generation_seconds = MPI_Wtime() - generation_start;

    const double ghost_start = MPI_Wtime();
    app.synchronize_ghost_atoms(fixture_spec.ghost_distance);
    const double ghost_seconds = MPI_Wtime() - ghost_start;

    SpparksAdapterConfig adapter_config{};
    adapter_config.requested_spacing = fixture_spec.requested_spacing;
    adapter_config.gas_periodic = fixture_spec.gas_periodic;
    adapter_config.reservoir_faces = fixture_spec.reservoir_faces;
    adapter_config.atom_ghost_distance = fixture_spec.ghost_distance;
    adapter_config.maximum_excluded_radius =
        fixture_spec.atom_radius + fixture_spec.precursor_radius;

    const double grid_start = MPI_Wtime();
    const auto grid_spec = gasaccess::make_spparks_grid_spec(
        domain, adapter_config);
    const auto decomposition_spec = gasaccess::make_spparks_decomposition_spec(
        domain, domain.world, adapter_config);
    DistributedGasGrid gas_grid(grid_spec, decomposition_spec);
    const double grid_seconds = MPI_Wtime() - grid_start;

    const double voxelization_start = MPI_Wtime();
    SpparksAtomBuffer atom_buffer;
    atom_buffer.assign(
        app,
        [&app](std::size_t atom_index) {
            return app.radius[atom_index];
        });
    const auto local_newly_solid = gas_grid.voxelize_owned_atoms(
        atom_buffer.atom_view(), fixture_spec.precursor_radius);
    const double voxelization_seconds = MPI_Wtime() - voxelization_start;

    const double classification_start = MPI_Wtime();
    DistributedExteriorClassifier classifier;
    const auto classification_summary = classifier.classify(gas_grid);
    const double classification_seconds = MPI_Wtime() - classification_start;

    const double query_start = MPI_Wtime();
    const DistributedGasAccessibilityQuery query(gas_grid);
    std::uint64_t local_accessible_atom_count = 0;
    std::uint64_t local_inner_probe_count = 0;
    std::uint64_t local_accessible_inner_probe_count = 0;
    for (int atom_index = 0; atom_index < app.nlocal; ++atom_index) {
        const Point3 atom_position{
            app.xyz[atom_index][0],
            app.xyz[atom_index][1],
            app.xyz[atom_index][2]};
        const bool accessible = query.is_site_accessible(atom_position);
        if (accessible) {
            ++local_accessible_atom_count;
        }
        const auto coordinate = gas_grid.locate_voxel(atom_position);
        if (!coordinate.has_value()) {
            throw std::logic_error("owned mock atom is outside the gas grid");
        }
        if (gasaccess::testing::is_inner_trench_probe(
                options.scenario, coordinate.value())) {
            ++local_inner_probe_count;
            if (accessible) {
                ++local_accessible_inner_probe_count;
            }
        }
    }
    const double query_seconds = MPI_Wtime() - query_start;
    const double total_seconds = MPI_Wtime() - total_start;

    const auto global_owned_atom_count = reduce_sum(
        static_cast<std::uint64_t>(app.nlocal));
    const auto global_ghost_atom_count = reduce_sum(
        static_cast<std::uint64_t>(app.nghost));
    const auto global_solid_count = reduce_sum(
        gas_grid.owned_gas_state_count(GasState::Solid));
    const auto global_outside_count = reduce_sum(
        gas_grid.owned_gas_state_count(GasState::OutsideAccessible));
    const auto global_closed_count = reduce_sum(
        gas_grid.owned_gas_state_count(GasState::ClosedVoid));
    const auto global_newly_solid_count = reduce_sum(local_newly_solid);
    const auto global_accessible_atom_count = reduce_sum(
        local_accessible_atom_count);
    const auto global_inner_probe_count = reduce_sum(local_inner_probe_count);
    const auto global_accessible_inner_probe_count = reduce_sum(
        local_accessible_inner_probe_count);
    const auto global_visited_count = reduce_sum(
        classification_summary.local_visited_voxel_count);
    const auto global_sent_frontier_count = reduce_sum(
        classification_summary.sent_frontier_entry_count);
    const auto global_received_frontier_count = reduce_sum(
        classification_summary.received_frontier_entry_count);
    const auto global_peak_rss = reduce_sum(peak_rss_bytes());

    const auto max_generation_seconds = reduce_max(generation_seconds);
    const auto max_ghost_seconds = reduce_max(ghost_seconds);
    const auto max_grid_seconds = reduce_max(grid_seconds);
    const auto max_voxelization_seconds = reduce_max(voxelization_seconds);
    const auto max_classification_seconds = reduce_max(classification_seconds);
    const auto max_query_seconds = reduce_max(query_seconds);
    const auto max_total_seconds = reduce_max(total_seconds);

    if (rank != 0) {
        return 0;
    }

    const auto expected_voxel_count = fixture_spec.dimensions.x
        * fixture_spec.dimensions.y * fixture_spec.dimensions.z;
    if (global_solid_count != global_owned_atom_count
        || global_newly_solid_count != global_owned_atom_count
        || global_solid_count + global_outside_count + global_closed_count
            != expected_voxel_count) {
        throw std::runtime_error("mock integration state counts are inconsistent");
    }
    if (options.scenario == MockScenario::MillionSlab) {
        const auto expected_atoms = std::uint64_t{128} * 128U * 64U;
        const auto expected_accessible = std::uint64_t{128} * 128U;
        if (global_owned_atom_count != expected_atoms
            || global_accessible_atom_count != expected_accessible) {
            throw std::runtime_error("million-slab acceptance result is incorrect");
        }
    } else {
        const auto expected_accessible_inner =
            options.scenario == MockScenario::OpenTrench
            ? global_inner_probe_count
            : std::uint64_t{0};
        if (global_inner_probe_count != 32U
            || global_accessible_inner_probe_count
                != expected_accessible_inner) {
            throw std::runtime_error("trench accessibility result is incorrect");
        }
    }

    std::cout << std::setprecision(9)
              << "scenario="
              << gasaccess::testing::mock_scenario_name(options.scenario) << '\n'
              << "mpi_ranks=" << process_count << '\n'
              << "process_grid=" << process_grid.x << 'x' << process_grid.y
              << 'x' << process_grid.z << '\n'
              << "owned_atom_count=" << global_owned_atom_count << '\n'
              << "ghost_atom_count_sum=" << global_ghost_atom_count << '\n'
              << "solid_voxel_count=" << global_solid_count << '\n'
              << "outside_voxel_count=" << global_outside_count << '\n'
              << "closed_void_voxel_count=" << global_closed_count << '\n'
              << "accessible_atom_count=" << global_accessible_atom_count << '\n'
              << "inner_probe_count=" << global_inner_probe_count << '\n'
              << "accessible_inner_probe_count="
              << global_accessible_inner_probe_count << '\n'
              << "classification_visited_voxel_count="
              << global_visited_count << '\n'
              << "classification_sent_frontier_count="
              << global_sent_frontier_count << '\n'
              << "classification_received_frontier_count="
              << global_received_frontier_count << '\n'
              << "atom_generation_seconds_max="
              << max_generation_seconds << '\n'
              << "ghost_exchange_seconds_max=" << max_ghost_seconds << '\n'
              << "grid_construction_seconds_max=" << max_grid_seconds << '\n'
              << "voxelization_seconds_max=" << max_voxelization_seconds << '\n'
              << "classification_seconds_max="
              << max_classification_seconds << '\n'
              << "query_seconds_max=" << max_query_seconds << '\n'
              << "query_rate_per_second="
              << (max_query_seconds > 0.0
                      ? static_cast<double>(global_owned_atom_count)
                            / max_query_seconds
                      : 0.0)
              << '\n'
              << "total_seconds_max=" << max_total_seconds << '\n'
              << "peak_rss_bytes_sum=" << global_peak_rss << '\n'
              << "acceptance=pass\n";
    return 0;
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
    try {
        const auto options = parse_options(argument_count, arguments);
        if (options.show_help) {
            if (rank == 0) {
                print_help();
            }
            MPI_Finalize();
            return 0;
        }
        const auto status = run(options, rank, process_count);
        MPI_Finalize();
        return status;
    } catch (const std::exception& error) {
        std::cerr << "rank " << rank << ": " << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
