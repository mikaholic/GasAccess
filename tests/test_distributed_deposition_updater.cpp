#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/deposition_updater.hpp"
#include "gasaccess/distributed_deposition_updater.hpp"
#include "gasaccess/distributed_exterior_classifier.hpp"
#include "gasaccess/exterior_classifier.hpp"
#include "gasaccess/gas_grid.hpp"
#include "gasaccess/mpi_gas_grid.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using gasaccess::Atom;
using gasaccess::AtomView;
using gasaccess::ConnectivityRepairMode;
using gasaccess::DepositionUpdater;
using gasaccess::DistributedDepositionUpdateResult;
using gasaccess::DistributedDepositionUpdater;
using gasaccess::DistributedExteriorClassifier;
using gasaccess::DistributedGasAccessibilityQuery;
using gasaccess::DistributedGasGrid;
using gasaccess::ExteriorClassifier;
using gasaccess::GasAccessibilityQuery;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridDimensions;
using gasaccess::GridSpec;
using gasaccess::MpiDecompositionSpec;
using gasaccess::Point3;
using gasaccess::VoxelCoord;

class TestFailure : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

void require_condition(
    bool condition,
    const char* expression,
    const char* file,
    int line)
{
    if (condition) {
        return;
    }
    std::ostringstream message;
    message << file << ':' << line << ": check failed: " << expression;
    throw TestFailure(message.str());
}

#define REQUIRE(condition) \
    require_condition((condition), #condition, __FILE__, __LINE__)

enum class ExpectedPath {
    Unspecified,
    NoGeometryChange,
    SafeWithoutRepair,
    TopologyFilterSafe,
    DistributedRepair
};

VoxelCoord process_location(int rank, const GridDimensions& process_grid)
{
    const auto unsigned_rank = static_cast<std::uint64_t>(rank);
    return {
        static_cast<std::int64_t>(unsigned_rank % process_grid.x),
        static_cast<std::int64_t>(
            (unsigned_rank / process_grid.x) % process_grid.y),
        static_cast<std::int64_t>(
            unsigned_rank / (process_grid.x * process_grid.y))};
}

GridSpec make_grid_spec()
{
    GridSpec grid_spec{};
    grid_spec.origin = {-1.0, 2.0, 4.0};
    grid_spec.spacing = {0.5, 0.75, 1.25};
    grid_spec.dimensions = {8, 8, 8};
    grid_spec.periodic = {true, true, false};
    grid_spec.reservoir_faces.z_high = true;
    return grid_spec;
}

MpiDecompositionSpec make_decomposition_spec(
    const GridSpec& grid_spec,
    const GridDimensions& process_grid,
    int rank)
{
    const auto location = process_location(rank, process_grid);
    const Point3 global_upper{
        grid_spec.origin.x
            + grid_spec.spacing.x * static_cast<double>(grid_spec.dimensions.x),
        grid_spec.origin.y
            + grid_spec.spacing.y * static_cast<double>(grid_spec.dimensions.y),
        grid_spec.origin.z
            + grid_spec.spacing.z * static_cast<double>(grid_spec.dimensions.z)};
    const Point3 process_lengths{
        (global_upper.x - grid_spec.origin.x)
            / static_cast<double>(process_grid.x),
        (global_upper.y - grid_spec.origin.y)
            / static_cast<double>(process_grid.y),
        (global_upper.z - grid_spec.origin.z)
            / static_cast<double>(process_grid.z)};

    MpiDecompositionSpec decomposition_spec{};
    decomposition_spec.communicator = MPI_COMM_WORLD;
    decomposition_spec.global_lower = grid_spec.origin;
    decomposition_spec.global_upper = global_upper;
    decomposition_spec.local_lower = {
        grid_spec.origin.x + static_cast<double>(location.x) * process_lengths.x,
        grid_spec.origin.y + static_cast<double>(location.y) * process_lengths.y,
        grid_spec.origin.z + static_cast<double>(location.z) * process_lengths.z};
    decomposition_spec.local_upper = {
        decomposition_spec.local_lower.x + process_lengths.x,
        decomposition_spec.local_lower.y + process_lengths.y,
        decomposition_spec.local_lower.z + process_lengths.z};
    decomposition_spec.process_grid = process_grid;
    decomposition_spec.process_location = location;
    decomposition_spec.atom_ghost_distance = 2.0;
    decomposition_spec.maximum_excluded_radius = 1.0;
    return decomposition_spec;
}

std::optional<VoxelCoord> normalize_coordinate(
    VoxelCoord coordinate,
    const GridSpec& grid_spec)
{
    const auto normalize_axis = [](std::int64_t& value,
                                   std::uint64_t dimension,
                                   bool periodic) {
        if (value < 0) {
            if (!periodic) {
                return false;
            }
            value = static_cast<std::int64_t>(dimension - 1);
        } else if (static_cast<std::uint64_t>(value) >= dimension) {
            if (!periodic) {
                return false;
            }
            value = 0;
        }
        return true;
    };
    if (!normalize_axis(coordinate.x, grid_spec.dimensions.x, grid_spec.periodic.x)
        || !normalize_axis(
            coordinate.y,
            grid_spec.dimensions.y,
            grid_spec.periodic.y)
        || !normalize_axis(
            coordinate.z,
            grid_spec.dimensions.z,
            grid_spec.periodic.z)) {
        return std::nullopt;
    }
    return coordinate;
}

std::vector<VoxelCoord> make_barrier(
    const std::vector<VoxelCoord>& openings)
{
    std::vector<VoxelCoord> solid_voxels;
    for (std::int64_t y = 0; y < 8; ++y) {
        for (std::int64_t x = 0; x < 8; ++x) {
            const VoxelCoord voxel_coord{x, y, 4};
            if (std::find(openings.begin(), openings.end(), voxel_coord)
                == openings.end()) {
                solid_voxels.push_back(voxel_coord);
            }
        }
    }
    return solid_voxels;
}

std::vector<Atom> atoms_for_voxels(
    const GasGrid& gas_grid,
    const std::vector<VoxelCoord>& voxel_coords)
{
    std::vector<Atom> atoms;
    atoms.reserve(voxel_coords.size());
    for (const auto& voxel_coord : voxel_coords) {
        atoms.push_back({gas_grid.voxel_center(voxel_coord), 0.0});
    }
    return atoms;
}

void set_initial_solids(
    GasGrid& serial_grid,
    DistributedGasGrid& distributed_grid,
    const std::vector<VoxelCoord>& solid_voxels)
{
    for (const auto& voxel_coord : solid_voxels) {
        serial_grid.set_gas_state(
            serial_grid.voxel_id(voxel_coord),
            GasState::Solid);
        if (distributed_grid.owns(voxel_coord)) {
            distributed_grid.set_owned_gas_state(voxel_coord, GasState::Solid);
        }
    }
}

void set_initial_solids(
    DistributedGasGrid& distributed_grid,
    const std::vector<VoxelCoord>& solid_voxels)
{
    for (const auto& voxel_coord : solid_voxels) {
        if (distributed_grid.owns(voxel_coord)) {
            distributed_grid.set_owned_gas_state(voxel_coord, GasState::Solid);
        }
    }
}

void verify_state_and_queries(
    const GasGrid& serial_grid,
    const DistributedGasGrid& incremental_grid,
    const DistributedGasGrid& full_grid,
    const DistributedDepositionUpdateResult& incremental_result)
{
    const auto& range = incremental_grid.owned_range();
    const auto& grid_spec = incremental_grid.global_grid_spec();
    const GasAccessibilityQuery serial_query(serial_grid);
    const DistributedGasAccessibilityQuery incremental_query(incremental_grid);
    const DistributedGasAccessibilityQuery full_query(full_grid);
    std::array<std::uint64_t, 4> counted_states{};

    for (auto z = range.begin.z; z < range.end.z; ++z) {
        for (auto y = range.begin.y; y < range.end.y; ++y) {
            for (auto x = range.begin.x; x < range.end.x; ++x) {
                const VoxelCoord voxel_coord{x, y, z};
                const auto serial_state = serial_grid.gas_state(voxel_coord);
                REQUIRE(incremental_grid.gas_state(voxel_coord) == serial_state);
                REQUIRE(full_grid.gas_state(voxel_coord) == serial_state);
                ++counted_states[static_cast<std::size_t>(serial_state)];

                const std::array<VoxelCoord, 6> neighbors{{
                    {x - 1, y, z},
                    {x + 1, y, z},
                    {x, y - 1, z},
                    {x, y + 1, z},
                    {x, y, z - 1},
                    {x, y, z + 1}}};
                for (const auto& candidate : neighbors) {
                    const auto neighbor = normalize_coordinate(candidate, grid_spec);
                    if (neighbor.has_value()) {
                        REQUIRE(incremental_grid.gas_state(*neighbor)
                            == serial_grid.gas_state(*neighbor));
                        REQUIRE(full_grid.gas_state(*neighbor)
                            == serial_grid.gas_state(*neighbor));
                    }
                }

                const auto center = incremental_grid.voxel_center(voxel_coord);
                const auto expected_accessible =
                    serial_query.is_site_accessible(center);
                REQUIRE(incremental_query.is_site_accessible(center)
                    == expected_accessible);
                REQUIRE(full_query.is_site_accessible(center)
                    == expected_accessible);
            }
        }
    }

    for (const auto gas_state : {
             GasState::Unclassified,
             GasState::Solid,
             GasState::OutsideAccessible,
             GasState::ClosedVoid}) {
        REQUIRE(incremental_grid.owned_gas_state_count(gas_state)
            == counted_states[static_cast<std::size_t>(gas_state)]);
        REQUIRE(full_grid.owned_gas_state_count(gas_state)
            == counted_states[static_cast<std::size_t>(gas_state)]);
    }
    REQUIRE(incremental_result.classification.local_solid_count
        == counted_states[static_cast<std::size_t>(GasState::Solid)]);
    REQUIRE(incremental_result.classification.local_outside_accessible_count
        == counted_states[
            static_cast<std::size_t>(GasState::OutsideAccessible)]);
    REQUIRE(incremental_result.classification.local_closed_void_count
        == counted_states[static_cast<std::size_t>(GasState::ClosedVoid)]);

    std::uint64_t global_solid = 0;
    std::uint64_t global_outside = 0;
    std::uint64_t global_closed = 0;
    const auto local_solid = counted_states[
        static_cast<std::size_t>(GasState::Solid)];
    const auto local_outside = counted_states[
        static_cast<std::size_t>(GasState::OutsideAccessible)];
    const auto local_closed = counted_states[
        static_cast<std::size_t>(GasState::ClosedVoid)];
    MPI_Allreduce(
        &local_solid,
        &global_solid,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        &local_outside,
        &global_outside,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        &local_closed,
        &global_closed,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    REQUIRE(global_solid == serial_grid.gas_state_count(GasState::Solid));
    REQUIRE(global_outside
        == serial_grid.gas_state_count(GasState::OutsideAccessible));
    REQUIRE(global_closed == serial_grid.gas_state_count(GasState::ClosedVoid));
}

void verify_expected_path(
    const DistributedDepositionUpdateResult& result,
    ExpectedPath expected_path)
{
    if (expected_path == ExpectedPath::Unspecified) {
        return;
    }
    if (expected_path == ExpectedPath::NoGeometryChange) {
        REQUIRE(!result.geometry_changed());
        REQUIRE(!result.used_distributed_repair());
        return;
    }
    REQUIRE(result.geometry_changed());
    if (expected_path == ExpectedPath::SafeWithoutRepair
        || expected_path == ExpectedPath::TopologyFilterSafe) {
        REQUIRE(!result.used_distributed_repair());
        if (expected_path == ExpectedPath::TopologyFilterSafe) {
            REQUIRE(result.topology_filter_performed);
            REQUIRE(result.topology_filter_declared_safe);
        }
    } else {
        REQUIRE(result.used_distributed_repair());
    }
}

void run_sequence(
    const GridSpec& grid_spec,
    const GridDimensions& process_grid,
    const std::vector<VoxelCoord>& initial_solids,
    const std::vector<std::vector<VoxelCoord>>& deposition_events,
    const std::vector<ExpectedPath>& expected_paths,
    int rank)
{
    REQUIRE(expected_paths.empty()
        || expected_paths.size() == deposition_events.size());
    GasGrid serial_grid(grid_spec);
    DistributedGasGrid incremental_grid(
        grid_spec,
        make_decomposition_spec(grid_spec, process_grid, rank));
    DistributedGasGrid full_grid(
        grid_spec,
        make_decomposition_spec(grid_spec, process_grid, rank));
    set_initial_solids(serial_grid, incremental_grid, initial_solids);
    set_initial_solids(full_grid, initial_solids);

    ExteriorClassifier().classify(serial_grid);
    DistributedExteriorClassifier incremental_classifier;
    DistributedExteriorClassifier full_classifier;
    incremental_classifier.classify(incremental_grid);
    full_classifier.classify(full_grid);

    DepositionUpdater serial_updater(
        0.0,
        ConnectivityRepairMode::AffectedRegion);
    DistributedDepositionUpdater incremental_updater(
        0.0,
        ConnectivityRepairMode::AffectedRegion);
    DistributedDepositionUpdater full_updater(
        0.0,
        ConnectivityRepairMode::FullReclassification);

    for (std::size_t event_index = 0;
         event_index < deposition_events.size();
         ++event_index) {
        const auto atoms = atoms_for_voxels(
            serial_grid,
            deposition_events[event_index]);
        const AtomView atom_view{atoms.data(), atoms.size()};
        const auto serial_result = serial_updater.apply_deposition(
            serial_grid,
            atom_view);
        const auto incremental_result = incremental_updater.apply_deposition(
            incremental_grid,
            atom_view);
        const auto full_result = full_updater.apply_deposition(
            full_grid,
            atom_view);

        REQUIRE(incremental_result.global_newly_solid_count
            == serial_result.newly_solid_count);
        REQUIRE(full_result.global_newly_solid_count
            == serial_result.newly_solid_count);
        REQUIRE(full_result.used_full_reclassification()
            == full_result.geometry_changed());
        REQUIRE(incremental_result.changed_owned_voxel_coords
            == full_result.changed_owned_voxel_coords);
        verify_expected_path(
            incremental_result,
            expected_paths.empty()
                ? ExpectedPath::Unspecified
                : expected_paths[event_index]);
        verify_state_and_queries(
            serial_grid,
            incremental_grid,
            full_grid,
            incremental_result);

        std::uint64_t global_incremental_closed = 0;
        std::uint64_t global_sent = 0;
        std::uint64_t global_received = 0;
        MPI_Allreduce(
            &incremental_result.local_repair_closed_voxel_count,
            &global_incremental_closed,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            MPI_COMM_WORLD);
        MPI_Allreduce(
            &incremental_result.sent_frontier_entry_count,
            &global_sent,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            MPI_COMM_WORLD);
        MPI_Allreduce(
            &incremental_result.received_frontier_entry_count,
            &global_received,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            MPI_COMM_WORLD);
        REQUIRE(global_incremental_closed
            == serial_result.repair_closed_voxel_count);
        REQUIRE(global_sent == global_received);
    }
}

GridDimensions local_pinch_layout(int size)
{
    if (size == 4) {
        return {2, 2, 1};
    }
    return {static_cast<std::uint64_t>(size), 1, 1};
}

void test_no_change_and_safe_updates(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const GridDimensions x_layout{static_cast<std::uint64_t>(size), 1, 1};
    run_sequence(
        grid_spec,
        x_layout,
        make_barrier({}),
        {{{0, 0, 4}}, {{2, 2, 2}}},
        {ExpectedPath::NoGeometryChange, ExpectedPath::SafeWithoutRepair},
        rank);

    run_sequence(
        grid_spec,
        local_pinch_layout(size),
        {},
        {{{2, 2, 3}}},
        {ExpectedPath::TopologyFilterSafe},
        rank);
}

void test_local_and_rank_boundary_pinch_off(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const VoxelCoord local_opening{2, 2, 4};
    run_sequence(
        grid_spec,
        local_pinch_layout(size),
        make_barrier({local_opening}),
        {{local_opening}},
        {ExpectedPath::DistributedRepair},
        rank);

    const VoxelCoord boundary_opening{3, 3, 4};
    const GridDimensions z_layout{1, 1, static_cast<std::uint64_t>(size)};
    run_sequence(
        grid_spec,
        z_layout,
        make_barrier({boundary_opening}),
        {{boundary_opening}},
        {ExpectedPath::DistributedRepair},
        rank);
}

void test_periodic_seam_and_multi_voxel_closure(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const GridDimensions x_layout{static_cast<std::uint64_t>(size), 1, 1};
    const VoxelCoord seam_opening{0, 3, 4};
    run_sequence(
        grid_spec,
        x_layout,
        make_barrier({seam_opening}),
        {{seam_opening}},
        {ExpectedPath::DistributedRepair},
        rank);

    const VoxelCoord first_opening{1, 2, 4};
    const VoxelCoord second_opening{6, 5, 4};
    run_sequence(
        grid_spec,
        x_layout,
        make_barrier({first_opening, second_opening}),
        {{first_opening, second_opening}},
        {ExpectedPath::DistributedRepair},
        rank);
}

void test_deterministic_deposition_sequence(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    std::vector<VoxelCoord> plane_voxels;
    for (std::int64_t y = 0; y < 8; ++y) {
        for (std::int64_t x = 0; x < 8; ++x) {
            plane_voxels.push_back({x, y, 4});
        }
    }
    std::mt19937_64 random_generator(193731U);
    std::shuffle(
        plane_voxels.begin(),
        plane_voxels.end(),
        random_generator);

    std::vector<std::vector<VoxelCoord>> events;
    events.reserve(plane_voxels.size());
    for (const auto& voxel_coord : plane_voxels) {
        events.push_back({voxel_coord});
    }
    const GridDimensions x_layout{static_cast<std::uint64_t>(size), 1, 1};
    run_sequence(grid_spec, x_layout, {}, events, {}, rank);
}

}  // namespace

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_failure_count = 0;
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"no-change and safe updates", [&]() {
             test_no_change_and_safe_updates(rank, size);
         }},
        {"local and rank-boundary pinch-off", [&]() {
             test_local_and_rank_boundary_pinch_off(rank, size);
         }},
        {"periodic seam and multi-voxel closure", [&]() {
             test_periodic_seam_and_multi_voxel_closure(rank, size);
         }},
        {"deterministic deposition sequence", [&]() {
             test_deterministic_deposition_sequence(rank, size);
         }}
    };

    for (const auto& test : tests) {
        try {
            test.second();
            if (rank == 0) {
                std::cout << "[PASS] " << test.first << '\n';
            }
        } catch (const std::exception& error) {
            ++local_failure_count;
            std::cerr << "[rank " << rank << "] [FAIL] " << test.first
                      << ": " << error.what() << '\n';
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int global_failure_count = 0;
    MPI_Allreduce(
        &local_failure_count,
        &global_failure_count,
        1,
        MPI_INT,
        MPI_SUM,
        MPI_COMM_WORLD);
    if (rank == 0 && global_failure_count == 0) {
        std::cout << tests.size() << " distributed update test groups passed on "
                  << size << " rank(s)\n";
    }
    MPI_Finalize();
    return global_failure_count == 0 ? 0 : 1;
}
