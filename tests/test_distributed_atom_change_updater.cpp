#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/atom_change_updater.hpp"
#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/distributed_atom_change_updater.hpp"
#include "gasaccess/distributed_exterior_classifier.hpp"
#include "gasaccess/exterior_classifier.hpp"
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

using gasaccess::AccessibilityRepairKind;
using gasaccess::Atom;
using gasaccess::AtomChangeBatch;
using gasaccess::AtomChangeRepairMode;
using gasaccess::AtomChangeUpdateResult;
using gasaccess::AtomChangeUpdater;
using gasaccess::AtomView;
using gasaccess::AtomVoxelizer;
using gasaccess::DistributedAtomChangeUpdateResult;
using gasaccess::DistributedAtomChangeUpdater;
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

bool voxel_coord_less(const VoxelCoord& lhs, const VoxelCoord& rhs) noexcept
{
    if (lhs.z != rhs.z) {
        return lhs.z < rhs.z;
    }
    if (lhs.y != rhs.y) {
        return lhs.y < rhs.y;
    }
    return lhs.x < rhs.x;
}

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

GridSpec make_grid_spec(
    std::uint64_t y = 8,
    std::uint64_t z = 8)
{
    GridSpec grid_spec{};
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {32, y, z};
    grid_spec.reservoir_faces.z_high = true;
    return grid_spec;
}

MpiDecompositionSpec make_decomposition_spec(
    const GridSpec& grid_spec,
    int rank,
    int size)
{
    const GridDimensions process_grid{
        static_cast<std::uint64_t>(size), 1, 1};
    const auto location = process_location(rank, process_grid);
    const Point3 global_upper{
        grid_spec.origin.x
            + grid_spec.spacing.x * static_cast<double>(grid_spec.dimensions.x),
        grid_spec.origin.y
            + grid_spec.spacing.y * static_cast<double>(grid_spec.dimensions.y),
        grid_spec.origin.z
            + grid_spec.spacing.z * static_cast<double>(grid_spec.dimensions.z)};
    const double local_x_length =
        (global_upper.x - grid_spec.origin.x) / static_cast<double>(size);

    MpiDecompositionSpec decomposition_spec{};
    decomposition_spec.communicator = MPI_COMM_WORLD;
    decomposition_spec.global_lower = grid_spec.origin;
    decomposition_spec.global_upper = global_upper;
    decomposition_spec.local_lower = {
        grid_spec.origin.x + static_cast<double>(location.x) * local_x_length,
        grid_spec.origin.y,
        grid_spec.origin.z};
    decomposition_spec.local_upper = {
        decomposition_spec.local_lower.x + local_x_length,
        global_upper.y,
        global_upper.z};
    decomposition_spec.process_grid = process_grid;
    decomposition_spec.process_location = location;
    decomposition_spec.atom_ghost_distance = 0.0;
    decomposition_spec.maximum_excluded_radius = 0.0;
    return decomposition_spec;
}

std::optional<VoxelCoord> normalize_coordinate(
    VoxelCoord coordinate,
    const GridSpec& grid_spec)
{
    const auto normalize_axis = [](
                                    std::int64_t& value,
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

std::vector<VoxelCoord> all_solid_except(
    const GridSpec& grid_spec,
    const std::vector<VoxelCoord>& open_voxels)
{
    std::vector<VoxelCoord> solids;
    for (std::int64_t z = 0;
         z < static_cast<std::int64_t>(grid_spec.dimensions.z);
         ++z) {
        for (std::int64_t y = 0;
             y < static_cast<std::int64_t>(grid_spec.dimensions.y);
             ++y) {
            for (std::int64_t x = 0;
                 x < static_cast<std::int64_t>(grid_spec.dimensions.x);
                 ++x) {
                const VoxelCoord voxel_coord{x, y, z};
                if (std::find(
                        open_voxels.begin(),
                        open_voxels.end(),
                        voxel_coord) == open_voxels.end()) {
                    solids.push_back(voxel_coord);
                }
            }
        }
    }
    return solids;
}

std::vector<VoxelCoord> plane(
    const GridSpec& grid_spec,
    std::int64_t z)
{
    std::vector<VoxelCoord> result;
    for (std::int64_t y = 0;
         y < static_cast<std::int64_t>(grid_spec.dimensions.y);
         ++y) {
        for (std::int64_t x = 0;
             x < static_cast<std::int64_t>(grid_spec.dimensions.x);
             ++x) {
            result.push_back({x, y, z});
        }
    }
    return result;
}

void initialize_grids(
    GasGrid& serial_grid,
    DistributedGasGrid& incremental_grid,
    DistributedGasGrid& full_grid,
    const std::vector<VoxelCoord>& atom_voxel_coords)
{
    const auto atoms = atoms_for_voxels(serial_grid, atom_voxel_coords);
    const AtomView atom_view{atoms.data(), atoms.size()};
    AtomVoxelizer(0.0).voxelize(serial_grid, atom_view);
    incremental_grid.voxelize_owned_atoms(atom_view, 0.0);
    full_grid.voxelize_owned_atoms(atom_view, 0.0);
    ExteriorClassifier{}.classify(serial_grid);
    DistributedExteriorClassifier{}.classify(incremental_grid);
    DistributedExteriorClassifier{}.classify(full_grid);
}

void require_sorted_unique(const std::vector<VoxelCoord>& voxel_coords)
{
    REQUIRE(std::is_sorted(
        voxel_coords.begin(),
        voxel_coords.end(),
        voxel_coord_less));
    REQUIRE(std::adjacent_find(voxel_coords.begin(), voxel_coords.end())
        == voxel_coords.end());
}

std::vector<VoxelCoord> owned_serial_changes(
    const GasGrid& serial_grid,
    const DistributedGasGrid& distributed_grid,
    const AtomChangeUpdateResult& serial_result)
{
    std::vector<VoxelCoord> result;
    for (const auto voxel_id : serial_result.changed_voxel_ids) {
        const auto voxel_coord = serial_grid.voxel_coord(voxel_id);
        if (distributed_grid.owns(voxel_coord)) {
            result.push_back(voxel_coord);
        }
    }
    std::sort(result.begin(), result.end(), voxel_coord_less);
    return result;
}

void verify_states_ghosts_and_queries(
    const GasGrid& serial_grid,
    const DistributedGasGrid& incremental_grid,
    const DistributedGasGrid& full_grid,
    const DistributedAtomChangeUpdateResult& incremental_result,
    const DistributedAtomChangeUpdateResult& full_result)
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
                const auto expected_state = serial_grid.gas_state(voxel_coord);
                REQUIRE(incremental_grid.gas_state(voxel_coord)
                    == expected_state);
                REQUIRE(full_grid.gas_state(voxel_coord) == expected_state);
                REQUIRE(incremental_grid.owned_blocker_count(voxel_coord)
                    == serial_grid.blocker_count(voxel_coord));
                REQUIRE(full_grid.owned_blocker_count(voxel_coord)
                    == serial_grid.blocker_count(voxel_coord));
                ++counted_states[static_cast<std::size_t>(expected_state)];

                const std::array<VoxelCoord, 6> neighbors{{
                    {x - 1, y, z},
                    {x + 1, y, z},
                    {x, y - 1, z},
                    {x, y + 1, z},
                    {x, y, z - 1},
                    {x, y, z + 1}}};
                for (const auto& candidate : neighbors) {
                    const auto neighbor = normalize_coordinate(
                        candidate,
                        grid_spec);
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

    for (const auto state : {
             GasState::Unclassified,
             GasState::Solid,
             GasState::OutsideAccessible,
             GasState::ClosedVoid}) {
        const auto expected = counted_states[static_cast<std::size_t>(state)];
        REQUIRE(incremental_grid.owned_gas_state_count(state) == expected);
        REQUIRE(full_grid.owned_gas_state_count(state) == expected);
    }
    REQUIRE(incremental_result.classification.local_solid_count
        == counted_states[static_cast<std::size_t>(GasState::Solid)]);
    REQUIRE(incremental_result.classification.local_outside_accessible_count
        == counted_states[
            static_cast<std::size_t>(GasState::OutsideAccessible)]);
    REQUIRE(incremental_result.classification.local_closed_void_count
        == counted_states[static_cast<std::size_t>(GasState::ClosedVoid)]);
    REQUIRE(full_result.classification.local_solid_count
        == counted_states[static_cast<std::size_t>(GasState::Solid)]);
    REQUIRE(full_result.classification.local_outside_accessible_count
        == counted_states[
            static_cast<std::size_t>(GasState::OutsideAccessible)]);
    REQUIRE(full_result.classification.local_closed_void_count
        == counted_states[static_cast<std::size_t>(GasState::ClosedVoid)]);

    const std::array<std::uint64_t, 3> local_counts{
        counted_states[static_cast<std::size_t>(GasState::Solid)],
        counted_states[static_cast<std::size_t>(GasState::OutsideAccessible)],
        counted_states[static_cast<std::size_t>(GasState::ClosedVoid)]};
    std::array<std::uint64_t, 3> global_counts{};
    MPI_Allreduce(
        local_counts.data(),
        global_counts.data(),
        static_cast<int>(global_counts.size()),
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    REQUIRE(global_counts[0]
        == serial_grid.gas_state_count(GasState::Solid));
    REQUIRE(global_counts[1]
        == serial_grid.gas_state_count(GasState::OutsideAccessible));
    REQUIRE(global_counts[2]
        == serial_grid.gas_state_count(GasState::ClosedVoid));
}

struct CoordinateBatch {
    std::vector<VoxelCoord> added{};
    std::vector<VoxelCoord> removed{};
};

using ResultVerifier = std::function<void(
    std::size_t,
    const DistributedAtomChangeUpdateResult&,
    const AtomChangeUpdateResult&)>;

void run_sequence(
    const GridSpec& grid_spec,
    const std::vector<VoxelCoord>& initial_atom_voxels,
    const std::vector<CoordinateBatch>& batches,
    int rank,
    int size,
    const ResultVerifier& result_verifier = {})
{
    GasGrid serial_grid(grid_spec);
    DistributedGasGrid incremental_grid(
        grid_spec,
        make_decomposition_spec(grid_spec, rank, size));
    DistributedGasGrid full_grid(
        grid_spec,
        make_decomposition_spec(grid_spec, rank, size));
    initialize_grids(
        serial_grid,
        incremental_grid,
        full_grid,
        initial_atom_voxels);

    AtomChangeUpdater serial_updater(0.0);
    DistributedAtomChangeUpdater incremental_updater(0.0);
    DistributedAtomChangeUpdater full_updater(
        0.0,
        AtomChangeRepairMode::FullReclassification);

    for (std::size_t event_index = 0;
         event_index < batches.size();
         ++event_index) {
        const auto added_atoms = atoms_for_voxels(
            serial_grid,
            batches[event_index].added);
        const auto removed_atoms = atoms_for_voxels(
            serial_grid,
            batches[event_index].removed);
        const AtomChangeBatch changes{
            {added_atoms.data(), added_atoms.size()},
            {removed_atoms.data(), removed_atoms.size()}};
        const auto serial_result = serial_updater.apply_atom_changes(
            serial_grid,
            changes);
        const auto incremental_result = incremental_updater.apply_atom_changes(
            incremental_grid,
            changes);
        const auto full_result = full_updater.apply_atom_changes(
            full_grid,
            changes);

        REQUIRE(incremental_result.global_blocker_count_changed_voxel_count
            == serial_result.blocker_count_changed_voxel_count);
        REQUIRE(full_result.global_blocker_count_changed_voxel_count
            == serial_result.blocker_count_changed_voxel_count);
        REQUIRE(incremental_result.global_newly_solid_count
            == serial_result.newly_solid_count);
        REQUIRE(incremental_result.global_newly_gas_count
            == serial_result.newly_gas_count);
        REQUIRE(full_result.global_newly_solid_count
            == serial_result.newly_solid_count);
        REQUIRE(full_result.global_newly_gas_count
            == serial_result.newly_gas_count);
        REQUIRE(incremental_result.repair_kind == serial_result.repair_kind);
        REQUIRE(full_result.repair_kind
            == (full_result.geometry_changed()
                ? AccessibilityRepairKind::FullReclassification
                : AccessibilityRepairKind::None));
        REQUIRE(full_result.used_full_reclassification()
            == full_result.geometry_changed());

        const auto expected_local_changes = owned_serial_changes(
            serial_grid,
            incremental_grid,
            serial_result);
        REQUIRE(incremental_result.changed_owned_voxel_coords
            == expected_local_changes);
        REQUIRE(full_result.changed_owned_voxel_coords
            == expected_local_changes);
        require_sorted_unique(
            incremental_result.changed_owned_voxel_coords);
        require_sorted_unique(full_result.changed_owned_voxel_coords);

        const std::array<std::uint64_t, 6> local_metrics{
            incremental_result.local_repair_closed_voxel_count,
            incremental_result.local_repair_opened_voxel_count,
            incremental_result.closing_sent_frontier_entry_count,
            incremental_result.closing_received_frontier_entry_count,
            incremental_result.opening_sent_frontier_entry_count,
            incremental_result.opening_received_frontier_entry_count};
        std::array<std::uint64_t, 6> global_metrics{};
        MPI_Allreduce(
            local_metrics.data(),
            global_metrics.data(),
            static_cast<int>(global_metrics.size()),
            MPI_UINT64_T,
            MPI_SUM,
            MPI_COMM_WORLD);
        REQUIRE(global_metrics[0] == serial_result.repair_closed_voxel_count);
        REQUIRE(global_metrics[1] == serial_result.repair_opened_voxel_count);
        REQUIRE(global_metrics[2] == global_metrics[3]);
        REQUIRE(global_metrics[4] == global_metrics[5]);

        verify_states_ghosts_and_queries(
            serial_grid,
            incremental_grid,
            full_grid,
            incremental_result,
            full_result);
        if (result_verifier) {
            result_verifier(event_index, incremental_result, serial_result);
        }
    }
}

void test_empty_cancellation_and_cross_rank_move(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const VoxelCoord old_position{1, 3, 3};
    const VoxelCoord new_position{30, 3, 3};
    run_sequence(
        grid_spec,
        {old_position},
        {{{old_position}, {old_position}},
         {{new_position}, {old_position}}},
        rank,
        size,
        [](std::size_t event_index,
           const DistributedAtomChangeUpdateResult& distributed_result,
           const AtomChangeUpdateResult&) {
            if (event_index == 0) {
                REQUIRE(!distributed_result.geometry_changed());
                REQUIRE(distributed_result.repair_kind
                    == AccessibilityRepairKind::None);
            } else {
                REQUIRE(distributed_result.global_newly_solid_count == 1);
                REQUIRE(distributed_result.global_newly_gas_count == 1);
                REQUIRE(distributed_result.repair_kind
                    == AccessibilityRepairKind::Mixed);
            }
        });
}

void test_collective_mixed_rollback(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const auto block_width = std::int64_t{32}
        / static_cast<std::int64_t>(size);
    std::vector<VoxelCoord> occupied_voxels;
    for (int owner_rank = 0; owner_rank < size; ++owner_rank) {
        occupied_voxels.push_back({
            static_cast<std::int64_t>(owner_rank) * block_width + 1,
            2,
            2});
    }
    GasGrid serial_grid(grid_spec);
    DistributedGasGrid distributed_grid(
        grid_spec,
        make_decomposition_spec(grid_spec, rank, size));
    const auto initial_atoms = atoms_for_voxels(
        serial_grid,
        occupied_voxels);
    AtomVoxelizer(0.0).voxelize(
        serial_grid,
        {initial_atoms.data(), initial_atoms.size()});
    distributed_grid.voxelize_owned_atoms(
        {initial_atoms.data(), initial_atoms.size()},
        0.0);
    ExteriorClassifier{}.classify(serial_grid);
    DistributedExteriorClassifier{}.classify(distributed_grid);

    const VoxelCoord invalid_removal{2, 2, 2};
    const auto rank_local_removal = rank == 0
        ? invalid_removal
        : occupied_voxels[static_cast<std::size_t>(rank)];
    const VoxelCoord rank_local_addition{
        static_cast<std::int64_t>(rank) * block_width + 2,
        4,
        2};
    const auto removed = atoms_for_voxels(
        serial_grid,
        {rank_local_removal});
    const auto added = atoms_for_voxels(
        serial_grid,
        {rank_local_addition});
    bool threw = false;
    try {
        static_cast<void>(DistributedAtomChangeUpdater(0.0).apply_atom_changes(
            distributed_grid,
            {{added.data(), added.size()},
             {removed.data(), removed.size()}}));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    REQUIRE(threw);
    for (const auto& occupied : occupied_voxels) {
        if (distributed_grid.owns(occupied)) {
            REQUIRE(distributed_grid.owned_blocker_count(occupied) == 1);
            REQUIRE(distributed_grid.gas_state(occupied) == GasState::Solid);
        }
    }
    if (distributed_grid.owns(rank_local_addition)) {
        REQUIRE(distributed_grid.owned_blocker_count(rank_local_addition) == 0);
        REQUIRE(distributed_grid.gas_state(rank_local_addition)
            == GasState::OutsideAccessible);
    }
}

std::vector<VoxelCoord> two_cavity_solids(const GridSpec& grid_spec)
{
    std::vector<VoxelCoord> open_voxels;
    for (const auto center_x : {std::int64_t{2}, std::int64_t{29}}) {
        for (std::int64_t z = 2; z <= 4; ++z) {
            for (std::int64_t y = 1; y <= 3; ++y) {
                for (std::int64_t x = center_x - 1;
                     x <= center_x + 1;
                     ++x) {
                    open_voxels.push_back({x, y, z});
                }
            }
        }
    }
    for (std::int64_t z = 5; z <= 7; ++z) {
        open_voxels.push_back({2, 2, z});
    }
    for (std::int64_t z = 6; z <= 7; ++z) {
        open_voxels.push_back({29, 2, z});
    }
    return all_solid_except(grid_spec, open_voxels);
}

void test_close_and_open_on_different_ranks(int rank, int size)
{
    const auto grid_spec = make_grid_spec(5, 8);
    run_sequence(
        grid_spec,
        two_cavity_solids(grid_spec),
        {{{{2, 2, 5}}, {{29, 2, 5}}}},
        rank,
        size,
        [size](std::size_t,
               const DistributedAtomChangeUpdateResult& distributed_result,
               const AtomChangeUpdateResult& serial_result) {
            REQUIRE(serial_result.repair_closed_voxel_count > 0);
            REQUIRE(serial_result.repair_opened_voxel_count > 0);
            REQUIRE(distributed_result.used_distributed_closing_repair());
            REQUIRE(distributed_result.used_distributed_opening_repair());
            if (size > 1) {
                std::uint64_t local_transition_rank =
                    distributed_result.local_newly_solid_count != 0
                        || distributed_result.local_newly_gas_count != 0
                    ? 1U
                    : 0U;
                std::uint64_t transition_ranks = 0;
                MPI_Allreduce(
                    &local_transition_rank,
                    &transition_ranks,
                    1,
                    MPI_UINT64_T,
                    MPI_SUM,
                    MPI_COMM_WORLD);
                REQUIRE(transition_ranks == 2);
            }
        });
}

void test_post_closing_boundary_seed_check(int rank, int size)
{
    auto grid_spec = make_grid_spec();
    grid_spec.periodic.x = false;
    const auto block_width = std::int64_t{32}
        / static_cast<std::int64_t>(size);
    const std::int64_t target_x = size == 1 ? 8 : block_width;
    const VoxelCoord target{target_x, 3, 3};
    const VoxelCoord channel_hole{15, 1, 4};
    auto initial_solids = plane(grid_spec, 4);
    initial_solids.erase(
        std::remove(initial_solids.begin(), initial_solids.end(), channel_hole),
        initial_solids.end());
    initial_solids.push_back(target);
    initial_solids.push_back({target.x + 1, target.y, target.z});
    initial_solids.push_back({target.x, target.y - 1, target.z});
    initial_solids.push_back({target.x, target.y + 1, target.z});
    initial_solids.push_back({target.x, target.y, target.z - 1});
    initial_solids.push_back({target.x, target.y, target.z + 1});

    run_sequence(
        grid_spec,
        initial_solids,
        {{{channel_hole}, {target}}},
        rank,
        size,
        [target](std::size_t,
                 const DistributedAtomChangeUpdateResult& distributed_result,
                 const AtomChangeUpdateResult& serial_result) {
            REQUIRE(serial_result.repair_closed_voxel_count > 0);
            REQUIRE(serial_result.repair_opened_voxel_count == 0);
            REQUIRE(!distributed_result.used_distributed_opening_repair());
            REQUIRE(distributed_result.global_newly_gas_count == 1);
            static_cast<void>(target);
        });
}

void test_channel_swap_reopens_transient_closure(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const VoxelCoord old_channel{15, 3, 4};
    const VoxelCoord new_channel{16, 3, 4};
    auto initial_solids = plane(grid_spec, 4);
    initial_solids.erase(
        std::remove(initial_solids.begin(), initial_solids.end(), old_channel),
        initial_solids.end());
    run_sequence(
        grid_spec,
        initial_solids,
        {{{old_channel}, {new_channel}}},
        rank,
        size,
        [size](std::size_t,
               const DistributedAtomChangeUpdateResult& distributed_result,
               const AtomChangeUpdateResult& serial_result) {
            REQUIRE(serial_result.repair_closed_voxel_count > 0);
            REQUIRE(serial_result.repair_opened_voxel_count > 0);
            REQUIRE(serial_result.changed_voxel_ids.size() == 2);
            REQUIRE(distributed_result.closing_participating_rank_count
                == static_cast<std::uint64_t>(size));
            REQUIRE(distributed_result.opening_participating_rank_count
                == static_cast<std::uint64_t>(size));
        });
}

std::vector<VoxelCoord> local_close_all_rank_open_solids(
    const GridSpec& grid_spec)
{
    std::vector<VoxelCoord> open_voxels;
    for (std::int64_t z = 0; z <= 7; ++z) {
        for (std::int64_t y = 0; y < 8; ++y) {
            for (std::int64_t x = 0; x < 32; ++x) {
                open_voxels.push_back({x, y, z});
            }
        }
    }
    for (std::int64_t x = 1; x <= 3; ++x) {
        for (std::int64_t y = 1; y <= 3; ++y) {
            open_voxels.push_back({x, y, 9});
        }
    }
    open_voxels.push_back({2, 2, 10});
    open_voxels.push_back({2, 2, 11});
    open_voxels.push_back({15, 6, 9});
    open_voxels.push_back({15, 6, 10});
    open_voxels.push_back({15, 6, 11});
    return all_solid_except(grid_spec, open_voxels);
}

void test_local_closing_and_all_rank_opening(int rank, int size)
{
    const auto grid_spec = make_grid_spec(8, 12);
    const VoxelCoord local_plug{2, 2, 10};
    const VoxelCoord global_plug{15, 6, 8};
    run_sequence(
        grid_spec,
        local_close_all_rank_open_solids(grid_spec),
        {{{local_plug}, {global_plug}}},
        rank,
        size,
        [size](std::size_t,
               const DistributedAtomChangeUpdateResult& distributed_result,
               const AtomChangeUpdateResult& serial_result) {
            REQUIRE(serial_result.repair_closed_voxel_count == 9);
            REQUIRE(serial_result.repair_opened_voxel_count
                > 32U * 8U * 12U / 2U);
            REQUIRE(distributed_result.closing_participating_rank_count == 1);
            REQUIRE(distributed_result.opening_participating_rank_count
                == static_cast<std::uint64_t>(size));
            REQUIRE(distributed_result.local_opening_visited_voxel_count > 0);
        });
}

std::vector<VoxelCoord> periodic_two_cavity_solids(const GridSpec& grid_spec)
{
    std::vector<VoxelCoord> open_voxels;
    for (const auto center_y : {std::int64_t{2}, std::int64_t{7}}) {
        for (const auto x : {
                 std::int64_t{31}, std::int64_t{0}, std::int64_t{1}}) {
            for (std::int64_t y = center_y - 1;
                 y <= center_y + 1;
                 ++y) {
                for (std::int64_t z = 2; z <= 4; ++z) {
                    open_voxels.push_back({x, y, z});
                }
            }
        }
    }
    for (std::int64_t z = 5; z <= 7; ++z) {
        open_voxels.push_back({0, 2, z});
    }
    for (std::int64_t z = 6; z <= 7; ++z) {
        open_voxels.push_back({0, 7, z});
    }
    return all_solid_except(grid_spec, open_voxels);
}

void test_both_passes_cross_periodic_seam(int rank, int size)
{
    auto grid_spec = make_grid_spec(10, 8);
    grid_spec.periodic.x = true;
    run_sequence(
        grid_spec,
        periodic_two_cavity_solids(grid_spec),
        {{{{0, 2, 5}}, {{0, 7, 5}}}},
        rank,
        size,
        [size](std::size_t,
               const DistributedAtomChangeUpdateResult& distributed_result,
               const AtomChangeUpdateResult& serial_result) {
            REQUIRE(serial_result.repair_closed_voxel_count > 0);
            REQUIRE(serial_result.repair_opened_voxel_count > 0);
            if (size > 1) {
                REQUIRE(distributed_result.closing_participating_rank_count
                    >= 2);
                REQUIRE(distributed_result.opening_participating_rank_count
                    >= 2);
                std::array<std::uint64_t, 2> local_sent{
                    distributed_result.closing_sent_frontier_entry_count,
                    distributed_result.opening_sent_frontier_entry_count};
                std::array<std::uint64_t, 2> global_sent{};
                MPI_Allreduce(
                    local_sent.data(),
                    global_sent.data(),
                    static_cast<int>(global_sent.size()),
                    MPI_UINT64_T,
                    MPI_SUM,
                    MPI_COMM_WORLD);
                REQUIRE(global_sent[0] > 0);
                REQUIRE(global_sent[1] > 0);
            }
        });
}

void test_deterministic_mixed_sequence(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    GasGrid coordinate_grid(grid_spec);
    std::vector<bool> occupied(
        static_cast<std::size_t>(coordinate_grid.voxel_count()),
        false);
    std::mt19937_64 random_generator(971413U);
    for (std::size_t index = 0; index < occupied.size(); ++index) {
        occupied[index] = random_generator() % 4U == 0;
    }

    std::vector<VoxelCoord> initial_atoms;
    for (gasaccess::VoxelId voxel_id = 0;
         voxel_id < coordinate_grid.voxel_count();
         ++voxel_id) {
        if (occupied[static_cast<std::size_t>(voxel_id)]) {
            initial_atoms.push_back(coordinate_grid.voxel_coord(voxel_id));
        }
    }

    std::vector<CoordinateBatch> batches;
    for (std::size_t event_index = 0; event_index < 24; ++event_index) {
        std::vector<gasaccess::VoxelId> occupied_ids;
        std::vector<gasaccess::VoxelId> empty_ids;
        for (gasaccess::VoxelId voxel_id = 0;
             voxel_id < coordinate_grid.voxel_count();
             ++voxel_id) {
            auto& destination = occupied[static_cast<std::size_t>(voxel_id)]
                ? occupied_ids
                : empty_ids;
            destination.push_back(voxel_id);
        }
        std::shuffle(
            occupied_ids.begin(),
            occupied_ids.end(),
            random_generator);
        std::shuffle(
            empty_ids.begin(),
            empty_ids.end(),
            random_generator);
        REQUIRE(!occupied_ids.empty());
        REQUIRE(!empty_ids.empty());
        CoordinateBatch batch{};
        batch.removed.push_back(
            coordinate_grid.voxel_coord(occupied_ids.front()));
        batch.added.push_back(
            coordinate_grid.voxel_coord(empty_ids.front()));
        occupied[static_cast<std::size_t>(occupied_ids.front())] = false;
        occupied[static_cast<std::size_t>(empty_ids.front())] = true;
        batches.push_back(std::move(batch));
    }

    run_sequence(grid_spec, initial_atoms, batches, rank, size);
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
        {"empty, cancellation, and cross-rank move", [&]() {
             test_empty_cancellation_and_cross_rank_move(rank, size);
         }},
        {"collective mixed rollback", [&]() {
             test_collective_mixed_rollback(rank, size);
         }},
        {"close and open on different ranks", [&]() {
             test_close_and_open_on_different_ranks(rank, size);
         }},
        {"post-closing boundary seed check", [&]() {
             test_post_closing_boundary_seed_check(rank, size);
         }},
        {"channel swap reopens transient closure", [&]() {
             test_channel_swap_reopens_transient_closure(rank, size);
         }},
        {"local closing and all-rank opening", [&]() {
             test_local_closing_and_all_rank_opening(rank, size);
         }},
        {"both passes cross periodic seam", [&]() {
             test_both_passes_cross_periodic_seam(rank, size);
         }},
        {"deterministic mixed sequence", [&]() {
             test_deterministic_mixed_sequence(rank, size);
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
        std::cout << tests.size()
                  << " distributed atom-change test groups passed on "
                  << size << " rank(s)\n";
    }
    MPI_Finalize();
    return global_failure_count == 0 ? 0 : 1;
}
