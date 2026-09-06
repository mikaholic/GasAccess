#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/desorption_updater.hpp"
#include "gasaccess/distributed_desorption_updater.hpp"
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

using gasaccess::Atom;
using gasaccess::AtomView;
using gasaccess::AtomVoxelizer;
using gasaccess::DesorptionRepairMode;
using gasaccess::DesorptionUpdateResult;
using gasaccess::DesorptionUpdater;
using gasaccess::DistributedDesorptionUpdateResult;
using gasaccess::DistributedDesorptionUpdater;
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

GridSpec make_grid_spec()
{
    GridSpec grid_spec{};
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {32, 8, 8};
    grid_spec.periodic.x = true;
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

std::vector<VoxelCoord> make_plane(
    const GridSpec& grid_spec,
    std::int64_t z)
{
    std::vector<VoxelCoord> plane;
    plane.reserve(static_cast<std::size_t>(
        grid_spec.dimensions.x * grid_spec.dimensions.y));
    for (std::int64_t y = 0;
         y < static_cast<std::int64_t>(grid_spec.dimensions.y);
         ++y) {
        for (std::int64_t x = 0;
             x < static_cast<std::int64_t>(grid_spec.dimensions.x);
             ++x) {
            plane.push_back({x, y, z});
        }
    }
    return plane;
}

std::vector<VoxelCoord> make_all_voxels_except(
    const GridSpec& grid_spec,
    const std::vector<VoxelCoord>& excluded)
{
    std::vector<VoxelCoord> voxel_coords;
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
                if (std::find(excluded.begin(), excluded.end(), voxel_coord)
                    == excluded.end()) {
                    voxel_coords.push_back(voxel_coord);
                }
            }
        }
    }
    return voxel_coords;
}

void initialize_from_atoms(
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
    const DesorptionUpdateResult& serial_result)
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
    const DistributedDesorptionUpdateResult& incremental_result,
    const DistributedDesorptionUpdateResult& full_result)
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
                REQUIRE(incremental_grid.owned_blocker_count(voxel_coord)
                    == serial_grid.blocker_count(voxel_coord));
                REQUIRE(full_grid.owned_blocker_count(voxel_coord)
                    == serial_grid.blocker_count(voxel_coord));
                ++counted_states[static_cast<std::size_t>(serial_state)];

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

using ResultVerifier = std::function<void(
    std::size_t,
    const DistributedDesorptionUpdateResult&,
    const DesorptionUpdateResult&)>;

void run_sequence(
    const GridSpec& grid_spec,
    const std::vector<VoxelCoord>& initial_atom_voxels,
    const std::vector<std::vector<VoxelCoord>>& removal_events,
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
    initialize_from_atoms(
        serial_grid,
        incremental_grid,
        full_grid,
        initial_atom_voxels);

    DesorptionUpdater serial_updater(0.0);
    DistributedDesorptionUpdater incremental_updater(0.0);
    DistributedDesorptionUpdater full_updater(
        0.0,
        DesorptionRepairMode::FullReclassification);

    for (std::size_t event_index = 0;
         event_index < removal_events.size();
         ++event_index) {
        const auto removed_atoms = atoms_for_voxels(
            serial_grid,
            removal_events[event_index]);
        const AtomView removed_atom_view{
            removed_atoms.data(), removed_atoms.size()};
        const auto serial_result = serial_updater.apply_desorption(
            serial_grid,
            removed_atom_view);
        const auto incremental_result = incremental_updater.apply_desorption(
            incremental_grid,
            removed_atom_view);
        const auto full_result = full_updater.apply_desorption(
            full_grid,
            removed_atom_view);

        REQUIRE(incremental_result.global_blocker_count_changed_voxel_count
            == serial_result.blocker_count_changed_voxel_count);
        REQUIRE(full_result.global_blocker_count_changed_voxel_count
            == serial_result.blocker_count_changed_voxel_count);
        REQUIRE(incremental_result.global_newly_gas_count
            == serial_result.newly_gas_count);
        REQUIRE(full_result.global_newly_gas_count
            == serial_result.newly_gas_count);
        REQUIRE(incremental_result.geometry_changed()
            == serial_result.geometry_changed());
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

        std::array<std::uint64_t, 4> local_metrics{
            incremental_result.local_repair_opened_voxel_count,
            incremental_result.local_opening_visited_voxel_count,
            incremental_result.sent_frontier_entry_count,
            incremental_result.received_frontier_entry_count};
        std::array<std::uint64_t, 4> global_metrics{};
        MPI_Allreduce(
            local_metrics.data(),
            global_metrics.data(),
            static_cast<int>(global_metrics.size()),
            MPI_UINT64_T,
            MPI_SUM,
            MPI_COMM_WORLD);
        REQUIRE(global_metrics[0] == serial_result.repair_opened_voxel_count);
        REQUIRE(global_metrics[1]
            == serial_result.opening_visited_voxel_count);
        REQUIRE(global_metrics[2] == global_metrics[3]);
        REQUIRE(incremental_result.local_repair_opened_voxel_count
            == incremental_result.local_opening_visited_voxel_count);
        REQUIRE(incremental_result.used_distributed_opening_repair()
            == serial_result.used_opening_region_repair());

        const std::uint64_t local_participated =
            incremental_result.local_opening_visited_voxel_count == 0
            ? 0U
            : 1U;
        std::uint64_t participating_ranks = 0;
        MPI_Allreduce(
            &local_participated,
            &participating_ranks,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            MPI_COMM_WORLD);
        REQUIRE(incremental_result.opening_participating_rank_count
            == participating_ranks);

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

void test_collective_validation_and_rollback(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const auto block_width = std::int64_t{32}
        / static_cast<std::int64_t>(size);
    std::vector<VoxelCoord> occupied_voxels;
    occupied_voxels.reserve(static_cast<std::size_t>(size));
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
    const auto atoms = atoms_for_voxels(serial_grid, occupied_voxels);
    AtomVoxelizer(0.0).voxelize(serial_grid, {atoms.data(), atoms.size()});
    distributed_grid.voxelize_owned_atoms({atoms.data(), atoms.size()}, 0.0);
    ExteriorClassifier{}.classify(serial_grid);
    DistributedExteriorClassifier{}.classify(distributed_grid);

    const VoxelCoord absent{2, 2, 2};
    const auto rank_local_removal = rank == 0
        ? absent
        : occupied_voxels[static_cast<std::size_t>(rank)];
    const auto invalid_atoms = atoms_for_voxels(
        serial_grid,
        {rank_local_removal});
    bool threw = false;
    try {
        static_cast<void>(DistributedDesorptionUpdater(0.0).apply_desorption(
            distributed_grid,
            {invalid_atoms.data(), invalid_atoms.size()}));
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
    if (distributed_grid.owns(absent)) {
        REQUIRE(distributed_grid.owned_blocker_count(absent) == 0);
        REQUIRE(distributed_grid.gas_state(absent)
            == GasState::OutsideAccessible);
    }

    DistributedGasGrid unclassified_grid(
        grid_spec,
        make_decomposition_spec(grid_spec, rank, size));
    threw = false;
    try {
        static_cast<void>(DistributedDesorptionUpdater(0.0).apply_desorption(
            unclassified_grid,
            {nullptr, 0}));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    REQUIRE(threw);
}

void test_overlap_empty_and_final_removal(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const VoxelCoord occupied{1, 2, 2};
    run_sequence(
        grid_spec,
        {occupied, occupied},
        {{}, {occupied}, {occupied}},
        rank,
        size,
        [](std::size_t event_index,
           const DistributedDesorptionUpdateResult& distributed_result,
           const DesorptionUpdateResult&) {
            if (event_index == 0) {
                REQUIRE(
                    distributed_result.global_blocker_count_changed_voxel_count
                    == 0);
                REQUIRE(!distributed_result.geometry_changed());
            } else if (event_index == 1) {
                REQUIRE(
                    distributed_result.global_blocker_count_changed_voxel_count
                    == 1);
                REQUIRE(distributed_result.global_newly_gas_count == 0);
                REQUIRE(!distributed_result.used_distributed_opening_repair());
            } else {
                REQUIRE(distributed_result.global_newly_gas_count == 1);
                REQUIRE(distributed_result.opening_participating_rank_count == 1);
            }
        });
}

void test_closed_and_local_openings(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const VoxelCoord enclosed{9, 3, 3};
    const std::vector<VoxelCoord> enclosure{
        enclosed,
        {8, 3, 3},
        {10, 3, 3},
        {9, 2, 3},
        {9, 4, 3},
        {9, 3, 2},
        {9, 3, 4}};
    run_sequence(
        grid_spec,
        enclosure,
        {{enclosed}},
        rank,
        size,
        [](std::size_t,
           const DistributedDesorptionUpdateResult& distributed_result,
           const DesorptionUpdateResult& serial_result) {
            REQUIRE(distributed_result.global_newly_gas_count == 1);
            REQUIRE(serial_result.repair_opened_voxel_count == 0);
            REQUIRE(!distributed_result.used_distributed_opening_repair());
        });

    const VoxelCoord local_opening{1, 3, 3};
    run_sequence(
        grid_spec,
        {local_opening},
        {{local_opening}},
        rank,
        size,
        [](std::size_t,
           const DistributedDesorptionUpdateResult& distributed_result,
           const DesorptionUpdateResult& serial_result) {
            REQUIRE(serial_result.repair_opened_voxel_count == 1);
            REQUIRE(distributed_result.opening_participating_rank_count == 1);
        });
}

void test_rank_boundary_opening(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const auto block_width = std::int64_t{32}
        / static_cast<std::int64_t>(size);
    const VoxelCoord boundary_opening{
        size == 1 ? block_width / 2 : block_width,
        3,
        3};
    run_sequence(
        grid_spec,
        {boundary_opening},
        {{boundary_opening}},
        rank,
        size,
        [](std::size_t,
           const DistributedDesorptionUpdateResult& distributed_result,
           const DesorptionUpdateResult& serial_result) {
            REQUIRE(serial_result.repair_opened_voxel_count == 1);
            REQUIRE(distributed_result.opening_participating_rank_count == 1);
        });
}

void test_periodic_seam_opening(int rank, int size)
{
    auto grid_spec = make_grid_spec();
    grid_spec.reservoir_faces = {};
    const VoxelCoord closed_side{0, 3, 3};
    const VoxelCoord source_plug{31, 3, 3};
    grid_spec.explicit_source_voxels.push_back(source_plug);
    const auto initial_atoms = make_all_voxels_except(
        grid_spec,
        {closed_side});
    run_sequence(
        grid_spec,
        initial_atoms,
        {{source_plug}},
        rank,
        size,
        [size](std::size_t,
               const DistributedDesorptionUpdateResult& distributed_result,
               const DesorptionUpdateResult& serial_result) {
            REQUIRE(serial_result.repair_opened_voxel_count == 2);
            const auto expected_participants = size == 1 ? 1U : 2U;
            REQUIRE(distributed_result.opening_participating_rank_count
                == expected_participants);
            std::uint64_t global_sent = 0;
            MPI_Allreduce(
                &distributed_result.sent_frontier_entry_count,
                &global_sent,
                1,
                MPI_UINT64_T,
                MPI_SUM,
                MPI_COMM_WORLD);
            if (size > 1) {
                REQUIRE(global_sent > 0);
            }
        });
}

void test_large_opening_reaches_every_rank(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    const VoxelCoord plug{15, 3, 6};
    const auto initial_atoms = make_plane(grid_spec, plug.z);
    constexpr std::uint64_t closed_cavity_voxels = 32U * 8U * 6U;
    constexpr std::uint64_t opened_voxels = closed_cavity_voxels + 1U;
    constexpr std::uint64_t total_voxels = 32U * 8U * 8U;
    static_assert(opened_voxels > total_voxels / 2U);

    run_sequence(
        grid_spec,
        initial_atoms,
        {{plug}},
        rank,
        size,
        [size](std::size_t,
               const DistributedDesorptionUpdateResult& distributed_result,
               const DesorptionUpdateResult& serial_result) {
            REQUIRE(serial_result.repair_opened_voxel_count == opened_voxels);
            REQUIRE(distributed_result.opening_participating_rank_count
                == static_cast<std::uint64_t>(size));
            REQUIRE(distributed_result.local_opening_visited_voxel_count > 0);
            REQUIRE(distributed_result.repair_communication_round_count > 0);
            std::uint64_t global_sent = 0;
            MPI_Allreduce(
                &distributed_result.sent_frontier_entry_count,
                &global_sent,
                1,
                MPI_UINT64_T,
                MPI_SUM,
                MPI_COMM_WORLD);
            if (size > 1) {
                REQUIRE(global_sent > 0);
            }
        });
}

void test_deterministic_desorption_sequence(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    auto initial_atoms = make_plane(grid_spec, 4);
    std::vector<VoxelCoord> removals = initial_atoms;
    std::mt19937_64 random_generator(872341U);
    std::shuffle(removals.begin(), removals.end(), random_generator);
    removals.resize(16);

    std::vector<std::vector<VoxelCoord>> events;
    events.reserve(removals.size());
    for (const auto& voxel_coord : removals) {
        events.push_back({voxel_coord});
    }
    run_sequence(grid_spec, initial_atoms, events, rank, size);
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
        {"collective validation and rollback", [&]() {
             test_collective_validation_and_rollback(rank, size);
         }},
        {"overlap, empty, and final removal", [&]() {
             test_overlap_empty_and_final_removal(rank, size);
         }},
        {"closed and local openings", [&]() {
             test_closed_and_local_openings(rank, size);
         }},
        {"rank-boundary opening", [&]() {
             test_rank_boundary_opening(rank, size);
         }},
        {"periodic-seam opening", [&]() {
             test_periodic_seam_opening(rank, size);
         }},
        {"large opening reaches every rank", [&]() {
             test_large_opening_reaches_every_rank(rank, size);
         }},
        {"deterministic desorption sequence", [&]() {
             test_deterministic_desorption_sequence(rank, size);
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
                  << " distributed desorption test groups passed on "
                  << size << " rank(s)\n";
    }
    MPI_Finalize();
    return global_failure_count == 0 ? 0 : 1;
}
