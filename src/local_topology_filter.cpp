#include "gasaccess/local_topology_filter.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace gasaccess {
namespace {

constexpr std::size_t local_voxel_capacity = 27;

std::uint64_t axis_distance(
    std::int64_t lhs,
    std::int64_t rhs,
    std::uint64_t dimension,
    bool periodic)
{
    const auto lhs_index = static_cast<std::uint64_t>(lhs);
    const auto rhs_index = static_cast<std::uint64_t>(rhs);
    const auto direct_distance = lhs_index > rhs_index
        ? lhs_index - rhs_index
        : rhs_index - lhs_index;
    if (!periodic) {
        return direct_distance;
    }
    return std::min(direct_distance, dimension - direct_distance);
}

bool contains_id(
    const std::array<VoxelId, local_voxel_capacity>& voxel_ids,
    std::size_t count,
    VoxelId voxel_id)
{
    return std::find(
        voxel_ids.begin(),
        voxel_ids.begin() + static_cast<std::ptrdiff_t>(count),
        voxel_id) != voxel_ids.begin() + static_cast<std::ptrdiff_t>(count);
}

void validate_removed_voxels(
    const GasGrid& gas_grid,
    RemovedVoxelView removed_voxel_view)
{
    if (removed_voxel_view.count != 0
        && removed_voxel_view.removed_voxels == nullptr) {
        throw std::invalid_argument("removed-voxel view has a null pointer");
    }

    for (std::size_t index = 0; index < removed_voxel_view.count; ++index) {
        const auto& removed_voxel = removed_voxel_view.removed_voxels[index];
        if (gas_grid.gas_state(removed_voxel.voxel_id) != GasState::Solid) {
            throw std::invalid_argument("removed voxel is not solid in the updated grid");
        }
        if (removed_voxel.previous_state != GasState::OutsideAccessible
            && removed_voxel.previous_state != GasState::ClosedVoid) {
            throw std::invalid_argument(
                "removed voxel was not classified as empty before adsorption");
        }
    }
}

}  // namespace

bool TopologyCheckResult::is_safe() const noexcept
{
    return decision == TopologyDecision::Safe;
}

TopologyCheckResult LocalTopologyFilter::evaluate(
    const GasGrid& gas_grid,
    RemovedVoxelView removed_voxel_view) const
{
    validate_removed_voxels(gas_grid, removed_voxel_view);
    if (removed_voxel_view.count == 0) {
        return {TopologyDecision::Safe, 0, 0};
    }
    if (removed_voxel_view.count != 1) {
        return {TopologyDecision::RequiresConnectivityRepair, 0, 0};
    }

    const auto& removed_voxel = removed_voxel_view.removed_voxels[0];
    if (removed_voxel.previous_state == GasState::ClosedVoid) {
        return {TopologyDecision::Safe, 0, 0};
    }
    if (gas_grid.is_reservoir_source(removed_voxel.voxel_id)) {
        return {TopologyDecision::RequiresConnectivityRepair, 0, 0};
    }

    std::array<VoxelId, 6> accessible_neighbors{};
    std::size_t accessible_neighbor_count = 0;
    const auto neighbor_list = gas_grid.neighbors(removed_voxel.voxel_id);
    for (std::size_t index = 0; index < neighbor_list.count; ++index) {
        const auto neighbor_id = neighbor_list.ids[index];
        if (gas_grid.gas_state(neighbor_id) == GasState::OutsideAccessible) {
            accessible_neighbors[accessible_neighbor_count] = neighbor_id;
            ++accessible_neighbor_count;
        }
    }

    if (accessible_neighbor_count <= 1) {
        return {
            TopologyDecision::Safe,
            accessible_neighbor_count,
            accessible_neighbor_count
        };
    }

    std::array<VoxelId, local_voxel_capacity> visited_voxels{};
    visited_voxels[0] = accessible_neighbors[0];
    std::size_t visited_count = 1;
    std::size_t frontier_index = 0;
    while (frontier_index < visited_count) {
        const auto current_voxel_id = visited_voxels[frontier_index];
        ++frontier_index;
        const auto current_neighbors = gas_grid.neighbors(current_voxel_id);
        for (std::size_t index = 0; index < current_neighbors.count; ++index) {
            const auto candidate_id = current_neighbors.ids[index];
            if (gas_grid.gas_state(candidate_id) != GasState::OutsideAccessible
                || !is_in_local_neighborhood(
                    gas_grid,
                    removed_voxel.voxel_id,
                    candidate_id)
                || contains_id(visited_voxels, visited_count, candidate_id)) {
                continue;
            }
            if (visited_count >= visited_voxels.size()) {
                throw std::logic_error("local topology search exceeded 3x3x3 capacity");
            }
            visited_voxels[visited_count] = candidate_id;
            ++visited_count;
        }
    }

    for (std::size_t index = 1; index < accessible_neighbor_count; ++index) {
        if (!contains_id(
            visited_voxels,
            visited_count,
            accessible_neighbors[index])) {
            return {
                TopologyDecision::RequiresConnectivityRepair,
                accessible_neighbor_count,
                visited_count
            };
        }
    }
    return {
        TopologyDecision::Safe,
        accessible_neighbor_count,
        visited_count
    };
}

bool LocalTopologyFilter::is_in_local_neighborhood(
    const GasGrid& gas_grid,
    VoxelId center_voxel_id,
    VoxelId candidate_voxel_id)
{
    const auto center = gas_grid.voxel_coord(center_voxel_id);
    const auto candidate = gas_grid.voxel_coord(candidate_voxel_id);
    const auto& grid_spec = gas_grid.grid_spec();
    return axis_distance(
        center.x,
        candidate.x,
        grid_spec.dimensions.x,
        grid_spec.periodic.x) <= 1
        && axis_distance(
            center.y,
            candidate.y,
            grid_spec.dimensions.y,
            grid_spec.periodic.y) <= 1
        && axis_distance(
            center.z,
            candidate.z,
            grid_spec.dimensions.z,
            grid_spec.periodic.z) <= 1;
}

}  // namespace gasaccess
