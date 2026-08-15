#include "gasaccess/deposition_updater.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace gasaccess {

bool DepositionUpdateResult::geometry_changed() const noexcept
{
    return newly_solid_count != 0;
}

bool DepositionUpdateResult::used_full_reclassification() const noexcept
{
    return full_reclassification_performed;
}

bool DepositionUpdateResult::used_affected_region_repair() const noexcept
{
    return affected_region_repair_performed;
}

DepositionUpdater::DepositionUpdater(
    double precursor_radius,
    ConnectivityRepairMode repair_mode)
    : atom_voxelizer_(precursor_radius),
      repair_mode_(repair_mode)
{
    if (repair_mode_ != ConnectivityRepairMode::AffectedRegion
        && repair_mode_ != ConnectivityRepairMode::FullReclassification) {
        throw std::invalid_argument("invalid connectivity repair mode");
    }
}

double DepositionUpdater::precursor_radius() const noexcept
{
    return atom_voxelizer_.precursor_radius();
}

ConnectivityRepairMode DepositionUpdater::repair_mode() const noexcept
{
    return repair_mode_;
}

DepositionUpdateResult DepositionUpdater::apply_deposition(
    GasGrid& gas_grid,
    AtomView deposited_atoms) const
{
    DepositionUpdateResult result{};
    std::vector<GasState> previous_states;
    previous_states.reserve(static_cast<std::size_t>(gas_grid.voxel_count()));

    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        const auto gas_state = gas_grid.gas_state(voxel_id);
        previous_states.push_back(gas_state);
        switch (gas_state) {
        case GasState::Solid:
            ++result.classification.solid_count;
            break;
        case GasState::OutsideAccessible:
            ++result.classification.outside_accessible_count;
            break;
        case GasState::ClosedVoid:
            ++result.classification.closed_void_count;
            break;
        case GasState::Unclassified:
            throw std::invalid_argument(
                "deposition update requires a fully classified gas grid");
        default:
            throw std::invalid_argument("gas grid contains an invalid state value");
        }
    }

    result.newly_solid_count = atom_voxelizer_.voxelize(gas_grid, deposited_atoms);
    if (!result.geometry_changed()) {
        return result;
    }

    std::vector<RemovedVoxel> removed_voxels;
    removed_voxels.reserve(static_cast<std::size_t>(result.newly_solid_count));
    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        const auto previous_state = previous_states[static_cast<std::size_t>(voxel_id)];
        if (previous_state != GasState::Solid
            && gas_grid.gas_state(voxel_id) == GasState::Solid) {
            removed_voxels.push_back({voxel_id, previous_state});
        }
    }
    if (removed_voxels.size()
        != static_cast<std::size_t>(result.newly_solid_count)) {
        throw std::logic_error("voxelizer solid-count result is inconsistent");
    }

    const auto update_removed_voxel_counts = [&]() {
        result.changed_voxel_ids.reserve(removed_voxels.size());
        for (const auto& removed_voxel : removed_voxels) {
            result.changed_voxel_ids.push_back(removed_voxel.voxel_id);
            ++result.classification.solid_count;
            if (removed_voxel.previous_state == GasState::OutsideAccessible) {
                --result.classification.outside_accessible_count;
            } else {
                --result.classification.closed_void_count;
            }
        }
    };

    if (repair_mode_ == ConnectivityRepairMode::FullReclassification) {
        result.full_reclassification_performed = true;
        result.classification = ExteriorClassifier{}.classify(gas_grid);
        for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
            if (previous_states[static_cast<std::size_t>(voxel_id)]
                != gas_grid.gas_state(voxel_id)) {
                result.changed_voxel_ids.push_back(voxel_id);
            }
        }
        return result;
    }

    const auto topology_result = local_topology_filter_.evaluate(
        gas_grid,
        {removed_voxels.data(), removed_voxels.size()});
    if (topology_result.is_safe()) {
        update_removed_voxel_counts();
        return result;
    }

    const auto repair_result = affected_region_repair_.repair(
        gas_grid,
        {removed_voxels.data(), removed_voxels.size()});
    result.affected_region_repair_performed = true;
    result.repair_visited_voxel_count = repair_result.visited_voxel_count;
    result.repair_closed_voxel_count = static_cast<VoxelId>(
        repair_result.newly_closed_voxel_ids.size());
    update_removed_voxel_counts();
    result.classification.outside_accessible_count -=
        result.repair_closed_voxel_count;
    result.classification.closed_void_count += result.repair_closed_voxel_count;
    result.changed_voxel_ids.insert(
        result.changed_voxel_ids.end(),
        repair_result.newly_closed_voxel_ids.begin(),
        repair_result.newly_closed_voxel_ids.end());
    std::sort(result.changed_voxel_ids.begin(), result.changed_voxel_ids.end());
    return result;
}

}  // namespace gasaccess
