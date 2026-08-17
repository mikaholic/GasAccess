#include "gasaccess/deposition_updater.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace gasaccess {
namespace {

ClassificationSummary current_classification(const GasGrid& gas_grid)
{
    if (gas_grid.gas_state_count(GasState::Unclassified) != 0) {
        throw std::invalid_argument(
            "deposition update requires a fully classified gas grid");
    }
    return {
        gas_grid.gas_state_count(GasState::Solid),
        gas_grid.gas_state_count(GasState::OutsideAccessible),
        gas_grid.gas_state_count(GasState::ClosedVoid)
    };
}

}  // namespace

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
    result.classification = current_classification(gas_grid);

    std::vector<GasState> previous_states;
    if (repair_mode_ == ConnectivityRepairMode::FullReclassification) {
        previous_states.reserve(static_cast<std::size_t>(gas_grid.voxel_count()));
        for (VoxelId voxel_id = 0;
             voxel_id < gas_grid.voxel_count();
             ++voxel_id) {
            previous_states.push_back(gas_grid.gas_state(voxel_id));
        }
    }

    result.newly_solid_count = atom_voxelizer_.voxelize(
        gas_grid,
        deposited_atoms,
        removed_voxels_);
    if (!result.geometry_changed()) {
        return result;
    }

    if (removed_voxels_.size()
        != static_cast<std::size_t>(result.newly_solid_count)) {
        throw std::logic_error("voxelizer solid-count result is inconsistent");
    }

    const auto append_removed_voxels = [&]() {
        result.changed_voxel_ids.reserve(removed_voxels_.size());
        for (const auto& removed_voxel : removed_voxels_) {
            result.changed_voxel_ids.push_back(removed_voxel.voxel_id);
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
        {removed_voxels_.data(), removed_voxels_.size()});
    if (topology_result.is_safe()) {
        append_removed_voxels();
        result.classification = current_classification(gas_grid);
        return result;
    }

    const auto repair_result = affected_region_repair_.repair(
        gas_grid,
        {removed_voxels_.data(), removed_voxels_.size()});
    result.affected_region_repair_performed = true;
    result.repair_visited_voxel_count = repair_result.visited_voxel_count;
    result.repair_closed_voxel_count = static_cast<VoxelId>(
        repair_result.newly_closed_voxel_ids.size());
    append_removed_voxels();
    result.changed_voxel_ids.insert(
        result.changed_voxel_ids.end(),
        repair_result.newly_closed_voxel_ids.begin(),
        repair_result.newly_closed_voxel_ids.end());
    std::sort(result.changed_voxel_ids.begin(), result.changed_voxel_ids.end());
    result.classification = current_classification(gas_grid);
    return result;
}

}  // namespace gasaccess
