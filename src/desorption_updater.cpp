#include "gasaccess/desorption_updater.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace gasaccess {
namespace {

ClassificationSummary current_classification(const GasGrid& gas_grid)
{
    if (gas_grid.gas_state_count(GasState::Unclassified) != 0) {
        throw std::invalid_argument(
            "desorption update requires a fully classified gas grid");
    }
    return {
        gas_grid.gas_state_count(GasState::Solid),
        gas_grid.gas_state_count(GasState::OutsideAccessible),
        gas_grid.gas_state_count(GasState::ClosedVoid)
    };
}

}  // namespace

bool DesorptionUpdateResult::geometry_changed() const noexcept
{
    return newly_gas_count != 0;
}

bool DesorptionUpdateResult::used_full_reclassification() const noexcept
{
    return full_reclassification_performed;
}

bool DesorptionUpdateResult::used_opening_region_repair() const noexcept
{
    return opening_region_repair_performed;
}

DesorptionUpdater::DesorptionUpdater(
    double precursor_radius,
    DesorptionRepairMode repair_mode)
    : atom_voxelizer_(precursor_radius),
      repair_mode_(repair_mode)
{
    if (repair_mode_ != DesorptionRepairMode::OpeningRegion
        && repair_mode_ != DesorptionRepairMode::FullReclassification) {
        throw std::invalid_argument("invalid desorption repair mode");
    }
}

double DesorptionUpdater::precursor_radius() const noexcept
{
    return atom_voxelizer_.precursor_radius();
}

DesorptionRepairMode DesorptionUpdater::repair_mode() const noexcept
{
    return repair_mode_;
}

DesorptionUpdateResult DesorptionUpdater::apply_desorption(
    GasGrid& gas_grid,
    AtomView removed_atoms) const
{
    DesorptionUpdateResult result{};
    result.classification = current_classification(gas_grid);

    std::vector<GasState> previous_states;
    if (repair_mode_ == DesorptionRepairMode::FullReclassification) {
        previous_states.reserve(static_cast<std::size_t>(gas_grid.voxel_count()));
        for (VoxelId voxel_id = 0;
             voxel_id < gas_grid.voxel_count();
             ++voxel_id) {
            previous_states.push_back(gas_grid.gas_state(voxel_id));
        }
    }

    const auto occupancy_result = atom_voxelizer_.apply_atom_changes(
        gas_grid,
        {{nullptr, 0}, removed_atoms},
        occupancy_changes_);
    result.blocker_count_changed_voxel_count =
        occupancy_result.blocker_count_changed_voxel_count;
    result.newly_gas_count = occupancy_result.newly_gas_count;
    if (!result.geometry_changed()) {
        return result;
    }

    newly_gas_voxel_ids_.clear();
    newly_gas_voxel_ids_.reserve(
        static_cast<std::size_t>(result.newly_gas_count));
    for (const auto& occupancy_change : occupancy_changes_) {
        if (occupancy_change.previous_blocker_count != 0
            && occupancy_change.blocker_count == 0) {
            newly_gas_voxel_ids_.push_back(occupancy_change.voxel_id);
        }
    }
    if (newly_gas_voxel_ids_.size()
        != static_cast<std::size_t>(result.newly_gas_count)) {
        throw std::logic_error("voxelizer gas-count result is inconsistent");
    }

    if (repair_mode_ == DesorptionRepairMode::FullReclassification) {
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

    for (const auto voxel_id : newly_gas_voxel_ids_) {
        gas_grid.set_gas_state(voxel_id, GasState::ClosedVoid);
    }
    const auto repair_result = opening_region_repair_.repair(
        gas_grid,
        {newly_gas_voxel_ids_.data(), newly_gas_voxel_ids_.size()});
    result.opening_visited_voxel_count = repair_result.visited_voxel_count;
    result.repair_opened_voxel_count = static_cast<VoxelId>(
        repair_result.newly_opened_voxel_ids.size());
    result.opening_region_repair_performed =
        result.opening_visited_voxel_count != 0;

    result.changed_voxel_ids = newly_gas_voxel_ids_;
    result.changed_voxel_ids.insert(
        result.changed_voxel_ids.end(),
        repair_result.newly_opened_voxel_ids.begin(),
        repair_result.newly_opened_voxel_ids.end());
    std::sort(result.changed_voxel_ids.begin(), result.changed_voxel_ids.end());
    result.changed_voxel_ids.erase(
        std::unique(
            result.changed_voxel_ids.begin(),
            result.changed_voxel_ids.end()),
        result.changed_voxel_ids.end());
    result.classification = current_classification(gas_grid);
    return result;
}

}  // namespace gasaccess
