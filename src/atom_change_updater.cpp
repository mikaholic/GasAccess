#include "gasaccess/atom_change_updater.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace gasaccess {
namespace {

ClassificationSummary current_classification(const GasGrid& gas_grid)
{
    if (gas_grid.gas_state_count(GasState::Unclassified) != 0) {
        throw std::invalid_argument(
            "atom-change update requires a fully classified gas grid");
    }
    return {
        gas_grid.gas_state_count(GasState::Solid),
        gas_grid.gas_state_count(GasState::OutsideAccessible),
        gas_grid.gas_state_count(GasState::ClosedVoid)};
}

AccessibilityRepairKind incremental_repair_kind(
    VoxelId newly_solid_count,
    VoxelId newly_gas_count) noexcept
{
    if (newly_solid_count != 0 && newly_gas_count != 0) {
        return AccessibilityRepairKind::Mixed;
    }
    if (newly_solid_count != 0) {
        return AccessibilityRepairKind::Closing;
    }
    if (newly_gas_count != 0) {
        return AccessibilityRepairKind::Opening;
    }
    return AccessibilityRepairKind::None;
}

}  // namespace

bool AtomChangeUpdateResult::geometry_changed() const noexcept
{
    return newly_solid_count != 0 || newly_gas_count != 0;
}

bool AtomChangeUpdateResult::used_full_reclassification() const noexcept
{
    return full_reclassification_performed;
}

bool AtomChangeUpdateResult::used_closing_repair() const noexcept
{
    return closing_repair_performed;
}

bool AtomChangeUpdateResult::used_opening_repair() const noexcept
{
    return opening_repair_performed;
}

AtomChangeUpdater::AtomChangeUpdater(
    double precursor_radius,
    AtomChangeRepairMode repair_mode)
    : atom_voxelizer_(precursor_radius),
      repair_mode_(repair_mode)
{
    if (repair_mode_ != AtomChangeRepairMode::Incremental
        && repair_mode_ != AtomChangeRepairMode::FullReclassification) {
        throw std::invalid_argument("invalid atom-change repair mode");
    }
}

double AtomChangeUpdater::precursor_radius() const noexcept
{
    return atom_voxelizer_.precursor_radius();
}

AtomChangeRepairMode AtomChangeUpdater::repair_mode() const noexcept
{
    return repair_mode_;
}

AtomChangeUpdateResult AtomChangeUpdater::apply_atom_changes(
    GasGrid& gas_grid,
    const AtomChangeBatch& atom_changes) const
{
    AtomChangeUpdateResult result{};
    result.classification = current_classification(gas_grid);

    if (repair_mode_ == AtomChangeRepairMode::FullReclassification) {
        previous_states_.clear();
        previous_states_.reserve(
            static_cast<std::size_t>(gas_grid.voxel_count()));
        for (VoxelId voxel_id = 0;
             voxel_id < gas_grid.voxel_count();
             ++voxel_id) {
            previous_states_.push_back(gas_grid.gas_state(voxel_id));
        }
    }

    const auto occupancy_result = atom_voxelizer_.apply_atom_changes(
        gas_grid,
        atom_changes,
        occupancy_changes_);
    result.blocker_count_changed_voxel_count =
        occupancy_result.blocker_count_changed_voxel_count;
    result.newly_solid_count = occupancy_result.newly_solid_count;
    result.newly_gas_count = occupancy_result.newly_gas_count;
    if (!result.geometry_changed()) {
        return result;
    }

    if (repair_mode_ == AtomChangeRepairMode::FullReclassification) {
        result.repair_kind = AccessibilityRepairKind::FullReclassification;
        result.full_reclassification_performed = true;
        result.classification = ExteriorClassifier{}.classify(gas_grid);
        for (VoxelId voxel_id = 0;
             voxel_id < gas_grid.voxel_count();
             ++voxel_id) {
            if (previous_states_[static_cast<std::size_t>(voxel_id)]
                != gas_grid.gas_state(voxel_id)) {
                result.changed_voxel_ids.push_back(voxel_id);
            }
        }
        return result;
    }

    result.repair_kind = incremental_repair_kind(
        result.newly_solid_count,
        result.newly_gas_count);
    newly_solid_voxels_.clear();
    newly_solid_voxels_.reserve(
        static_cast<std::size_t>(result.newly_solid_count));
    newly_gas_voxel_ids_.clear();
    newly_gas_voxel_ids_.reserve(
        static_cast<std::size_t>(result.newly_gas_count));
    for (const auto& change : occupancy_changes_) {
        if (change.previous_blocker_count == 0
            && change.blocker_count != 0) {
            newly_solid_voxels_.push_back({
                change.voxel_id,
                change.previous_state});
            result.changed_voxel_ids.push_back(change.voxel_id);
        } else if (change.previous_blocker_count != 0
                   && change.blocker_count == 0) {
            newly_gas_voxel_ids_.push_back(change.voxel_id);
            result.changed_voxel_ids.push_back(change.voxel_id);
        }
    }
    if (newly_solid_voxels_.size()
            != static_cast<std::size_t>(result.newly_solid_count)
        || newly_gas_voxel_ids_.size()
            != static_cast<std::size_t>(result.newly_gas_count)) {
        throw std::logic_error("atom-change transition counts are inconsistent");
    }

    if (!newly_solid_voxels_.empty()) {
        auto closing_result = closing_repair_.repair(
            gas_grid,
            {newly_solid_voxels_.data(), newly_solid_voxels_.size()});
        result.closing_visited_voxel_count =
            closing_result.visited_voxel_count;
        result.repair_closed_voxel_count = static_cast<VoxelId>(
            closing_result.newly_closed_voxel_ids.size());
        result.closing_repair_performed =
            result.closing_visited_voxel_count != 0;
        closing_changed_voxel_ids_ = std::move(
            closing_result.newly_closed_voxel_ids);
    } else {
        closing_changed_voxel_ids_.clear();
    }

    if (!newly_gas_voxel_ids_.empty()) {
        for (const auto voxel_id : newly_gas_voxel_ids_) {
            gas_grid.set_gas_state(voxel_id, GasState::ClosedVoid);
        }
        const auto opening_result = opening_repair_.repair(
            gas_grid,
            {newly_gas_voxel_ids_.data(), newly_gas_voxel_ids_.size()});
        result.opening_visited_voxel_count =
            opening_result.visited_voxel_count;
        result.repair_opened_voxel_count = static_cast<VoxelId>(
            opening_result.newly_opened_voxel_ids.size());
        result.opening_repair_performed =
            result.opening_visited_voxel_count != 0;
        for (const auto voxel_id : opening_result.newly_opened_voxel_ids) {
            if (!std::binary_search(
                    closing_changed_voxel_ids_.begin(),
                    closing_changed_voxel_ids_.end(),
                    voxel_id)) {
                result.changed_voxel_ids.push_back(voxel_id);
            }
        }
    }

    for (const auto voxel_id : closing_changed_voxel_ids_) {
        if (gas_grid.gas_state(voxel_id) == GasState::ClosedVoid) {
            result.changed_voxel_ids.push_back(voxel_id);
        }
    }

    std::sort(result.changed_voxel_ids.begin(), result.changed_voxel_ids.end());
    result.changed_voxel_ids.erase(
        std::unique(
            result.changed_voxel_ids.begin(),
            result.changed_voxel_ids.end()),
        result.changed_voxel_ids.end());
    result.classification = current_classification(gas_grid);
    return result;
}

AtomChangeUpdateResult AtomChangeUpdater::apply_deposition(
    GasGrid& gas_grid,
    AtomView added_atoms) const
{
    return apply_atom_changes(gas_grid, {added_atoms, {nullptr, 0}});
}

AtomChangeUpdateResult AtomChangeUpdater::apply_desorption(
    GasGrid& gas_grid,
    AtomView removed_atoms) const
{
    return apply_atom_changes(gas_grid, {{nullptr, 0}, removed_atoms});
}

}  // namespace gasaccess
