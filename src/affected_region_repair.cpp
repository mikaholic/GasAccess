#include "gasaccess/affected_region_repair.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace gasaccess {
namespace {

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

void increment_checked(VoxelId& value)
{
    if (value == std::numeric_limits<VoxelId>::max()) {
        throw std::overflow_error("affected-region visit count overflow");
    }
    ++value;
}

}  // namespace

AffectedRegionRepairResult AffectedRegionRepair::repair(
    GasGrid& gas_grid,
    RemovedVoxelView removed_voxel_view)
{
    validate_removed_voxels(gas_grid, removed_voxel_view);
    prepare_workspace(gas_grid);
    begin_repair_epoch();
    seed_voxels_.clear();

    for (std::size_t index = 0; index < removed_voxel_view.count; ++index) {
        const auto& removed_voxel = removed_voxel_view.removed_voxels[index];
        if (removed_voxel.previous_state != GasState::OutsideAccessible) {
            continue;
        }
        const auto neighbors = gas_grid.neighbors(removed_voxel.voxel_id);
        for (std::size_t neighbor_index = 0;
             neighbor_index < neighbors.count;
             ++neighbor_index) {
            const auto neighbor_id = neighbors.ids[neighbor_index];
            if (gas_grid.gas_state(neighbor_id) == GasState::OutsideAccessible) {
                seed_voxels_.push_back(neighbor_id);
            }
        }
    }

    AffectedRegionRepairResult result{};
    for (const auto seed_voxel_id : seed_voxels_) {
        if (gas_grid.gas_state(seed_voxel_id) != GasState::OutsideAccessible
            || is_confirmed_outside(seed_voxel_id)) {
            continue;
        }
        if (gas_grid.is_reservoir_source(seed_voxel_id)) {
            confirmed_outside_epochs_[static_cast<std::size_t>(seed_voxel_id)] =
                repair_epoch_;
            continue;
        }

        begin_search_epoch();
        frontier_.clear();
        frontier_.push_back(seed_voxel_id);
        search_epochs_[static_cast<std::size_t>(seed_voxel_id)] = search_epoch_;

        bool reaches_outside_source = false;
        std::size_t frontier_index = 0;
        while (frontier_index < frontier_.size() && !reaches_outside_source) {
            const auto current_voxel_id = frontier_[frontier_index];
            ++frontier_index;
            increment_checked(result.visited_voxel_count);
            if (gas_grid.is_reservoir_source(current_voxel_id)) {
                reaches_outside_source = true;
                break;
            }

            const auto neighbors = gas_grid.neighbors(current_voxel_id);
            for (std::size_t neighbor_index = 0;
                 neighbor_index < neighbors.count;
                 ++neighbor_index) {
                const auto neighbor_id = neighbors.ids[neighbor_index];
                if (gas_grid.gas_state(neighbor_id)
                    != GasState::OutsideAccessible) {
                    continue;
                }
                if (is_confirmed_outside(neighbor_id)) {
                    reaches_outside_source = true;
                    break;
                }
                const auto neighbor_index_value =
                    static_cast<std::size_t>(neighbor_id);
                if (search_epochs_[neighbor_index_value] == search_epoch_) {
                    continue;
                }
                search_epochs_[neighbor_index_value] = search_epoch_;
                frontier_.push_back(neighbor_id);
            }
        }

        if (reaches_outside_source) {
            confirm_region_outside();
            continue;
        }

        for (const auto voxel_id : frontier_) {
            gas_grid.set_gas_state(voxel_id, GasState::ClosedVoid);
            result.newly_closed_voxel_ids.push_back(voxel_id);
        }
    }

    std::sort(
        result.newly_closed_voxel_ids.begin(),
        result.newly_closed_voxel_ids.end());
    return result;
}

void AffectedRegionRepair::prepare_workspace(const GasGrid& gas_grid)
{
    const auto voxel_count = static_cast<std::size_t>(gas_grid.voxel_count());
    if (search_epochs_.size() == voxel_count
        && confirmed_outside_epochs_.size() == voxel_count) {
        return;
    }
    search_epochs_.assign(voxel_count, 0);
    confirmed_outside_epochs_.assign(voxel_count, 0);
    search_epoch_ = 0;
    repair_epoch_ = 0;
}

void AffectedRegionRepair::begin_repair_epoch()
{
    ++repair_epoch_;
    if (repair_epoch_ != 0) {
        return;
    }
    std::fill(confirmed_outside_epochs_.begin(), confirmed_outside_epochs_.end(), 0);
    repair_epoch_ = 1;
}

void AffectedRegionRepair::begin_search_epoch()
{
    ++search_epoch_;
    if (search_epoch_ != 0) {
        return;
    }
    std::fill(search_epochs_.begin(), search_epochs_.end(), 0);
    search_epoch_ = 1;
}

bool AffectedRegionRepair::is_confirmed_outside(VoxelId voxel_id) const noexcept
{
    return confirmed_outside_epochs_[static_cast<std::size_t>(voxel_id)]
        == repair_epoch_;
}

void AffectedRegionRepair::confirm_region_outside() noexcept
{
    for (const auto voxel_id : frontier_) {
        confirmed_outside_epochs_[static_cast<std::size_t>(voxel_id)] = repair_epoch_;
    }
}

}  // namespace gasaccess
