#include "gasaccess/opening_region_repair.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace gasaccess {
namespace {

void validate_newly_gas_voxels(
    const GasGrid& gas_grid,
    NewlyGasVoxelView newly_gas_voxel_view)
{
    if (newly_gas_voxel_view.count != 0
        && newly_gas_voxel_view.voxel_ids == nullptr) {
        throw std::invalid_argument("newly-gas voxel view has a null pointer");
    }
    if (gas_grid.gas_state_count(GasState::Unclassified) != 0) {
        throw std::invalid_argument(
            "opening repair requires a fully classified gas grid");
    }
    for (std::size_t index = 0; index < newly_gas_voxel_view.count; ++index) {
        const auto voxel_id = newly_gas_voxel_view.voxel_ids[index];
        if (gas_grid.gas_state(voxel_id) != GasState::ClosedVoid
            || gas_grid.blocker_count(voxel_id) != 0) {
            throw std::invalid_argument(
                "newly gas voxel is not an unblocked ClosedVoid voxel");
        }
    }
}

void increment_checked(VoxelId& value)
{
    if (value == std::numeric_limits<VoxelId>::max()) {
        throw std::overflow_error("opening-repair visit count overflow");
    }
    ++value;
}

}  // namespace

OpeningRegionRepairResult OpeningRegionRepair::repair(
    GasGrid& gas_grid,
    NewlyGasVoxelView newly_gas_voxel_view)
{
    validate_newly_gas_voxels(gas_grid, newly_gas_voxel_view);
    prepare_workspace(gas_grid);
    begin_repair_epoch();
    seed_voxels_.clear();

    for (std::size_t index = 0; index < newly_gas_voxel_view.count; ++index) {
        const auto voxel_id = newly_gas_voxel_view.voxel_ids[index];
        bool is_seed = gas_grid.is_reservoir_source(voxel_id);
        if (!is_seed) {
            const auto neighbors = gas_grid.neighbors(voxel_id);
            for (std::size_t neighbor_index = 0;
                 neighbor_index < neighbors.count;
                 ++neighbor_index) {
                if (gas_grid.gas_state(neighbors.ids[neighbor_index])
                    == GasState::OutsideAccessible) {
                    is_seed = true;
                    break;
                }
            }
        }
        if (is_seed) {
            seed_voxels_.push_back(voxel_id);
        }
    }

    OpeningRegionRepairResult result{};
    frontier_.clear();
    for (const auto seed_voxel_id : seed_voxels_) {
        if (was_visited(seed_voxel_id)) {
            continue;
        }
        mark_visited(seed_voxel_id);
        frontier_.push_back(seed_voxel_id);
    }

    std::size_t frontier_index = 0;
    while (frontier_index < frontier_.size()) {
        const auto voxel_id = frontier_[frontier_index];
        ++frontier_index;
        increment_checked(result.visited_voxel_count);
        gas_grid.set_gas_state(voxel_id, GasState::OutsideAccessible);
        result.newly_opened_voxel_ids.push_back(voxel_id);

        const auto neighbors = gas_grid.neighbors(voxel_id);
        for (std::size_t neighbor_index = 0;
             neighbor_index < neighbors.count;
             ++neighbor_index) {
            const auto neighbor_id = neighbors.ids[neighbor_index];
            if (gas_grid.gas_state(neighbor_id) != GasState::ClosedVoid
                || was_visited(neighbor_id)) {
                continue;
            }
            mark_visited(neighbor_id);
            frontier_.push_back(neighbor_id);
        }
    }

    std::sort(
        result.newly_opened_voxel_ids.begin(),
        result.newly_opened_voxel_ids.end());
    return result;
}

void OpeningRegionRepair::prepare_workspace(const GasGrid& gas_grid)
{
    const auto voxel_count = static_cast<std::size_t>(gas_grid.voxel_count());
    if (visit_epochs_.size() == voxel_count) {
        return;
    }
    visit_epochs_.assign(voxel_count, 0);
    repair_epoch_ = 0;
}

void OpeningRegionRepair::begin_repair_epoch()
{
    ++repair_epoch_;
    if (repair_epoch_ != 0) {
        return;
    }
    std::fill(visit_epochs_.begin(), visit_epochs_.end(), 0);
    repair_epoch_ = 1;
}

bool OpeningRegionRepair::was_visited(VoxelId voxel_id) const noexcept
{
    return visit_epochs_[static_cast<std::size_t>(voxel_id)] == repair_epoch_;
}

void OpeningRegionRepair::mark_visited(VoxelId voxel_id) noexcept
{
    visit_epochs_[static_cast<std::size_t>(voxel_id)] = repair_epoch_;
}

}  // namespace gasaccess
