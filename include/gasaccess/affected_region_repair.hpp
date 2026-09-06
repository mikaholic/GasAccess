#ifndef GASACCESS_AFFECTED_REGION_REPAIR_HPP
#define GASACCESS_AFFECTED_REGION_REPAIR_HPP

#include "gasaccess/local_topology_filter.hpp"

#include <cstdint>
#include <vector>

namespace gasaccess {

struct AffectedRegionRepairResult {
    // Total BFS work; a voxel may be counted again if separate early-stopped
    // searches overlap before reaching a confirmed exterior path.
    VoxelId visited_voxel_count = 0;
    std::vector<VoxelId> newly_closed_voxel_ids{};
};

class AffectedRegionRepair {
public:
    // Repairs stale OutsideAccessible labels after all removed voxels have
    // already become Solid. Traversal storage is retained between calls.
    AffectedRegionRepairResult repair(
        GasGrid& gas_grid,
        RemovedVoxelView removed_voxel_view);

private:
    void prepare_workspace(const GasGrid& gas_grid);
    void begin_repair_epoch();
    void begin_search_epoch();
    bool is_confirmed_outside(VoxelId voxel_id) const noexcept;
    void confirm_region_outside() noexcept;

    std::vector<std::uint32_t> search_epochs_{};
    std::vector<std::uint32_t> confirmed_outside_epochs_{};
    std::vector<VoxelId> seed_voxels_{};
    std::vector<VoxelId> frontier_{};
    std::uint32_t search_epoch_ = 0;
    std::uint32_t repair_epoch_ = 0;
};

}  // namespace gasaccess

#endif
