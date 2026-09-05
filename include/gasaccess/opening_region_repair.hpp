#ifndef GASACCESS_OPENING_REGION_REPAIR_HPP
#define GASACCESS_OPENING_REGION_REPAIR_HPP

#include "gasaccess/gas_grid.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace gasaccess {

struct NewlyGasVoxelView {
    const VoxelId* voxel_ids = nullptr;
    std::size_t count = 0;
};

struct OpeningRegionRepairResult {
    VoxelId visited_voxel_count = 0;
    std::vector<VoxelId> newly_opened_voxel_ids{};
};

class OpeningRegionRepair {
public:
    // Newly gas voxels must already have blocker count zero and be initialized
    // as ClosedVoid. Only components connected to a reservoir source or an
    // existing OutsideAccessible voxel are opened. Traversal storage is
    // retained between calls.
    OpeningRegionRepairResult repair(
        GasGrid& gas_grid,
        NewlyGasVoxelView newly_gas_voxel_view);

private:
    void prepare_workspace(const GasGrid& gas_grid);
    void begin_repair_epoch();
    bool was_visited(VoxelId voxel_id) const noexcept;
    void mark_visited(VoxelId voxel_id) noexcept;

    std::vector<std::uint32_t> visit_epochs_{};
    std::vector<VoxelId> seed_voxels_{};
    std::vector<VoxelId> frontier_{};
    std::uint32_t repair_epoch_ = 0;
};

}  // namespace gasaccess

#endif
