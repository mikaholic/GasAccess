#ifndef GASACCESS_VOXEL_CHANGE_HPP
#define GASACCESS_VOXEL_CHANGE_HPP

#include "gasaccess/gas_grid.hpp"

#include <cstddef>

namespace gasaccess {

// Records a gas voxel removed by solidification and its state immediately
// before the occupancy change.
struct RemovedVoxel {
    VoxelId voxel_id = 0;
    GasState previous_state = GasState::Unclassified;
};

struct RemovedVoxelView {
    const RemovedVoxel* removed_voxels = nullptr;
    std::size_t count = 0;
};

}  // namespace gasaccess

#endif
