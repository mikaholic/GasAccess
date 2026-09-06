#ifndef GASACCESS_DISTRIBUTED_CLOSING_REGION_REPAIR_HPP
#define GASACCESS_DISTRIBUTED_CLOSING_REGION_REPAIR_HPP

#include "gasaccess/mpi_gas_grid.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace gasaccess {

struct NewlySolidVoxelCoordView {
    const DistributedRemovedVoxel* voxels = nullptr;
    std::size_t count = 0;
};

struct DistributedClosingRegionRepairResult {
    std::uint64_t local_visited_voxel_count = 0;
    std::uint64_t local_closed_voxel_count = 0;
    std::uint64_t participating_rank_count = 0;
    std::uint64_t seed_search_count = 0;
    std::uint64_t communication_round_count = 0;
    std::uint64_t sent_frontier_entry_count = 0;
    std::uint64_t received_frontier_entry_count = 0;
    std::vector<VoxelCoord> newly_closed_owned_voxel_coords{};
};

class DistributedClosingRegionRepair {
public:
    // Collective over the grid communicator. The supplied voxels must already
    // be owned Solid voxels and retain their pre-change gas states. Newly gas
    // voxels from the same atomic batch may still be Unclassified; the closing
    // pass deliberately excludes them and lets a following opening pass
    // restore connectivity in the final geometry.
    DistributedClosingRegionRepairResult repair(
        DistributedGasGrid& gas_grid,
        NewlySolidVoxelCoordView newly_solid_voxel_view);

private:
    void gather_newly_solid_voxels(
        const DistributedGasGrid& gas_grid,
        NewlySolidVoxelCoordView newly_solid_voxel_view);
    void prepare_workspace(const DistributedGasGrid& gas_grid);
    void begin_repair_epoch();
    void begin_search_epoch();

    std::vector<DistributedRemovedVoxel> global_newly_solid_voxels_{};
    std::vector<std::uint32_t> search_epochs_{};
    std::vector<std::uint32_t> confirmed_outside_epochs_{};
    std::vector<VoxelCoord> seed_voxels_{};
    std::vector<VoxelCoord> frontier_{};
    std::vector<VoxelCoord> visited_voxels_{};
    std::array<std::vector<std::uint64_t>, 6> send_buffers_{};
    std::array<std::vector<std::uint64_t>, 6> receive_buffers_{};
    std::uint32_t search_epoch_ = 0;
    std::uint32_t repair_epoch_ = 0;
};

}  // namespace gasaccess

#endif
