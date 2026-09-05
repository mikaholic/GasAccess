#ifndef GASACCESS_DISTRIBUTED_OPENING_REGION_REPAIR_HPP
#define GASACCESS_DISTRIBUTED_OPENING_REGION_REPAIR_HPP

#include "gasaccess/mpi_gas_grid.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace gasaccess {

struct NewlyGasVoxelCoordView {
    const VoxelCoord* voxel_coords = nullptr;
    std::size_t count = 0;
};

struct DistributedOpeningRegionRepairResult {
    std::uint64_t local_visited_voxel_count = 0;
    std::uint64_t local_opened_voxel_count = 0;
    std::uint64_t participating_rank_count = 0;
    std::uint64_t communication_round_count = 0;
    std::uint64_t sent_frontier_entry_count = 0;
    std::uint64_t received_frontier_entry_count = 0;
    std::vector<VoxelCoord> newly_opened_owned_voxel_coords{};
};

class DistributedOpeningRegionRepair {
public:
    // Collective over the grid communicator. Newly gas voxels must be owned,
    // unblocked ClosedVoid voxels. The repair promotes only components that
    // touch a reservoir source or the existing OutsideAccessible region.
    // Face ghosts are deliberately not synchronized here so an updater can
    // perform one final exchange after all state changes are complete.
    DistributedOpeningRegionRepairResult repair(
        DistributedGasGrid& gas_grid,
        NewlyGasVoxelCoordView newly_gas_voxel_view);

private:
    void prepare_workspace(const DistributedGasGrid& gas_grid);
    void begin_repair_epoch();

    std::vector<std::uint32_t> visit_epochs_{};
    std::vector<VoxelCoord> seed_voxels_{};
    std::vector<VoxelCoord> frontier_{};
    std::array<std::vector<std::uint64_t>, 6> send_buffers_{};
    std::array<std::vector<std::uint64_t>, 6> receive_buffers_{};
    std::uint32_t repair_epoch_ = 0;
};

}  // namespace gasaccess

#endif
