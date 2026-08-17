#ifndef GASACCESS_LOCAL_TOPOLOGY_FILTER_HPP
#define GASACCESS_LOCAL_TOPOLOGY_FILTER_HPP

#include "gasaccess/voxel_change.hpp"

#include <cstdint>

namespace gasaccess {

enum class TopologyDecision : std::uint8_t {
    Safe = 0,
    RequiresConnectivityRepair = 1
};

struct TopologyCheckResult {
    TopologyDecision decision = TopologyDecision::RequiresConnectivityRepair;
    std::size_t accessible_neighbor_count = 0;
    std::size_t visited_voxel_count = 0;

    bool is_safe() const noexcept;
};

class LocalTopologyFilter {
public:
    // The grid is the post-deposition grid, while previous_state describes
    // each newly solid voxel immediately before deposition. The initial
    // implementation proves only single-voxel removals; larger changes return
    // RequiresConnectivityRepair conservatively.
    TopologyCheckResult evaluate(
        const GasGrid& gas_grid,
        RemovedVoxelView removed_voxel_view) const;

private:
    static bool is_in_local_neighborhood(
        const GasGrid& gas_grid,
        VoxelId center_voxel_id,
        VoxelId candidate_voxel_id);
};

}  // namespace gasaccess

#endif
