#ifndef GASACCESS_DISTRIBUTED_EXTERIOR_CLASSIFIER_HPP
#define GASACCESS_DISTRIBUTED_EXTERIOR_CLASSIFIER_HPP

#include "gasaccess/mpi_gas_grid.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace gasaccess {

struct DistributedClassificationSummary {
    std::uint64_t local_solid_count = 0;
    std::uint64_t local_outside_accessible_count = 0;
    std::uint64_t local_closed_void_count = 0;
    std::uint64_t local_visited_voxel_count = 0;
    std::uint64_t communication_round_count = 0;
    std::uint64_t sent_frontier_entry_count = 0;
    std::uint64_t received_frontier_entry_count = 0;
};

class DistributedExteriorClassifier {
public:
    // Performs a full distributed classification. GasState::Solid is treated
    // as occupancy; every other owned cached state is recomputed. The final
    // face-ghost states are synchronized before this function returns.
    DistributedClassificationSummary classify(DistributedGasGrid& gas_grid);

private:
    std::vector<VoxelCoord> frontier_{};
    std::size_t frontier_head_ = 0;
    std::array<std::vector<std::uint64_t>, 6> send_buffers_{};
    std::array<std::vector<std::uint64_t>, 6> receive_buffers_{};
};

}  // namespace gasaccess

#endif
