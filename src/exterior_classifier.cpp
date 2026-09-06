#include "gasaccess/exterior_classifier.hpp"

#include <cstddef>
#include <vector>

namespace gasaccess {

ClassificationSummary ExteriorClassifier::classify(GasGrid& gas_grid) const
{
    ClassificationSummary summary{};

    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        if (gas_grid.gas_state(voxel_id) == GasState::Solid) {
            ++summary.solid_count;
        } else {
            gas_grid.set_gas_state(voxel_id, GasState::ClosedVoid);
            ++summary.closed_void_count;
        }
    }

    std::vector<VoxelId> frontier;
    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        if (gas_grid.gas_state(voxel_id) == GasState::ClosedVoid
            && gas_grid.is_reservoir_source(voxel_id)) {
            gas_grid.set_gas_state(voxel_id, GasState::OutsideAccessible);
            frontier.push_back(voxel_id);
            ++summary.outside_accessible_count;
            --summary.closed_void_count;
        }
    }

    std::size_t frontier_index = 0;
    while (frontier_index < frontier.size()) {
        const auto current_voxel_id = frontier[frontier_index];
        ++frontier_index;

        const auto neighbor_list = gas_grid.neighbors(current_voxel_id);
        for (std::size_t neighbor_index = 0;
             neighbor_index < neighbor_list.count;
             ++neighbor_index) {
            const auto neighbor_id = neighbor_list.ids[neighbor_index];
            if (gas_grid.gas_state(neighbor_id) != GasState::ClosedVoid) {
                continue;
            }

            gas_grid.set_gas_state(neighbor_id, GasState::OutsideAccessible);
            frontier.push_back(neighbor_id);
            ++summary.outside_accessible_count;
            --summary.closed_void_count;
        }
    }

    return summary;
}

}  // namespace gasaccess
