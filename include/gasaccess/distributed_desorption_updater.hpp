#ifndef GASACCESS_DISTRIBUTED_DESORPTION_UPDATER_HPP
#define GASACCESS_DISTRIBUTED_DESORPTION_UPDATER_HPP

#include "gasaccess/desorption_updater.hpp"
#include "gasaccess/distributed_exterior_classifier.hpp"
#include "gasaccess/distributed_opening_region_repair.hpp"

#include <cstdint>
#include <vector>

namespace gasaccess {

struct DistributedDesorptionUpdateResult {
    std::uint64_t local_blocker_count_changed_voxel_count = 0;
    std::uint64_t global_blocker_count_changed_voxel_count = 0;
    std::uint64_t local_newly_gas_count = 0;
    std::uint64_t global_newly_gas_count = 0;
    std::vector<VoxelCoord> changed_owned_voxel_coords{};
    DistributedClassificationSummary classification{};
    bool full_reclassification_performed = false;
    bool distributed_opening_repair_performed = false;
    std::uint64_t local_opening_visited_voxel_count = 0;
    std::uint64_t local_repair_opened_voxel_count = 0;
    std::uint64_t opening_participating_rank_count = 0;
    std::uint64_t repair_communication_round_count = 0;
    std::uint64_t sent_frontier_entry_count = 0;
    std::uint64_t received_frontier_entry_count = 0;

    bool geometry_changed() const noexcept;
    bool used_full_reclassification() const noexcept;
    bool used_distributed_opening_repair() const noexcept;
};

class DistributedDesorptionUpdater {
public:
    explicit DistributedDesorptionUpdater(
        double precursor_radius,
        DesorptionRepairMode repair_mode =
            DesorptionRepairMode::OpeningRegion);

    double precursor_radius() const noexcept;
    DesorptionRepairMode repair_mode() const noexcept;

    // Collective over the grid communicator. Removed atoms must be supplied
    // at their old positions and synchronized far enough to cover every
    // affected owned voxel. Every rank must call in the same order, including
    // ranks with an empty local atom view. The grid must be fully classified.
    DistributedDesorptionUpdateResult apply_desorption(
        DistributedGasGrid& gas_grid,
        AtomView removed_atoms);

private:
    double precursor_radius_ = 0.0;
    DesorptionRepairMode repair_mode_ = DesorptionRepairMode::OpeningRegion;
    DistributedExteriorClassifier exterior_classifier_{};
    DistributedOpeningRegionRepair opening_region_repair_{};
    std::vector<DistributedVoxelOccupancyChange> occupancy_changes_{};
    std::vector<VoxelCoord> newly_gas_voxel_coords_{};
    std::vector<GasState> previous_states_{};
};

}  // namespace gasaccess

#endif
