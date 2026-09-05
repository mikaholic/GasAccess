#ifndef GASACCESS_DISTRIBUTED_ATOM_CHANGE_UPDATER_HPP
#define GASACCESS_DISTRIBUTED_ATOM_CHANGE_UPDATER_HPP

#include "gasaccess/atom_change_updater.hpp"
#include "gasaccess/distributed_closing_region_repair.hpp"
#include "gasaccess/distributed_exterior_classifier.hpp"
#include "gasaccess/distributed_opening_region_repair.hpp"

#include <cstdint>
#include <vector>

namespace gasaccess {

struct DistributedAtomChangeUpdateResult {
    std::uint64_t local_blocker_count_changed_voxel_count = 0;
    std::uint64_t global_blocker_count_changed_voxel_count = 0;
    std::uint64_t local_newly_solid_count = 0;
    std::uint64_t global_newly_solid_count = 0;
    std::uint64_t local_newly_gas_count = 0;
    std::uint64_t global_newly_gas_count = 0;
    std::vector<VoxelCoord> changed_owned_voxel_coords{};
    DistributedClassificationSummary classification{};
    AccessibilityRepairKind repair_kind = AccessibilityRepairKind::None;
    bool full_reclassification_performed = false;
    bool distributed_closing_repair_performed = false;
    bool distributed_opening_repair_performed = false;
    std::uint64_t local_closing_visited_voxel_count = 0;
    std::uint64_t local_opening_visited_voxel_count = 0;
    std::uint64_t local_repair_closed_voxel_count = 0;
    std::uint64_t local_repair_opened_voxel_count = 0;
    std::uint64_t closing_participating_rank_count = 0;
    std::uint64_t opening_participating_rank_count = 0;
    std::uint64_t closing_seed_search_count = 0;
    std::uint64_t closing_communication_round_count = 0;
    std::uint64_t opening_communication_round_count = 0;
    std::uint64_t closing_sent_frontier_entry_count = 0;
    std::uint64_t closing_received_frontier_entry_count = 0;
    std::uint64_t opening_sent_frontier_entry_count = 0;
    std::uint64_t opening_received_frontier_entry_count = 0;

    bool geometry_changed() const noexcept;
    bool used_full_reclassification() const noexcept;
    bool used_distributed_closing_repair() const noexcept;
    bool used_distributed_opening_repair() const noexcept;
};

class DistributedAtomChangeUpdater {
public:
    explicit DistributedAtomChangeUpdater(
        double precursor_radius,
        AtomChangeRepairMode repair_mode =
            AtomChangeRepairMode::Incremental);

    double precursor_radius() const noexcept;
    AtomChangeRepairMode repair_mode() const noexcept;

    // Collective over the grid communicator. Added atoms use final records;
    // removed atoms use old positions and radii. All ranks call in the same
    // order, including ranks with empty local views. Occupancy is committed
    // from the net batch before either connectivity pass observes it.
    DistributedAtomChangeUpdateResult apply_atom_changes(
        DistributedGasGrid& gas_grid,
        const AtomChangeBatch& atom_changes);

private:
    double precursor_radius_ = 0.0;
    AtomChangeRepairMode repair_mode_ = AtomChangeRepairMode::Incremental;
    DistributedExteriorClassifier exterior_classifier_{};
    DistributedClosingRegionRepair closing_repair_{};
    DistributedOpeningRegionRepair opening_repair_{};
    std::vector<DistributedVoxelOccupancyChange> occupancy_changes_{};
    std::vector<DistributedRemovedVoxel> newly_solid_voxels_{};
    std::vector<VoxelCoord> newly_gas_voxel_coords_{};
    std::vector<VoxelCoord> closing_changed_voxel_coords_{};
    std::vector<GasState> previous_states_{};
};

}  // namespace gasaccess

#endif
