#ifndef GASACCESS_ATOM_CHANGE_UPDATER_HPP
#define GASACCESS_ATOM_CHANGE_UPDATER_HPP

#include "gasaccess/affected_region_repair.hpp"
#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/exterior_classifier.hpp"
#include "gasaccess/opening_region_repair.hpp"

#include <cstdint>
#include <vector>

namespace gasaccess {

enum class AccessibilityRepairKind : std::uint8_t {
    None = 0,
    Closing = 1,
    Opening = 2,
    Mixed = 3,
    FullReclassification = 4
};

enum class AtomChangeRepairMode : std::uint8_t {
    Incremental = 0,
    FullReclassification = 1
};

struct AtomChangeUpdateResult {
    VoxelId blocker_count_changed_voxel_count = 0;
    VoxelId newly_solid_count = 0;
    VoxelId newly_gas_count = 0;
    std::vector<VoxelId> changed_voxel_ids{};
    ClassificationSummary classification{};
    AccessibilityRepairKind repair_kind = AccessibilityRepairKind::None;
    bool full_reclassification_performed = false;
    bool closing_repair_performed = false;
    bool opening_repair_performed = false;
    VoxelId closing_visited_voxel_count = 0;
    VoxelId opening_visited_voxel_count = 0;
    VoxelId repair_closed_voxel_count = 0;
    VoxelId repair_opened_voxel_count = 0;

    bool geometry_changed() const noexcept;
    bool used_full_reclassification() const noexcept;
    bool used_closing_repair() const noexcept;
    bool used_opening_repair() const noexcept;
};

class AtomChangeUpdater {
public:
    explicit AtomChangeUpdater(
        double precursor_radius,
        AtomChangeRepairMode repair_mode =
            AtomChangeRepairMode::Incremental);

    double precursor_radius() const noexcept;
    AtomChangeRepairMode repair_mode() const noexcept;

    // Applies additions and removals as one atomic occupancy transaction. The
    // grid must be fully classified. Removed atoms use their old positions and
    // radii; an atom move is remove-old plus add-new in the same batch.
    AtomChangeUpdateResult apply_atom_changes(
        GasGrid& gas_grid,
        const AtomChangeBatch& atom_changes) const;

    // Convenience wrappers for callers that use one persistent updater for
    // every event kind.
    AtomChangeUpdateResult apply_deposition(
        GasGrid& gas_grid,
        AtomView added_atoms) const;
    AtomChangeUpdateResult apply_desorption(
        GasGrid& gas_grid,
        AtomView removed_atoms) const;

private:
    AtomVoxelizer atom_voxelizer_;
    AtomChangeRepairMode repair_mode_ = AtomChangeRepairMode::Incremental;
    mutable AffectedRegionRepair closing_repair_{};
    mutable OpeningRegionRepair opening_repair_{};
    mutable std::vector<VoxelOccupancyChange> occupancy_changes_{};
    mutable std::vector<RemovedVoxel> newly_solid_voxels_{};
    mutable std::vector<VoxelId> newly_gas_voxel_ids_{};
    mutable std::vector<VoxelId> closing_changed_voxel_ids_{};
    mutable std::vector<GasState> previous_states_{};
};

}  // namespace gasaccess

#endif
