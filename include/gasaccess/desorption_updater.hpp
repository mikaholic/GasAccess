#ifndef GASACCESS_DESORPTION_UPDATER_HPP
#define GASACCESS_DESORPTION_UPDATER_HPP

#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/exterior_classifier.hpp"
#include "gasaccess/opening_region_repair.hpp"

#include <cstdint>
#include <vector>

namespace gasaccess {

enum class DesorptionRepairMode : std::uint8_t {
    OpeningRegion = 0,
    FullReclassification = 1
};

struct DesorptionUpdateResult {
    VoxelId blocker_count_changed_voxel_count = 0;
    VoxelId newly_gas_count = 0;
    std::vector<VoxelId> changed_voxel_ids{};
    ClassificationSummary classification{};
    bool full_reclassification_performed = false;
    bool opening_region_repair_performed = false;
    VoxelId opening_visited_voxel_count = 0;
    VoxelId repair_opened_voxel_count = 0;

    bool geometry_changed() const noexcept;
    bool used_full_reclassification() const noexcept;
    bool used_opening_region_repair() const noexcept;
};

class DesorptionUpdater {
public:
    explicit DesorptionUpdater(
        double precursor_radius,
        DesorptionRepairMode repair_mode =
            DesorptionRepairMode::OpeningRegion);

    double precursor_radius() const noexcept;
    DesorptionRepairMode repair_mode() const noexcept;

    // The grid must already be fully classified and removed atoms must be
    // supplied at their old positions and radii. The caller remains the atom
    // source of truth and must report every removed atom exactly once.
    DesorptionUpdateResult apply_desorption(
        GasGrid& gas_grid,
        AtomView removed_atoms) const;

private:
    AtomVoxelizer atom_voxelizer_;
    DesorptionRepairMode repair_mode_ = DesorptionRepairMode::OpeningRegion;
    // Scratch storage is reused by logically const update calls. A single
    // updater instance must not be called concurrently.
    mutable OpeningRegionRepair opening_region_repair_;
    mutable std::vector<VoxelOccupancyChange> occupancy_changes_{};
    mutable std::vector<VoxelId> newly_gas_voxel_ids_{};
};

}  // namespace gasaccess

#endif
