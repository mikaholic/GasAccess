#ifndef GASACCESS_ADSORPTION_UPDATER_HPP
#define GASACCESS_ADSORPTION_UPDATER_HPP

#include "gasaccess/affected_region_repair.hpp"
#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/exterior_classifier.hpp"
#include "gasaccess/local_topology_filter.hpp"

#include <cstdint>
#include <vector>

namespace gasaccess {

enum class ConnectivityRepairMode : std::uint8_t {
    AffectedRegion = 0,
    FullReclassification = 1
};

struct AdsorptionUpdateResult {
    VoxelId newly_solid_count = 0;
    std::vector<VoxelId> changed_voxel_ids{};
    ClassificationSummary classification{};
    bool full_reclassification_performed = false;
    bool affected_region_repair_performed = false;
    VoxelId repair_visited_voxel_count = 0;
    VoxelId repair_closed_voxel_count = 0;

    bool geometry_changed() const noexcept;
    bool used_full_reclassification() const noexcept;
    bool used_affected_region_repair() const noexcept;
};

class AdsorptionUpdater {
public:
    explicit AdsorptionUpdater(
        double precursor_radius,
        ConnectivityRepairMode repair_mode = ConnectivityRepairMode::AffectedRegion);

    double precursor_radius() const noexcept;
    ConnectivityRepairMode repair_mode() const noexcept;

    // The grid must already be fully classified. Adsorbed atoms are supplied
    // by the caller; this method does not model adsorption physics. Changed
    // voxel identifiers are unique and sorted in ascending order.
    AdsorptionUpdateResult apply_adsorption(
        GasGrid& gas_grid,
        AtomView adsorbed_atoms) const;

private:
    AtomVoxelizer atom_voxelizer_;
    ConnectivityRepairMode repair_mode_ = ConnectivityRepairMode::AffectedRegion;
    LocalTopologyFilter local_topology_filter_;
    // Scratch storage is reused by logically const update calls. A single
    // updater instance must not be called concurrently.
    mutable AffectedRegionRepair affected_region_repair_;
    mutable std::vector<RemovedVoxel> removed_voxels_{};
};

}  // namespace gasaccess

#endif
