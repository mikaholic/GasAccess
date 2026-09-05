#ifndef GASACCESS_ATOM_VOXELIZER_HPP
#define GASACCESS_ATOM_VOXELIZER_HPP

#include "gasaccess/voxel_change.hpp"

#include <cstddef>
#include <vector>

namespace gasaccess {

struct Atom {
    // Position and radius use the same physical units as GridSpec.
    Point3 position{};
    double radius = 0.0;
};

struct AtomView {
    // Non-owning view; the atom storage must remain valid for voxelize().
    const Atom* atoms = nullptr;
    std::size_t count = 0;
};

// A single atomic-lattice update. Added and removed atoms are combined before
// any voxel occupancy is mutated, so overlap and cancellation are handled by
// their net blocker-count change.
struct AtomChangeBatch {
    AtomView added_atoms{};
    AtomView removed_atoms{};
};

struct VoxelOccupancyChange {
    VoxelId voxel_id = 0;
    VoxelBlockerCount previous_blocker_count = 0;
    VoxelBlockerCount blocker_count = 0;
    GasState previous_state = GasState::Unclassified;
};

struct AtomChangeOccupancyResult {
    VoxelId blocker_count_changed_voxel_count = 0;
    VoxelId newly_solid_count = 0;
    VoxelId newly_gas_count = 0;

    bool geometry_changed() const noexcept
    {
        return newly_solid_count != 0 || newly_gas_count != 0;
    }
};

class AtomVoxelizer {
public:
    explicit AtomVoxelizer(double precursor_radius);

    double precursor_radius() const noexcept;

    // Voxelization is additive: existing solid voxels remain solid. All inputs
    // are validated before mutation. The return value is the number of voxels
    // newly changed to GasState::Solid.
    VoxelId voxelize(GasGrid& gas_grid, AtomView atom_view) const;

    // Additionally records every newly solid voxel and its previous state.
    // The caller-owned output is cleared after input validation and may retain
    // capacity for reuse across deposition events.
    VoxelId voxelize(
        GasGrid& gas_grid,
        AtomView atom_view,
        std::vector<RemovedVoxel>& removed_voxels) const;

    // Applies additions and removals transactionally. Every atom in the batch
    // must correspond exactly once to the caller's atom-list update. Invalid
    // input, blocker-count underflow, or overflow leaves the grid unchanged.
    AtomChangeOccupancyResult apply_atom_changes(
        GasGrid& gas_grid,
        const AtomChangeBatch& atom_changes) const;
    AtomChangeOccupancyResult apply_atom_changes(
        GasGrid& gas_grid,
        const AtomChangeBatch& atom_changes,
        std::vector<VoxelOccupancyChange>& voxel_changes) const;

private:
    VoxelId voxelize_impl(
        GasGrid& gas_grid,
        AtomView atom_view,
        std::vector<RemovedVoxel>* removed_voxels) const;

    double precursor_radius_ = 0.0;
};

}  // namespace gasaccess

#endif
