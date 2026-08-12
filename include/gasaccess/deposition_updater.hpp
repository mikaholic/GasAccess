#ifndef GASACCESS_DEPOSITION_UPDATER_HPP
#define GASACCESS_DEPOSITION_UPDATER_HPP

#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/exterior_classifier.hpp"

#include <vector>

namespace gasaccess {

struct DepositionUpdateResult {
    VoxelId newly_solid_count = 0;
    std::vector<VoxelId> changed_voxel_ids{};
    ClassificationSummary classification{};

    bool geometry_changed() const noexcept;
};

class DepositionUpdater {
public:
    explicit DepositionUpdater(double precursor_radius);

    double precursor_radius() const noexcept;

    // The grid must already be fully classified. Deposited atoms are supplied
    // by the caller; this method does not model deposition physics. Changed
    // voxel identifiers are unique and sorted in ascending order.
    DepositionUpdateResult apply_deposition(
        GasGrid& gas_grid,
        AtomView deposited_atoms) const;

private:
    AtomVoxelizer atom_voxelizer_;
};

}  // namespace gasaccess

#endif
