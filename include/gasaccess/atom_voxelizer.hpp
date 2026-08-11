#ifndef GASACCESS_ATOM_VOXELIZER_HPP
#define GASACCESS_ATOM_VOXELIZER_HPP

#include "gasaccess/gas_grid.hpp"

#include <cstddef>

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

class AtomVoxelizer {
public:
    explicit AtomVoxelizer(double precursor_radius);

    double precursor_radius() const noexcept;

    // Voxelization is additive: existing solid voxels remain solid. All inputs
    // are validated before mutation. The return value is the number of voxels
    // newly changed to GasState::Solid.
    VoxelId voxelize(GasGrid& gas_grid, AtomView atom_view) const;

private:
    double precursor_radius_ = 0.0;
};

}  // namespace gasaccess

#endif
