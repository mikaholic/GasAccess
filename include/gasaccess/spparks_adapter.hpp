#ifndef GASACCESS_SPPARKS_ADAPTER_HPP
#define GASACCESS_SPPARKS_ADAPTER_HPP

#include "gasaccess/mpi_gas_grid.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace gasaccess {

struct SpparksAdapterConfig {
    GridSpacing requested_spacing{};
    PeriodicAxes gas_periodic{};
    ReservoirFaces reservoir_faces{};
    std::vector<VoxelCoord> explicit_source_voxels{};
    double atom_ghost_distance = 0.0;
    double maximum_excluded_radius = 0.0;
};

namespace detail {

inline std::uint64_t positive_spparks_dimension(
    int dimension,
    const char* field_name)
{
    if (dimension <= 0) {
        throw std::invalid_argument(
            std::string("SPPARKS ") + field_name + " must be positive");
    }
    return static_cast<std::uint64_t>(dimension);
}

}  // namespace detail

// DomainType follows the public field contract of SPPARKS_NS::Domain. Keeping
// this adapter templated lets GasAccess test the mapping without constructing a
// complete SPPARKS runtime, while the same function accepts the real class.
template <typename DomainType>
GridSpec make_spparks_grid_spec(
    const DomainType& domain,
    const SpparksAdapterConfig& config)
{
    if (domain.box_exist == 0) {
        throw std::invalid_argument("SPPARKS simulation box does not exist");
    }
    if (domain.dimension != 3) {
        throw std::invalid_argument("GasAccess requires a 3D SPPARKS domain");
    }

    const GridDimensions process_grid{
        detail::positive_spparks_dimension(domain.procgrid[0], "procgrid[0]"),
        detail::positive_spparks_dimension(domain.procgrid[1], "procgrid[1]"),
        detail::positive_spparks_dimension(domain.procgrid[2], "procgrid[2]")};
    const Point3 global_lower{
        domain.boxxlo,
        domain.boxylo,
        domain.boxzlo};
    const Point3 global_upper{
        domain.boxxhi,
        domain.boxyhi,
        domain.boxzhi};
    const auto geometry = make_aligned_grid_geometry(
        global_lower,
        global_upper,
        process_grid,
        config.requested_spacing);

    GridSpec grid_spec{};
    grid_spec.origin = global_lower;
    grid_spec.spacing = geometry.spacing;
    grid_spec.dimensions = geometry.dimensions;
    grid_spec.periodic = config.gas_periodic;
    grid_spec.reservoir_faces = config.reservoir_faces;
    grid_spec.explicit_source_voxels = config.explicit_source_voxels;
    return grid_spec;
}

template <typename DomainType>
MpiDecompositionSpec make_spparks_decomposition_spec(
    const DomainType& domain,
    MPI_Comm world,
    const SpparksAdapterConfig& config)
{
    MpiDecompositionSpec decomposition_spec{};
    decomposition_spec.communicator = world;
    decomposition_spec.global_lower = {
        domain.boxxlo,
        domain.boxylo,
        domain.boxzlo};
    decomposition_spec.global_upper = {
        domain.boxxhi,
        domain.boxyhi,
        domain.boxzhi};
    decomposition_spec.local_lower = {
        domain.subxlo,
        domain.subylo,
        domain.subzlo};
    decomposition_spec.local_upper = {
        domain.subxhi,
        domain.subyhi,
        domain.subzhi};
    decomposition_spec.process_grid = {
        detail::positive_spparks_dimension(domain.procgrid[0], "procgrid[0]"),
        detail::positive_spparks_dimension(domain.procgrid[1], "procgrid[1]"),
        detail::positive_spparks_dimension(domain.procgrid[2], "procgrid[2]")};
    decomposition_spec.process_location = {
        static_cast<std::int64_t>(domain.myloc[0]),
        static_cast<std::int64_t>(domain.myloc[1]),
        static_cast<std::int64_t>(domain.myloc[2])};
    decomposition_spec.atom_ghost_distance = config.atom_ghost_distance;
    decomposition_spec.maximum_excluded_radius =
        config.maximum_excluded_radius;
    return decomposition_spec;
}

class SpparksAtomBuffer {
public:
    template <typename AppType, typename RadiusAccessor>
    void assign(
        const AppType& app,
        RadiusAccessor radius_accessor,
        bool include_ghost_atoms = true)
    {
        if (app.nlocal < 0 || app.nghost < 0) {
            throw std::invalid_argument("SPPARKS atom counts must be nonnegative");
        }
        const auto owned_count = static_cast<std::size_t>(app.nlocal);
        const auto ghost_count = include_ghost_atoms
            ? static_cast<std::size_t>(app.nghost)
            : std::size_t{0};
        if (ghost_count > std::numeric_limits<std::size_t>::max() - owned_count) {
            throw std::overflow_error("SPPARKS atom count exceeds size_t");
        }
        const auto atom_count = owned_count + ghost_count;
        if (atom_count != 0 && app.xyz == nullptr) {
            throw std::invalid_argument("SPPARKS xyz array is null");
        }

        atoms_.clear();
        atoms_.reserve(atom_count);
        for (std::size_t atom_index = 0;
             atom_index < atom_count;
             ++atom_index) {
            if (app.xyz[atom_index] == nullptr) {
                throw std::invalid_argument("SPPARKS xyz row is null");
            }
            atoms_.push_back({
                {app.xyz[atom_index][0],
                 app.xyz[atom_index][1],
                 app.xyz[atom_index][2]},
                radius_accessor(atom_index)});
        }
    }

    AtomView atom_view() const noexcept
    {
        return {atoms_.data(), atoms_.size()};
    }

    const std::vector<Atom>& atoms() const noexcept
    {
        return atoms_;
    }

private:
    std::vector<Atom> atoms_{};
};

}  // namespace gasaccess

#endif
