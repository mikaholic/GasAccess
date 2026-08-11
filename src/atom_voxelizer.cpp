#include "gasaccess/atom_voxelizer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace gasaccess {
namespace {

struct AxisCandidates {
    std::uint64_t start = 0;
    std::uint64_t count = 0;
    std::uint64_t dimension = 0;
};

AxisCandidates make_axis_candidates(
    std::uint64_t center_index,
    std::uint64_t dimension,
    double excluded_radius,
    double spacing,
    bool periodic)
{
    const double reach_value = std::ceil(excluded_radius / spacing) + 1.0;
    if (!std::isfinite(reach_value)
        || reach_value >= static_cast<double>(dimension)) {
        return {0, dimension, dimension};
    }

    const auto reach = static_cast<std::uint64_t>(reach_value);
    if (periodic) {
        if (reach >= dimension / 2) {
            return {0, dimension, dimension};
        }

        const auto start = center_index >= reach
            ? center_index - reach
            : dimension - (reach - center_index);
        return {start, 2 * reach + 1, dimension};
    }

    const auto start = center_index >= reach ? center_index - reach : 0;
    const auto distance_to_high = dimension - 1 - center_index;
    const auto end = reach >= distance_to_high
        ? dimension - 1
        : center_index + reach;
    return {start, end - start + 1, dimension};
}

std::uint64_t axis_index(const AxisCandidates& candidates, std::uint64_t offset)
{
    const auto unwrapped_index = candidates.start + offset;
    return unwrapped_index >= candidates.dimension
        ? unwrapped_index - candidates.dimension
        : unwrapped_index;
}

double wrap_position(double position, double origin, double length)
{
    double relative = std::fmod(position - origin, length);
    if (relative < 0.0) {
        relative += length;
    }
    return origin + relative;
}

double axis_distance(double lhs, double rhs, double length, bool periodic)
{
    const double direct_distance = std::abs(lhs - rhs);
    if (!periodic) {
        return direct_distance;
    }
    return std::min(direct_distance, std::abs(length - direct_distance));
}

Point3 normalized_atom_position(const GasGrid& gas_grid, const Atom& atom)
{
    if (!std::isfinite(atom.position.x)
        || !std::isfinite(atom.position.y)
        || !std::isfinite(atom.position.z)) {
        throw std::invalid_argument("atom position must be finite");
    }
    if (!gas_grid.locate_voxel(atom.position)) {
        throw std::invalid_argument("atom position is outside a non-periodic grid boundary");
    }

    const auto& grid_spec = gas_grid.grid_spec();
    const Point3 lengths{
        grid_spec.spacing * static_cast<double>(grid_spec.dimensions.x),
        grid_spec.spacing * static_cast<double>(grid_spec.dimensions.y),
        grid_spec.spacing * static_cast<double>(grid_spec.dimensions.z)
    };

    Point3 position = atom.position;
    if (grid_spec.periodic.x) {
        position.x = wrap_position(position.x, grid_spec.origin.x, lengths.x);
    }
    if (grid_spec.periodic.y) {
        position.y = wrap_position(position.y, grid_spec.origin.y, lengths.y);
    }
    if (grid_spec.periodic.z) {
        position.z = wrap_position(position.z, grid_spec.origin.z, lengths.z);
    }
    return position;
}

double validate_excluded_radius(const Atom& atom, double precursor_radius)
{
    if (!std::isfinite(atom.radius) || atom.radius < 0.0) {
        throw std::invalid_argument("atom radius must be finite and nonnegative");
    }

    const double excluded_radius = atom.radius + precursor_radius;
    if (!std::isfinite(excluded_radius)) {
        throw std::invalid_argument("atom and precursor radii have a non-finite sum");
    }
    if (!std::isfinite(excluded_radius * excluded_radius)) {
        throw std::invalid_argument("excluded radius is too large for stable distance tests");
    }
    return excluded_radius;
}

}  // namespace

AtomVoxelizer::AtomVoxelizer(double precursor_radius)
    : precursor_radius_(precursor_radius)
{
    if (!std::isfinite(precursor_radius_) || precursor_radius_ < 0.0) {
        throw std::invalid_argument("precursor radius must be finite and nonnegative");
    }
}

double AtomVoxelizer::precursor_radius() const noexcept
{
    return precursor_radius_;
}

VoxelId AtomVoxelizer::voxelize(GasGrid& gas_grid, AtomView atom_view) const
{
    if (atom_view.count != 0 && atom_view.atoms == nullptr) {
        throw std::invalid_argument("atom view has a null pointer with nonzero count");
    }

    for (std::size_t atom_index = 0; atom_index < atom_view.count; ++atom_index) {
        const auto& atom = atom_view.atoms[atom_index];
        static_cast<void>(validate_excluded_radius(atom, precursor_radius_));
        static_cast<void>(normalized_atom_position(gas_grid, atom));
    }

    const auto& grid_spec = gas_grid.grid_spec();
    const Point3 lengths{
        grid_spec.spacing * static_cast<double>(grid_spec.dimensions.x),
        grid_spec.spacing * static_cast<double>(grid_spec.dimensions.y),
        grid_spec.spacing * static_cast<double>(grid_spec.dimensions.z)
    };

    VoxelId newly_solid_count = 0;
    for (std::size_t atom_index = 0; atom_index < atom_view.count; ++atom_index) {
        const auto& atom = atom_view.atoms[atom_index];
        const double excluded_radius = validate_excluded_radius(atom, precursor_radius_);
        const double excluded_radius_squared = excluded_radius * excluded_radius;
        const Point3 atom_position = normalized_atom_position(gas_grid, atom);
        const auto atom_voxel = gas_grid.locate_voxel(atom_position);
        if (!atom_voxel) {
            throw std::logic_error("normalized atom position is outside the grid");
        }

        const auto x_candidates = make_axis_candidates(
            static_cast<std::uint64_t>(atom_voxel->x),
            grid_spec.dimensions.x,
            excluded_radius,
            grid_spec.spacing,
            grid_spec.periodic.x);
        const auto y_candidates = make_axis_candidates(
            static_cast<std::uint64_t>(atom_voxel->y),
            grid_spec.dimensions.y,
            excluded_radius,
            grid_spec.spacing,
            grid_spec.periodic.y);
        const auto z_candidates = make_axis_candidates(
            static_cast<std::uint64_t>(atom_voxel->z),
            grid_spec.dimensions.z,
            excluded_radius,
            grid_spec.spacing,
            grid_spec.periodic.z);

        for (std::uint64_t z_offset = 0; z_offset < z_candidates.count; ++z_offset) {
            const auto z = axis_index(z_candidates, z_offset);
            const double z_center = grid_spec.origin.z
                + (static_cast<double>(z) + 0.5) * grid_spec.spacing;
            const double z_distance = axis_distance(
                z_center,
                atom_position.z,
                lengths.z,
                grid_spec.periodic.z);
            if (z_distance > excluded_radius) {
                continue;
            }
            const double z_distance_squared = z_distance * z_distance;

            for (std::uint64_t y_offset = 0; y_offset < y_candidates.count; ++y_offset) {
                const auto y = axis_index(y_candidates, y_offset);
                const double y_center = grid_spec.origin.y
                    + (static_cast<double>(y) + 0.5) * grid_spec.spacing;
                const double y_distance = axis_distance(
                    y_center,
                    atom_position.y,
                    lengths.y,
                    grid_spec.periodic.y);
                if (y_distance > excluded_radius) {
                    continue;
                }
                const double yz_distance_squared = y_distance * y_distance
                    + z_distance_squared;
                if (yz_distance_squared > excluded_radius_squared) {
                    continue;
                }

                for (std::uint64_t x_offset = 0;
                     x_offset < x_candidates.count;
                     ++x_offset) {
                    const auto x = axis_index(x_candidates, x_offset);
                    const double x_center = grid_spec.origin.x
                        + (static_cast<double>(x) + 0.5) * grid_spec.spacing;
                    const double x_distance = axis_distance(
                        x_center,
                        atom_position.x,
                        lengths.x,
                        grid_spec.periodic.x);
                    if (x_distance > excluded_radius) {
                        continue;
                    }

                    const double distance_squared = x_distance * x_distance
                        + yz_distance_squared;
                    if (distance_squared > excluded_radius_squared) {
                        continue;
                    }

                    const auto id = gas_grid.voxel_id({
                        static_cast<std::int64_t>(x),
                        static_cast<std::int64_t>(y),
                        static_cast<std::int64_t>(z)
                    });
                    if (gas_grid.gas_state(id) != GasState::Solid) {
                        gas_grid.set_gas_state(id, GasState::Solid);
                        ++newly_solid_count;
                    }
                }
            }
        }
    }
    return newly_solid_count;
}

}  // namespace gasaccess
