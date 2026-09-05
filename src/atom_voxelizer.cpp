#include "gasaccess/atom_voxelizer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

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
        grid_spec.spacing.x * static_cast<double>(grid_spec.dimensions.x),
        grid_spec.spacing.y * static_cast<double>(grid_spec.dimensions.y),
        grid_spec.spacing.z * static_cast<double>(grid_spec.dimensions.z)
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

void validate_atom_view(
    const GasGrid& gas_grid,
    AtomView atom_view,
    double precursor_radius)
{
    if (atom_view.count != 0 && atom_view.atoms == nullptr) {
        throw std::invalid_argument("atom view has a null pointer with nonzero count");
    }
    for (std::size_t atom_index = 0; atom_index < atom_view.count; ++atom_index) {
        const auto& atom = atom_view.atoms[atom_index];
        static_cast<void>(validate_excluded_radius(atom, precursor_radius));
        static_cast<void>(normalized_atom_position(gas_grid, atom));
    }
}

template <typename Visitor>
void for_each_covered_voxel(
    const GasGrid& gas_grid,
    AtomView atom_view,
    double precursor_radius,
    Visitor&& visitor)
{
    const auto& grid_spec = gas_grid.grid_spec();
    const Point3 lengths{
        grid_spec.spacing.x * static_cast<double>(grid_spec.dimensions.x),
        grid_spec.spacing.y * static_cast<double>(grid_spec.dimensions.y),
        grid_spec.spacing.z * static_cast<double>(grid_spec.dimensions.z)
    };

    for (std::size_t atom_index = 0; atom_index < atom_view.count; ++atom_index) {
        const auto& atom = atom_view.atoms[atom_index];
        const double excluded_radius = validate_excluded_radius(atom, precursor_radius);
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
            grid_spec.spacing.x,
            grid_spec.periodic.x);
        const auto y_candidates = make_axis_candidates(
            static_cast<std::uint64_t>(atom_voxel->y),
            grid_spec.dimensions.y,
            excluded_radius,
            grid_spec.spacing.y,
            grid_spec.periodic.y);
        const auto z_candidates = make_axis_candidates(
            static_cast<std::uint64_t>(atom_voxel->z),
            grid_spec.dimensions.z,
            excluded_radius,
            grid_spec.spacing.z,
            grid_spec.periodic.z);

        for (std::uint64_t z_offset = 0; z_offset < z_candidates.count; ++z_offset) {
            const auto z = axis_index(z_candidates, z_offset);
            const double z_center = grid_spec.origin.z
                + (static_cast<double>(z) + 0.5) * grid_spec.spacing.z;
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
                    + (static_cast<double>(y) + 0.5) * grid_spec.spacing.y;
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
                        + (static_cast<double>(x) + 0.5) * grid_spec.spacing.x;
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
                    visitor(gas_grid.voxel_id({
                        static_cast<std::int64_t>(x),
                        static_cast<std::int64_t>(y),
                        static_cast<std::int64_t>(z)
                    }));
                }
            }
        }
    }
}

struct VoxelDelta {
    VoxelId voxel_id = 0;
    int delta = 0;
};

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
    return voxelize_impl(gas_grid, atom_view, nullptr);
}

VoxelId AtomVoxelizer::voxelize(
    GasGrid& gas_grid,
    AtomView atom_view,
    std::vector<RemovedVoxel>& removed_voxels) const
{
    return voxelize_impl(gas_grid, atom_view, &removed_voxels);
}

VoxelId AtomVoxelizer::voxelize_impl(
    GasGrid& gas_grid,
    AtomView atom_view,
    std::vector<RemovedVoxel>* removed_voxels) const
{
    validate_atom_view(gas_grid, atom_view, precursor_radius_);

    if (removed_voxels != nullptr) {
        removed_voxels->clear();
    }

    VoxelId newly_solid_count = 0;
    for_each_covered_voxel(
        gas_grid,
        atom_view,
        precursor_radius_,
        [&](VoxelId id) {
            const auto previous_count = gas_grid.blocker_count(id);
            if (previous_count == std::numeric_limits<VoxelBlockerCount>::max()) {
                throw std::overflow_error("voxel blocker count overflow");
            }
            if (previous_count == 0) {
                if (removed_voxels != nullptr) {
                    removed_voxels->push_back({id, gas_grid.gas_state(id)});
                }
                ++newly_solid_count;
            }
            gas_grid.set_blocker_count(id, previous_count + 1);
        });
    return newly_solid_count;
}

AtomChangeOccupancyResult AtomVoxelizer::apply_atom_changes(
    GasGrid& gas_grid,
    const AtomChangeBatch& atom_changes) const
{
    std::vector<VoxelOccupancyChange> voxel_changes;
    return apply_atom_changes(gas_grid, atom_changes, voxel_changes);
}

AtomChangeOccupancyResult AtomVoxelizer::apply_atom_changes(
    GasGrid& gas_grid,
    const AtomChangeBatch& atom_changes,
    std::vector<VoxelOccupancyChange>& voxel_changes) const
{
    validate_atom_view(gas_grid, atom_changes.added_atoms, precursor_radius_);
    validate_atom_view(gas_grid, atom_changes.removed_atoms, precursor_radius_);

    std::vector<VoxelDelta> deltas;
    for_each_covered_voxel(
        gas_grid,
        atom_changes.added_atoms,
        precursor_radius_,
        [&](VoxelId id) { deltas.push_back({id, 1}); });
    for_each_covered_voxel(
        gas_grid,
        atom_changes.removed_atoms,
        precursor_radius_,
        [&](VoxelId id) { deltas.push_back({id, -1}); });
    std::sort(
        deltas.begin(),
        deltas.end(),
        [](const VoxelDelta& lhs, const VoxelDelta& rhs) {
            return lhs.voxel_id < rhs.voxel_id;
        });

    std::vector<VoxelOccupancyChange> prepared_changes;
    prepared_changes.reserve(deltas.size());
    AtomChangeOccupancyResult result{};
    for (std::size_t begin = 0; begin < deltas.size();) {
        std::size_t end = begin;
        std::uint64_t addition_count = 0;
        std::uint64_t removal_count = 0;
        while (end < deltas.size()
               && deltas[end].voxel_id == deltas[begin].voxel_id) {
            if (deltas[end].delta > 0) {
                ++addition_count;
            } else {
                ++removal_count;
            }
            ++end;
        }
        const auto id = deltas[begin].voxel_id;
        const auto previous_count = gas_grid.blocker_count(id);
        VoxelBlockerCount updated_count = previous_count;
        if (addition_count >= removal_count) {
            const auto increase = addition_count - removal_count;
            if (increase
                > std::numeric_limits<VoxelBlockerCount>::max()
                    - previous_count) {
                throw std::overflow_error(
                    "atom addition overflows a voxel blocker count");
            }
            updated_count = static_cast<VoxelBlockerCount>(
                previous_count + increase);
        } else {
            const auto decrease = removal_count - addition_count;
            if (decrease > previous_count) {
                throw std::underflow_error(
                    "atom removal underflows a voxel blocker count");
            }
            updated_count = static_cast<VoxelBlockerCount>(
                previous_count - decrease);
        }
        if (updated_count != previous_count) {
            prepared_changes.push_back({
                id,
                previous_count,
                updated_count,
                gas_grid.gas_state(id)});
            ++result.blocker_count_changed_voxel_count;
            if (previous_count == 0) {
                ++result.newly_solid_count;
            } else if (updated_count == 0) {
                ++result.newly_gas_count;
            }
        }
        begin = end;
    }

    for (const auto& change : prepared_changes) {
        gas_grid.set_blocker_count(change.voxel_id, change.blocker_count);
    }
    voxel_changes.swap(prepared_changes);
    return result;
}

}  // namespace gasaccess
