#include "gasaccess/gas_grid.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace gasaccess {
namespace {

VoxelId checked_multiply(VoxelId lhs, VoxelId rhs)
{
    if (lhs != 0 && rhs > std::numeric_limits<VoxelId>::max() / lhs) {
        throw std::overflow_error("voxel count overflows VoxelId");
    }
    return lhs * rhs;
}

void validate_axis_geometry(double origin, double spacing, std::uint64_t dimension)
{
    const double length = spacing * static_cast<double>(dimension);
    if (!std::isfinite(length) || length <= 0.0) {
        throw std::invalid_argument("grid axis has a non-finite physical length");
    }
    if (!std::isfinite(origin + length)) {
        throw std::invalid_argument("grid axis upper boundary is not finite");
    }

    const double first_center = origin + 0.5 * spacing;
    if (first_center == origin) {
        throw std::invalid_argument("grid spacing is too small at the axis origin");
    }
    if (dimension > 1) {
        const double second_center = origin + 1.5 * spacing;
        const double last_center = origin
            + (static_cast<double>(dimension) - 0.5) * spacing;
        const double previous_center = origin
            + (static_cast<double>(dimension) - 1.5) * spacing;
        if (first_center == second_center || previous_center == last_center) {
            throw std::invalid_argument("adjacent voxel centers are not distinguishable");
        }
    }
}

}  // namespace

GasGrid::GasGrid(GridSpec grid_spec)
    : grid_spec_(std::move(grid_spec)),
      voxel_count_(validate_and_count_voxels(grid_spec_))
{
    validate_boundary_conditions(grid_spec_);

    const auto state_count = static_cast<std::size_t>(voxel_count_);
    if (state_count > states_.max_size()) {
        throw std::length_error("voxel count exceeds state container capacity");
    }
    states_.assign(state_count, GasState::Unclassified);
    initialize_explicit_sources();
}

const GridSpec& GasGrid::grid_spec() const noexcept
{
    return grid_spec_;
}

VoxelId GasGrid::voxel_count() const noexcept
{
    return voxel_count_;
}

bool GasGrid::contains(const VoxelCoord& voxel_coord) const noexcept
{
    if (voxel_coord.x < 0 || voxel_coord.y < 0 || voxel_coord.z < 0) {
        return false;
    }

    return static_cast<std::uint64_t>(voxel_coord.x) < grid_spec_.dimensions.x
        && static_cast<std::uint64_t>(voxel_coord.y) < grid_spec_.dimensions.y
        && static_cast<std::uint64_t>(voxel_coord.z) < grid_spec_.dimensions.z;
}

VoxelId GasGrid::voxel_id(const VoxelCoord& voxel_coord) const
{
    if (!contains(voxel_coord)) {
        throw std::out_of_range("voxel coordinate is outside the grid");
    }

    const auto x = static_cast<VoxelId>(voxel_coord.x);
    const auto y = static_cast<VoxelId>(voxel_coord.y);
    const auto z = static_cast<VoxelId>(voxel_coord.z);
    return (z * grid_spec_.dimensions.y + y) * grid_spec_.dimensions.x + x;
}

VoxelCoord GasGrid::voxel_coord(VoxelId voxel_id) const
{
    if (voxel_id >= voxel_count_) {
        throw std::out_of_range("voxel identifier is outside the grid");
    }

    const auto x = voxel_id % grid_spec_.dimensions.x;
    const auto yz = voxel_id / grid_spec_.dimensions.x;
    const auto y = yz % grid_spec_.dimensions.y;
    const auto z = yz / grid_spec_.dimensions.y;

    return {
        static_cast<std::int64_t>(x),
        static_cast<std::int64_t>(y),
        static_cast<std::int64_t>(z)
    };
}

std::optional<VoxelCoord> GasGrid::locate_voxel(const Point3& point) const noexcept
{
    const auto x = locate_axis(
        point.x,
        grid_spec_.origin.x,
        grid_spec_.dimensions.x,
        grid_spec_.periodic.x);
    if (!x) {
        return std::nullopt;
    }

    const auto y = locate_axis(
        point.y,
        grid_spec_.origin.y,
        grid_spec_.dimensions.y,
        grid_spec_.periodic.y);
    if (!y) {
        return std::nullopt;
    }

    const auto z = locate_axis(
        point.z,
        grid_spec_.origin.z,
        grid_spec_.dimensions.z,
        grid_spec_.periodic.z);
    if (!z) {
        return std::nullopt;
    }

    return VoxelCoord{*x, *y, *z};
}

Point3 GasGrid::voxel_center(const VoxelCoord& voxel_coord) const
{
    if (!contains(voxel_coord)) {
        throw std::out_of_range("voxel coordinate is outside the grid");
    }

    return {
        grid_spec_.origin.x
            + (static_cast<double>(voxel_coord.x) + 0.5) * grid_spec_.spacing,
        grid_spec_.origin.y
            + (static_cast<double>(voxel_coord.y) + 0.5) * grid_spec_.spacing,
        grid_spec_.origin.z
            + (static_cast<double>(voxel_coord.z) + 0.5) * grid_spec_.spacing
    };
}

Point3 GasGrid::voxel_center(VoxelId voxel_id) const
{
    return voxel_center(voxel_coord(voxel_id));
}

NeighborList GasGrid::neighbors(VoxelId voxel_id) const
{
    const auto center = voxel_coord(voxel_id);
    const std::array<VoxelCoord, 6> candidates{{
        {center.x - 1, center.y, center.z},
        {center.x + 1, center.y, center.z},
        {center.x, center.y - 1, center.z},
        {center.x, center.y + 1, center.z},
        {center.x, center.y, center.z - 1},
        {center.x, center.y, center.z + 1}
    }};

    NeighborList result{};
    for (auto candidate : candidates) {
        if (candidate.x < 0) {
            if (!grid_spec_.periodic.x) {
                continue;
            }
            candidate.x = static_cast<std::int64_t>(grid_spec_.dimensions.x - 1);
        } else if (static_cast<std::uint64_t>(candidate.x) >= grid_spec_.dimensions.x) {
            if (!grid_spec_.periodic.x) {
                continue;
            }
            candidate.x = 0;
        }

        if (candidate.y < 0) {
            if (!grid_spec_.periodic.y) {
                continue;
            }
            candidate.y = static_cast<std::int64_t>(grid_spec_.dimensions.y - 1);
        } else if (static_cast<std::uint64_t>(candidate.y) >= grid_spec_.dimensions.y) {
            if (!grid_spec_.periodic.y) {
                continue;
            }
            candidate.y = 0;
        }

        if (candidate.z < 0) {
            if (!grid_spec_.periodic.z) {
                continue;
            }
            candidate.z = static_cast<std::int64_t>(grid_spec_.dimensions.z - 1);
        } else if (static_cast<std::uint64_t>(candidate.z) >= grid_spec_.dimensions.z) {
            if (!grid_spec_.periodic.z) {
                continue;
            }
            candidate.z = 0;
        }

        if (candidate == center) {
            continue;
        }
        add_unique_neighbor(result, candidate);
    }
    return result;
}

bool GasGrid::is_reservoir_source(const VoxelCoord& voxel_coord) const
{
    if (!contains(voxel_coord)) {
        throw std::out_of_range("voxel coordinate is outside the grid");
    }

    if (is_boundary_source(voxel_coord)) {
        return true;
    }

    const auto id = voxel_id(voxel_coord);
    return std::binary_search(
        explicit_source_ids_.begin(),
        explicit_source_ids_.end(),
        id);
}

bool GasGrid::is_reservoir_source(VoxelId voxel_id) const
{
    const auto coordinate = voxel_coord(voxel_id);
    if (is_boundary_source(coordinate)) {
        return true;
    }
    return std::binary_search(
        explicit_source_ids_.begin(),
        explicit_source_ids_.end(),
        voxel_id);
}

const std::vector<VoxelId>& GasGrid::explicit_source_ids() const noexcept
{
    return explicit_source_ids_;
}

GasState GasGrid::gas_state(VoxelId voxel_id) const
{
    if (voxel_id >= voxel_count_) {
        throw std::out_of_range("voxel identifier is outside the grid");
    }
    return states_[static_cast<std::size_t>(voxel_id)];
}

GasState GasGrid::gas_state(const VoxelCoord& voxel_coord_value) const
{
    return gas_state(voxel_id(voxel_coord_value));
}

void GasGrid::set_gas_state(VoxelId voxel_id, GasState gas_state)
{
    if (voxel_id >= voxel_count_) {
        throw std::out_of_range("voxel identifier is outside the grid");
    }
    states_[static_cast<std::size_t>(voxel_id)] = gas_state;
}

void GasGrid::fill_gas_state(GasState gas_state) noexcept
{
    std::fill(states_.begin(), states_.end(), gas_state);
}

VoxelId GasGrid::validate_and_count_voxels(const GridSpec& grid_spec)
{
    if (!std::isfinite(grid_spec.spacing) || grid_spec.spacing <= 0.0) {
        throw std::invalid_argument("grid spacing must be finite and positive");
    }
    if (!std::isfinite(grid_spec.origin.x)
        || !std::isfinite(grid_spec.origin.y)
        || !std::isfinite(grid_spec.origin.z)) {
        throw std::invalid_argument("grid origin must be finite");
    }
    if (grid_spec.dimensions.x == 0
        || grid_spec.dimensions.y == 0
        || grid_spec.dimensions.z == 0) {
        throw std::invalid_argument("grid dimensions must be nonzero");
    }

    const auto max_coordinate = static_cast<std::uint64_t>(
        std::numeric_limits<std::int64_t>::max());
    if (grid_spec.dimensions.x > max_coordinate
        || grid_spec.dimensions.y > max_coordinate
        || grid_spec.dimensions.z > max_coordinate) {
        throw std::overflow_error("grid dimension exceeds VoxelCoord range");
    }

    validate_axis_geometry(
        grid_spec.origin.x,
        grid_spec.spacing,
        grid_spec.dimensions.x);
    validate_axis_geometry(
        grid_spec.origin.y,
        grid_spec.spacing,
        grid_spec.dimensions.y);
    validate_axis_geometry(
        grid_spec.origin.z,
        grid_spec.spacing,
        grid_spec.dimensions.z);

    const auto xy = checked_multiply(grid_spec.dimensions.x, grid_spec.dimensions.y);
    const auto xyz = checked_multiply(xy, grid_spec.dimensions.z);
    if (xyz > static_cast<VoxelId>(std::numeric_limits<std::size_t>::max())) {
        throw std::length_error("voxel count exceeds local addressable memory range");
    }
    return xyz;
}

void GasGrid::validate_boundary_conditions(const GridSpec& grid_spec)
{
    if (grid_spec.periodic.x
        && (grid_spec.reservoir_faces.x_low || grid_spec.reservoir_faces.x_high)) {
        throw std::invalid_argument("periodic x faces cannot be reservoir faces");
    }
    if (grid_spec.periodic.y
        && (grid_spec.reservoir_faces.y_low || grid_spec.reservoir_faces.y_high)) {
        throw std::invalid_argument("periodic y faces cannot be reservoir faces");
    }
    if (grid_spec.periodic.z
        && (grid_spec.reservoir_faces.z_low || grid_spec.reservoir_faces.z_high)) {
        throw std::invalid_argument("periodic z faces cannot be reservoir faces");
    }
}

std::optional<std::int64_t> GasGrid::locate_axis(
    double position,
    double origin,
    std::uint64_t dimension,
    bool periodic) const noexcept
{
    if (!std::isfinite(position)) {
        return std::nullopt;
    }

    const double length = grid_spec_.spacing * static_cast<double>(dimension);
    double relative = position - origin;
    if (!std::isfinite(relative)) {
        return std::nullopt;
    }

    if (periodic) {
        relative = std::fmod(relative, length);
        if (!std::isfinite(relative)) {
            return std::nullopt;
        }
        if (relative < 0.0) {
            relative += length;
        }
        if (relative >= length) {
            relative = 0.0;
        }
    } else if (relative < 0.0 || relative >= length) {
        return std::nullopt;
    }

    const double raw_index = std::floor(relative / grid_spec_.spacing);
    if (raw_index < 0.0 || raw_index >= static_cast<double>(dimension)) {
        return std::nullopt;
    }
    return static_cast<std::int64_t>(raw_index);
}

bool GasGrid::is_boundary_source(const VoxelCoord& voxel_coord) const noexcept
{
    const auto& faces = grid_spec_.reservoir_faces;
    const auto& dimensions = grid_spec_.dimensions;

    return (faces.x_low && voxel_coord.x == 0)
        || (faces.x_high
            && static_cast<std::uint64_t>(voxel_coord.x) == dimensions.x - 1)
        || (faces.y_low && voxel_coord.y == 0)
        || (faces.y_high
            && static_cast<std::uint64_t>(voxel_coord.y) == dimensions.y - 1)
        || (faces.z_low && voxel_coord.z == 0)
        || (faces.z_high
            && static_cast<std::uint64_t>(voxel_coord.z) == dimensions.z - 1);
}

void GasGrid::initialize_explicit_sources()
{
    explicit_source_ids_.reserve(grid_spec_.explicit_source_voxels.size());
    for (const auto& voxel_coord_value : grid_spec_.explicit_source_voxels) {
        if (!contains(voxel_coord_value)) {
            throw std::invalid_argument("explicit source voxel is outside the grid");
        }
        explicit_source_ids_.push_back(voxel_id(voxel_coord_value));
    }

    std::sort(explicit_source_ids_.begin(), explicit_source_ids_.end());
    explicit_source_ids_.erase(
        std::unique(explicit_source_ids_.begin(), explicit_source_ids_.end()),
        explicit_source_ids_.end());
}

void GasGrid::add_unique_neighbor(
    NeighborList& neighbor_list,
    const VoxelCoord& voxel_coord_value) const
{
    const auto id = voxel_id(voxel_coord_value);
    for (std::size_t index = 0; index < neighbor_list.count; ++index) {
        if (neighbor_list.ids[index] == id) {
            return;
        }
    }

    if (neighbor_list.count < neighbor_list.ids.size()) {
        neighbor_list.ids[neighbor_list.count] = id;
        ++neighbor_list.count;
    }
}

}  // namespace gasaccess
