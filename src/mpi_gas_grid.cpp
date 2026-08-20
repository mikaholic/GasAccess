#include "gasaccess/mpi_gas_grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace gasaccess {
namespace {

constexpr int halo_tag_base = 4210;

std::uint64_t checked_multiply(std::uint64_t lhs, std::uint64_t rhs)
{
    if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
        throw std::overflow_error("MPI voxel extent multiplication overflow");
    }
    return lhs * rhs;
}

std::uint64_t checked_add(std::uint64_t lhs, std::uint64_t rhs)
{
    if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
        throw std::overflow_error("MPI voxel extent addition overflow");
    }
    return lhs + rhs;
}

bool nearly_equal(double lhs, double rhs) noexcept
{
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * scale;
    return std::abs(lhs - rhs) <= tolerance;
}

void require_finite_point(const Point3& point, const char* field_name)
{
    if (!std::isfinite(point.x)
        || !std::isfinite(point.y)
        || !std::isfinite(point.z)) {
        throw std::invalid_argument(std::string(field_name) + " must be finite");
    }
}

void require_positive_spacing(const GridSpacing& spacing)
{
    if (!std::isfinite(spacing.x) || spacing.x <= 0.0
        || !std::isfinite(spacing.y) || spacing.y <= 0.0
        || !std::isfinite(spacing.z) || spacing.z <= 0.0) {
        throw std::invalid_argument(
            "all requested grid spacing components must be finite and positive");
    }
}

void require_process_grid(const GridDimensions& process_grid)
{
    if (process_grid.x == 0 || process_grid.y == 0 || process_grid.z == 0) {
        throw std::invalid_argument("all process-grid dimensions must be nonzero");
    }
}

std::uint64_t choose_aligned_dimension(
    double length,
    std::uint64_t process_count,
    double requested_spacing)
{
    const long double ideal_block_count = static_cast<long double>(length)
        / (static_cast<long double>(requested_spacing)
            * static_cast<long double>(process_count));
    if (!std::isfinite(ideal_block_count) || ideal_block_count <= 0.0L) {
        throw std::invalid_argument("requested spacing cannot resolve the box axis");
    }

    const auto maximum_coordinate = static_cast<std::uint64_t>(
        std::numeric_limits<std::int64_t>::max());
    const long double maximum_block_count = static_cast<long double>(
        maximum_coordinate / process_count);
    if (ideal_block_count > maximum_block_count) {
        throw std::overflow_error("aligned voxel dimension exceeds VoxelCoord range");
    }

    const auto lower_block_count = ideal_block_count < 1.0L
        ? std::uint64_t{1}
        : static_cast<std::uint64_t>(std::floor(ideal_block_count));
    const auto upper_block_count = ideal_block_count <= 1.0L
        ? std::uint64_t{1}
        : static_cast<std::uint64_t>(std::ceil(ideal_block_count));

    const auto lower_dimension = checked_multiply(lower_block_count, process_count);
    const auto upper_dimension = checked_multiply(upper_block_count, process_count);
    const double lower_spacing = length / static_cast<double>(lower_dimension);
    const double upper_spacing = length / static_cast<double>(upper_dimension);
    const double lower_error = std::abs(lower_spacing - requested_spacing);
    const double upper_error = std::abs(upper_spacing - requested_spacing);
    return upper_error <= lower_error ? upper_dimension : lower_dimension;
}

std::uint64_t rank_for_location(
    const VoxelCoord& process_location,
    const GridDimensions& process_grid)
{
    const auto x = static_cast<std::uint64_t>(process_location.x);
    const auto y = static_cast<std::uint64_t>(process_location.y);
    const auto z = static_cast<std::uint64_t>(process_location.z);
    return (z * process_grid.y + y) * process_grid.x + x;
}

void check_mpi(int error_code, const char* operation)
{
    if (error_code == MPI_SUCCESS) {
        return;
    }

    std::array<char, MPI_MAX_ERROR_STRING> error_text{};
    int error_length = 0;
    static_cast<void>(MPI_Error_string(
        error_code,
        error_text.data(),
        &error_length));
    throw std::runtime_error(
        std::string(operation) + " failed: "
        + std::string(error_text.data(), static_cast<std::size_t>(error_length)));
}

void require_mpi_active()
{
    int initialized = 0;
    int finalized = 0;
    check_mpi(MPI_Initialized(&initialized), "MPI_Initialized");
    check_mpi(MPI_Finalized(&finalized), "MPI_Finalized");
    if (initialized == 0 || finalized != 0) {
        throw std::logic_error(
            "MPI must be initialized and not finalized before creating an MPI grid");
    }
}

bool is_valid_gas_state(GasState gas_state) noexcept
{
    return gas_state == GasState::Unclassified
        || gas_state == GasState::Solid
        || gas_state == GasState::OutsideAccessible
        || gas_state == GasState::ClosedVoid;
}

std::optional<std::int64_t> locate_axis(
    double position,
    double origin,
    double spacing,
    std::uint64_t dimension,
    bool periodic) noexcept
{
    if (!std::isfinite(position)) {
        return std::nullopt;
    }

    const double length = spacing * static_cast<double>(dimension);
    double relative = position - origin;
    if (!std::isfinite(relative)) {
        return std::nullopt;
    }
    if (periodic) {
        relative = std::fmod(relative, length);
        if (relative < 0.0) {
            relative += length;
        }
        if (relative >= length) {
            relative = 0.0;
        }
    } else if (relative < 0.0 || relative >= length) {
        return std::nullopt;
    }

    const double raw_index = std::floor(relative / spacing);
    if (raw_index < 0.0 || raw_index >= static_cast<double>(dimension)) {
        return std::nullopt;
    }
    return static_cast<std::int64_t>(raw_index);
}

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
    return periodic
        ? std::min(direct_distance, std::abs(length - direct_distance))
        : direct_distance;
}

std::int64_t checked_int64(std::uint64_t value, const char* field_name)
{
    const auto maximum = static_cast<std::uint64_t>(
        std::numeric_limits<std::int64_t>::max());
    if (value > maximum) {
        throw std::overflow_error(std::string(field_name) + " exceeds int64 range");
    }
    return static_cast<std::int64_t>(value);
}

}  // namespace

AlignedGridGeometry make_aligned_grid_geometry(
    const Point3& global_lower,
    const Point3& global_upper,
    const GridDimensions& process_grid,
    const GridSpacing& requested_spacing)
{
    require_finite_point(global_lower, "global lower bound");
    require_finite_point(global_upper, "global upper bound");
    require_process_grid(process_grid);
    require_positive_spacing(requested_spacing);

    const Point3 lengths{
        global_upper.x - global_lower.x,
        global_upper.y - global_lower.y,
        global_upper.z - global_lower.z
    };
    if (!std::isfinite(lengths.x) || lengths.x <= 0.0
        || !std::isfinite(lengths.y) || lengths.y <= 0.0
        || !std::isfinite(lengths.z) || lengths.z <= 0.0) {
        throw std::invalid_argument(
            "global upper bounds must be finite and greater than lower bounds");
    }

    AlignedGridGeometry geometry{};
    geometry.dimensions = {
        choose_aligned_dimension(
            lengths.x,
            process_grid.x,
            requested_spacing.x),
        choose_aligned_dimension(
            lengths.y,
            process_grid.y,
            requested_spacing.y),
        choose_aligned_dimension(
            lengths.z,
            process_grid.z,
            requested_spacing.z)
    };
    geometry.spacing = {
        lengths.x / static_cast<double>(geometry.dimensions.x),
        lengths.y / static_cast<double>(geometry.dimensions.y),
        lengths.z / static_cast<double>(geometry.dimensions.z)
    };
    return geometry;
}

void validate_atom_ghost_coverage(
    double atom_ghost_distance,
    double maximum_excluded_radius)
{
    if (!std::isfinite(atom_ghost_distance) || atom_ghost_distance < 0.0) {
        throw std::invalid_argument(
            "atom ghost distance must be finite and nonnegative");
    }
    if (!std::isfinite(maximum_excluded_radius)
        || maximum_excluded_radius < 0.0) {
        throw std::invalid_argument(
            "maximum excluded radius must be finite and nonnegative");
    }
    if (atom_ghost_distance < maximum_excluded_radius
        && !nearly_equal(atom_ghost_distance, maximum_excluded_radius)) {
        throw std::invalid_argument(
            "atom ghost distance is smaller than the maximum excluded radius");
    }
}

MpiDecomposition::MpiDecomposition(
    const GridSpec& global_grid_spec,
    MpiDecompositionSpec decomposition_spec)
    : spec_(std::move(decomposition_spec))
{
    require_mpi_active();
    if (spec_.communicator == MPI_COMM_NULL) {
        throw std::invalid_argument("MPI communicator must not be MPI_COMM_NULL");
    }
    check_mpi(MPI_Comm_rank(spec_.communicator, &rank_), "MPI_Comm_rank");
    check_mpi(MPI_Comm_size(spec_.communicator, &size_), "MPI_Comm_size");

    require_finite_point(spec_.global_lower, "global lower bound");
    require_finite_point(spec_.global_upper, "global upper bound");
    require_finite_point(spec_.local_lower, "local lower bound");
    require_finite_point(spec_.local_upper, "local upper bound");
    require_process_grid(spec_.process_grid);
    require_positive_spacing(global_grid_spec.spacing);
    validate_atom_ghost_coverage(
        spec_.atom_ghost_distance,
        spec_.maximum_excluded_radius);

    const auto process_count = checked_multiply(
        checked_multiply(spec_.process_grid.x, spec_.process_grid.y),
        spec_.process_grid.z);
    if (process_count != static_cast<std::uint64_t>(size_)) {
        throw std::invalid_argument(
            "process-grid dimensions do not match MPI communicator size");
    }
    if (spec_.process_location.x < 0
        || spec_.process_location.y < 0
        || spec_.process_location.z < 0
        || static_cast<std::uint64_t>(spec_.process_location.x)
            >= spec_.process_grid.x
        || static_cast<std::uint64_t>(spec_.process_location.y)
            >= spec_.process_grid.y
        || static_cast<std::uint64_t>(spec_.process_location.z)
            >= spec_.process_grid.z) {
        throw std::invalid_argument("MPI process location is outside the process grid");
    }
    if (rank_for_location(spec_.process_location, spec_.process_grid)
        != static_cast<std::uint64_t>(rank_)) {
        throw std::invalid_argument(
            "MPI rank does not match the SPPARKS x-fastest process-grid ordering");
    }

    if (global_grid_spec.dimensions.x == 0
        || global_grid_spec.dimensions.y == 0
        || global_grid_spec.dimensions.z == 0) {
        throw std::invalid_argument("all global voxel dimensions must be nonzero");
    }
    const auto maximum_coordinate = static_cast<std::uint64_t>(
        std::numeric_limits<std::int64_t>::max());
    if (global_grid_spec.dimensions.x > maximum_coordinate
        || global_grid_spec.dimensions.y > maximum_coordinate
        || global_grid_spec.dimensions.z > maximum_coordinate) {
        throw std::overflow_error("global voxel dimension exceeds VoxelCoord range");
    }
    if ((global_grid_spec.periodic.x
            && (global_grid_spec.reservoir_faces.x_low
                || global_grid_spec.reservoir_faces.x_high))
        || (global_grid_spec.periodic.y
            && (global_grid_spec.reservoir_faces.y_low
                || global_grid_spec.reservoir_faces.y_high))
        || (global_grid_spec.periodic.z
            && (global_grid_spec.reservoir_faces.z_low
                || global_grid_spec.reservoir_faces.z_high))) {
        throw std::invalid_argument(
            "periodic GasAccess faces cannot be reservoir faces");
    }
    for (const auto& source : global_grid_spec.explicit_source_voxels) {
        if (source.x < 0 || source.y < 0 || source.z < 0
            || static_cast<std::uint64_t>(source.x)
                >= global_grid_spec.dimensions.x
            || static_cast<std::uint64_t>(source.y)
                >= global_grid_spec.dimensions.y
            || static_cast<std::uint64_t>(source.z)
                >= global_grid_spec.dimensions.z) {
            throw std::invalid_argument("explicit source voxel is outside the grid");
        }
    }
    if (global_grid_spec.dimensions.x % spec_.process_grid.x != 0
        || global_grid_spec.dimensions.y % spec_.process_grid.y != 0
        || global_grid_spec.dimensions.z % spec_.process_grid.z != 0) {
        throw std::invalid_argument(
            "global voxel dimensions must be divisible by process-grid dimensions");
    }
    if (!nearly_equal(global_grid_spec.origin.x, spec_.global_lower.x)
        || !nearly_equal(global_grid_spec.origin.y, spec_.global_lower.y)
        || !nearly_equal(global_grid_spec.origin.z, spec_.global_lower.z)) {
        throw std::invalid_argument("grid origin does not match the global lower bound");
    }

    const Point3 grid_upper{
        global_grid_spec.origin.x
            + global_grid_spec.spacing.x
                * static_cast<double>(global_grid_spec.dimensions.x),
        global_grid_spec.origin.y
            + global_grid_spec.spacing.y
                * static_cast<double>(global_grid_spec.dimensions.y),
        global_grid_spec.origin.z
            + global_grid_spec.spacing.z
                * static_cast<double>(global_grid_spec.dimensions.z)
    };
    if (!std::isfinite(grid_upper.x)
        || !std::isfinite(grid_upper.y)
        || !std::isfinite(grid_upper.z)
        || !nearly_equal(grid_upper.x, spec_.global_upper.x)
        || !nearly_equal(grid_upper.y, spec_.global_upper.y)
        || !nearly_equal(grid_upper.z, spec_.global_upper.z)) {
        throw std::invalid_argument("grid extent does not match the global box");
    }

    owned_range_.dimensions = {
        global_grid_spec.dimensions.x / spec_.process_grid.x,
        global_grid_spec.dimensions.y / spec_.process_grid.y,
        global_grid_spec.dimensions.z / spec_.process_grid.z
    };
    const auto begin_x = checked_multiply(
        static_cast<std::uint64_t>(spec_.process_location.x),
        owned_range_.dimensions.x);
    const auto begin_y = checked_multiply(
        static_cast<std::uint64_t>(spec_.process_location.y),
        owned_range_.dimensions.y);
    const auto begin_z = checked_multiply(
        static_cast<std::uint64_t>(spec_.process_location.z),
        owned_range_.dimensions.z);
    const auto end_x = checked_add(begin_x, owned_range_.dimensions.x);
    const auto end_y = checked_add(begin_y, owned_range_.dimensions.y);
    const auto end_z = checked_add(begin_z, owned_range_.dimensions.z);
    owned_range_.begin = {
        checked_int64(begin_x, "owned x begin"),
        checked_int64(begin_y, "owned y begin"),
        checked_int64(begin_z, "owned z begin")
    };
    owned_range_.end = {
        checked_int64(end_x, "owned x end"),
        checked_int64(end_y, "owned y end"),
        checked_int64(end_z, "owned z end")
    };

    const Point3 expected_local_lower{
        global_grid_spec.origin.x
            + static_cast<double>(begin_x) * global_grid_spec.spacing.x,
        global_grid_spec.origin.y
            + static_cast<double>(begin_y) * global_grid_spec.spacing.y,
        global_grid_spec.origin.z
            + static_cast<double>(begin_z) * global_grid_spec.spacing.z
    };
    const Point3 expected_local_upper{
        global_grid_spec.origin.x
            + static_cast<double>(end_x) * global_grid_spec.spacing.x,
        global_grid_spec.origin.y
            + static_cast<double>(end_y) * global_grid_spec.spacing.y,
        global_grid_spec.origin.z
            + static_cast<double>(end_z) * global_grid_spec.spacing.z
    };
    if (!nearly_equal(spec_.local_lower.x, expected_local_lower.x)
        || !nearly_equal(spec_.local_lower.y, expected_local_lower.y)
        || !nearly_equal(spec_.local_lower.z, expected_local_lower.z)
        || !nearly_equal(spec_.local_upper.x, expected_local_upper.x)
        || !nearly_equal(spec_.local_upper.y, expected_local_upper.y)
        || !nearly_equal(spec_.local_upper.z, expected_local_upper.z)) {
        throw std::invalid_argument(
            "SPPARKS subdomain faces do not align with global voxel faces");
    }

    const auto neighbor_for = [&](Face face) {
        VoxelCoord location = spec_.process_location;
        const bool low = face == Face::XLow
            || face == Face::YLow
            || face == Face::ZLow;
        std::int64_t* component = nullptr;
        std::uint64_t dimension = 0;
        bool periodic = false;
        switch (face) {
        case Face::XLow:
        case Face::XHigh:
            component = &location.x;
            dimension = spec_.process_grid.x;
            periodic = global_grid_spec.periodic.x;
            break;
        case Face::YLow:
        case Face::YHigh:
            component = &location.y;
            dimension = spec_.process_grid.y;
            periodic = global_grid_spec.periodic.y;
            break;
        case Face::ZLow:
        case Face::ZHigh:
            component = &location.z;
            dimension = spec_.process_grid.z;
            periodic = global_grid_spec.periodic.z;
            break;
        }

        *component += low ? -1 : 1;
        if (*component < 0) {
            if (!periodic) {
                return MPI_PROC_NULL;
            }
            *component = static_cast<std::int64_t>(dimension - 1);
        } else if (static_cast<std::uint64_t>(*component) >= dimension) {
            if (!periodic) {
                return MPI_PROC_NULL;
            }
            *component = 0;
        }
        return static_cast<int>(rank_for_location(location, spec_.process_grid));
    };

    for (std::size_t index = 0; index < neighbor_ranks_.size(); ++index) {
        neighbor_ranks_[index] = neighbor_for(static_cast<Face>(index));
    }
}

const MpiDecompositionSpec& MpiDecomposition::spec() const noexcept
{
    return spec_;
}

const OwnedVoxelRange& MpiDecomposition::owned_range() const noexcept
{
    return owned_range_;
}

int MpiDecomposition::rank() const noexcept
{
    return rank_;
}

int MpiDecomposition::size() const noexcept
{
    return size_;
}

int MpiDecomposition::neighbor_rank(Face face) const noexcept
{
    return neighbor_ranks_[face_index(face)];
}

bool MpiDecomposition::owns(const VoxelCoord& global_voxel_coord) const noexcept
{
    return global_voxel_coord.x >= owned_range_.begin.x
        && global_voxel_coord.x < owned_range_.end.x
        && global_voxel_coord.y >= owned_range_.begin.y
        && global_voxel_coord.y < owned_range_.end.y
        && global_voxel_coord.z >= owned_range_.begin.z
        && global_voxel_coord.z < owned_range_.end.z;
}

std::size_t MpiDecomposition::face_index(Face face) noexcept
{
    return static_cast<std::size_t>(face);
}

DistributedGasGrid::DistributedGasGrid(
    GridSpec global_grid_spec,
    MpiDecompositionSpec decomposition_spec)
    : global_grid_spec_(std::move(global_grid_spec)),
      decomposition_(global_grid_spec_, std::move(decomposition_spec))
{
    const auto& owned_dimensions = decomposition_.owned_range().dimensions;
    storage_dimensions_ = {
        checked_add(owned_dimensions.x, 2),
        checked_add(owned_dimensions.y, 2),
        checked_add(owned_dimensions.z, 2)
    };
    owned_voxel_count_ = checked_multiply(
        checked_multiply(owned_dimensions.x, owned_dimensions.y),
        owned_dimensions.z);
    const auto storage_count = checked_multiply(
        checked_multiply(storage_dimensions_.x, storage_dimensions_.y),
        storage_dimensions_.z);
    if (storage_count > static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max())) {
        throw std::length_error("local MPI voxel storage exceeds size_t range");
    }
    states_.assign(
        static_cast<std::size_t>(storage_count),
        GasState::Unclassified);
    owned_state_counts_[gas_state_index(GasState::Unclassified)] =
        owned_voxel_count_;

    for (std::size_t index = 0; index < send_buffers_.size(); ++index) {
        const auto count = face_element_count(static_cast<Face>(index));
        if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::length_error("gas-state face exceeds MPI int count range");
        }
        send_buffers_[index].resize(count);
        receive_buffers_[index].resize(count);
    }
}

const GridSpec& DistributedGasGrid::global_grid_spec() const noexcept
{
    return global_grid_spec_;
}

const MpiDecomposition& DistributedGasGrid::decomposition() const noexcept
{
    return decomposition_;
}

const OwnedVoxelRange& DistributedGasGrid::owned_range() const noexcept
{
    return decomposition_.owned_range();
}

std::uint64_t DistributedGasGrid::owned_voxel_count() const noexcept
{
    return owned_voxel_count_;
}

bool DistributedGasGrid::owns(const VoxelCoord& global_voxel_coord) const noexcept
{
    return decomposition_.owns(global_voxel_coord);
}

std::optional<VoxelCoord> DistributedGasGrid::locate_voxel(
    const Point3& point) const noexcept
{
    const auto x = locate_axis(
        point.x,
        global_grid_spec_.origin.x,
        global_grid_spec_.spacing.x,
        global_grid_spec_.dimensions.x,
        global_grid_spec_.periodic.x);
    const auto y = locate_axis(
        point.y,
        global_grid_spec_.origin.y,
        global_grid_spec_.spacing.y,
        global_grid_spec_.dimensions.y,
        global_grid_spec_.periodic.y);
    const auto z = locate_axis(
        point.z,
        global_grid_spec_.origin.z,
        global_grid_spec_.spacing.z,
        global_grid_spec_.dimensions.z,
        global_grid_spec_.periodic.z);
    if (!x || !y || !z) {
        return std::nullopt;
    }
    return VoxelCoord{*x, *y, *z};
}

Point3 DistributedGasGrid::voxel_center(
    const VoxelCoord& global_voxel_coord) const
{
    if (global_voxel_coord.x < 0
        || global_voxel_coord.y < 0
        || global_voxel_coord.z < 0
        || static_cast<std::uint64_t>(global_voxel_coord.x)
            >= global_grid_spec_.dimensions.x
        || static_cast<std::uint64_t>(global_voxel_coord.y)
            >= global_grid_spec_.dimensions.y
        || static_cast<std::uint64_t>(global_voxel_coord.z)
            >= global_grid_spec_.dimensions.z) {
        throw std::out_of_range("global voxel coordinate is outside the grid");
    }
    return {
        global_grid_spec_.origin.x
            + (static_cast<double>(global_voxel_coord.x) + 0.5)
                * global_grid_spec_.spacing.x,
        global_grid_spec_.origin.y
            + (static_cast<double>(global_voxel_coord.y) + 0.5)
                * global_grid_spec_.spacing.y,
        global_grid_spec_.origin.z
            + (static_cast<double>(global_voxel_coord.z) + 0.5)
                * global_grid_spec_.spacing.z
    };
}

GasState DistributedGasGrid::gas_state(
    const VoxelCoord& global_voxel_coord) const
{
    const auto coordinate = local_coord(global_voxel_coord);
    if (!coordinate) {
        throw std::out_of_range(
            "voxel is neither owned nor present in a face ghost layer");
    }
    return states_[state_index(*coordinate)];
}

std::uint64_t DistributedGasGrid::owned_gas_state_count(
    GasState gas_state_value) const
{
    if (!is_valid_gas_state(gas_state_value)) {
        throw std::invalid_argument("invalid gas state value");
    }
    return owned_state_counts_[gas_state_index(gas_state_value)];
}

void DistributedGasGrid::set_owned_gas_state(
    const VoxelCoord& global_voxel_coord,
    GasState gas_state_value)
{
    if (!is_valid_gas_state(gas_state_value)) {
        throw std::invalid_argument("invalid gas state value");
    }
    if (!owns(global_voxel_coord)) {
        throw std::out_of_range("cannot set a gas state for a non-owned voxel");
    }
    const auto coordinate = local_coord(global_voxel_coord);
    const auto index = state_index(*coordinate);
    const auto previous_state = states_[index];
    if (previous_state == gas_state_value) {
        return;
    }
    --owned_state_counts_[gas_state_index(previous_state)];
    ++owned_state_counts_[gas_state_index(gas_state_value)];
    states_[index] = gas_state_value;
}

void DistributedGasGrid::fill_owned_gas_state(GasState gas_state_value)
{
    if (!is_valid_gas_state(gas_state_value)) {
        throw std::invalid_argument("invalid gas state value");
    }
    owned_state_counts_.fill(0);
    owned_state_counts_[gas_state_index(gas_state_value)] = owned_voxel_count_;
    const auto& dimensions = owned_range().dimensions;
    for (std::uint64_t z = 1; z <= dimensions.z; ++z) {
        for (std::uint64_t y = 1; y <= dimensions.y; ++y) {
            for (std::uint64_t x = 1; x <= dimensions.x; ++x) {
                states_[state_index({x, y, z})] = gas_state_value;
            }
        }
    }
}

std::uint64_t DistributedGasGrid::voxelize_owned_atoms(
    AtomView atom_view,
    double precursor_radius)
{
    return voxelize_owned_atoms_impl(atom_view, precursor_radius, nullptr);
}

std::uint64_t DistributedGasGrid::voxelize_owned_atoms(
    AtomView atom_view,
    double precursor_radius,
    std::vector<DistributedRemovedVoxel>& removed_voxels)
{
    return voxelize_owned_atoms_impl(
        atom_view,
        precursor_radius,
        &removed_voxels);
}

std::uint64_t DistributedGasGrid::voxelize_owned_atoms_impl(
    AtomView atom_view,
    double precursor_radius,
    std::vector<DistributedRemovedVoxel>* removed_voxels)
{
    if (atom_view.count != 0 && atom_view.atoms == nullptr) {
        throw std::invalid_argument("atom view has a null pointer with nonzero count");
    }
    if (!std::isfinite(precursor_radius) || precursor_radius < 0.0) {
        throw std::invalid_argument("precursor radius must be finite and nonnegative");
    }

    for (std::size_t atom_index_value = 0;
         atom_index_value < atom_view.count;
         ++atom_index_value) {
        const auto& atom = atom_view.atoms[atom_index_value];
        require_finite_point(atom.position, "atom position");
        if (!std::isfinite(atom.radius) || atom.radius < 0.0) {
            throw std::invalid_argument("atom radius must be finite and nonnegative");
        }
        const double excluded_radius = atom.radius + precursor_radius;
        if (!std::isfinite(excluded_radius)
            || !std::isfinite(excluded_radius * excluded_radius)) {
            throw std::invalid_argument("excluded radius is not finite");
        }
        validate_atom_ghost_coverage(
            decomposition_.spec().atom_ghost_distance,
            excluded_radius);
        if (excluded_radius > decomposition_.spec().maximum_excluded_radius
            && !nearly_equal(
                excluded_radius,
                decomposition_.spec().maximum_excluded_radius)) {
            throw std::invalid_argument(
                "atom excluded radius exceeds the declared maximum");
        }
    }

    if (removed_voxels != nullptr) {
        removed_voxels->clear();
    }

    const Point3 lengths{
        global_grid_spec_.spacing.x
            * static_cast<double>(global_grid_spec_.dimensions.x),
        global_grid_spec_.spacing.y
            * static_cast<double>(global_grid_spec_.dimensions.y),
        global_grid_spec_.spacing.z
            * static_cast<double>(global_grid_spec_.dimensions.z)
    };
    std::uint64_t newly_solid_count = 0;
    for (std::size_t atom_index_value = 0;
         atom_index_value < atom_view.count;
         ++atom_index_value) {
        const auto& atom = atom_view.atoms[atom_index_value];
        Point3 atom_position = atom.position;
        if (global_grid_spec_.periodic.x) {
            atom_position.x = wrap_position(
                atom_position.x,
                global_grid_spec_.origin.x,
                lengths.x);
        }
        if (global_grid_spec_.periodic.y) {
            atom_position.y = wrap_position(
                atom_position.y,
                global_grid_spec_.origin.y,
                lengths.y);
        }
        if (global_grid_spec_.periodic.z) {
            atom_position.z = wrap_position(
                atom_position.z,
                global_grid_spec_.origin.z,
                lengths.z);
        }
        const auto atom_voxel = locate_voxel(atom_position);
        if (!atom_voxel) {
            continue;
        }

        const double excluded_radius = atom.radius + precursor_radius;
        const double excluded_radius_squared = excluded_radius * excluded_radius;
        const auto x_candidates = make_axis_candidates(
            static_cast<std::uint64_t>(atom_voxel->x),
            global_grid_spec_.dimensions.x,
            excluded_radius,
            global_grid_spec_.spacing.x,
            global_grid_spec_.periodic.x);
        const auto y_candidates = make_axis_candidates(
            static_cast<std::uint64_t>(atom_voxel->y),
            global_grid_spec_.dimensions.y,
            excluded_radius,
            global_grid_spec_.spacing.y,
            global_grid_spec_.periodic.y);
        const auto z_candidates = make_axis_candidates(
            static_cast<std::uint64_t>(atom_voxel->z),
            global_grid_spec_.dimensions.z,
            excluded_radius,
            global_grid_spec_.spacing.z,
            global_grid_spec_.periodic.z);

        for (std::uint64_t z_offset = 0;
             z_offset < z_candidates.count;
             ++z_offset) {
            const auto z = axis_index(z_candidates, z_offset);
            const double z_center = global_grid_spec_.origin.z
                + (static_cast<double>(z) + 0.5) * global_grid_spec_.spacing.z;
            const double z_distance = axis_distance(
                z_center,
                atom_position.z,
                lengths.z,
                global_grid_spec_.periodic.z);
            if (z_distance > excluded_radius) {
                continue;
            }
            const double z_distance_squared = z_distance * z_distance;

            for (std::uint64_t y_offset = 0;
                 y_offset < y_candidates.count;
                 ++y_offset) {
                const auto y = axis_index(y_candidates, y_offset);
                const double y_center = global_grid_spec_.origin.y
                    + (static_cast<double>(y) + 0.5)
                        * global_grid_spec_.spacing.y;
                const double y_distance = axis_distance(
                    y_center,
                    atom_position.y,
                    lengths.y,
                    global_grid_spec_.periodic.y);
                const double yz_distance_squared = y_distance * y_distance
                    + z_distance_squared;
                if (y_distance > excluded_radius
                    || yz_distance_squared > excluded_radius_squared) {
                    continue;
                }

                for (std::uint64_t x_offset = 0;
                     x_offset < x_candidates.count;
                     ++x_offset) {
                    const auto x = axis_index(x_candidates, x_offset);
                    const VoxelCoord coordinate{
                        static_cast<std::int64_t>(x),
                        static_cast<std::int64_t>(y),
                        static_cast<std::int64_t>(z)
                    };
                    if (!owns(coordinate)) {
                        continue;
                    }
                    const double x_center = global_grid_spec_.origin.x
                        + (static_cast<double>(x) + 0.5)
                            * global_grid_spec_.spacing.x;
                    const double x_distance = axis_distance(
                        x_center,
                        atom_position.x,
                        lengths.x,
                        global_grid_spec_.periodic.x);
                    const double distance_squared = x_distance * x_distance
                        + yz_distance_squared;
                    if (x_distance > excluded_radius
                        || distance_squared > excluded_radius_squared) {
                        continue;
                    }
                    if (gas_state(coordinate) != GasState::Solid) {
                        if (removed_voxels != nullptr) {
                            removed_voxels->push_back({
                                coordinate,
                                gas_state(coordinate)});
                        }
                        set_owned_gas_state(coordinate, GasState::Solid);
                        ++newly_solid_count;
                    }
                }
            }
        }
    }
    return newly_solid_count;
}

void DistributedGasGrid::exchange_ghost_states()
{
    std::array<MPI_Request, 12> requests{};
    int request_count = 0;

    for (std::size_t index = 0; index < send_buffers_.size(); ++index) {
        const auto face = static_cast<Face>(index);
        if (decomposition_.neighbor_rank(face) != MPI_PROC_NULL) {
            pack_face(face);
        }
    }

    for (std::size_t index = 0; index < send_buffers_.size(); ++index) {
        const auto face = static_cast<Face>(index);
        const int neighbor = decomposition_.neighbor_rank(face);
        if (neighbor == MPI_PROC_NULL) {
            fill_ghost_face(face, GasState::Unclassified);
            continue;
        }
        if (neighbor == decomposition_.rank()) {
            const auto& source = send_buffers_[face_index(opposite_face(face))];
            std::copy(
                source.begin(),
                source.end(),
                receive_buffers_[index].begin());
            continue;
        }

        const int count = static_cast<int>(receive_buffers_[index].size());
        check_mpi(MPI_Irecv(
            receive_buffers_[index].data(),
            count,
            MPI_BYTE,
            neighbor,
            halo_tag_base + static_cast<int>(face_index(opposite_face(face))),
            decomposition_.spec().communicator,
            &requests[static_cast<std::size_t>(request_count)]),
            "MPI_Irecv gas halo");
        ++request_count;
    }

    for (std::size_t index = 0; index < send_buffers_.size(); ++index) {
        const auto face = static_cast<Face>(index);
        const int neighbor = decomposition_.neighbor_rank(face);
        if (neighbor == MPI_PROC_NULL || neighbor == decomposition_.rank()) {
            continue;
        }
        const int count = static_cast<int>(send_buffers_[index].size());
        check_mpi(MPI_Isend(
            send_buffers_[index].data(),
            count,
            MPI_BYTE,
            neighbor,
            halo_tag_base + static_cast<int>(index),
            decomposition_.spec().communicator,
            &requests[static_cast<std::size_t>(request_count)]),
            "MPI_Isend gas halo");
        ++request_count;
    }

    if (request_count != 0) {
        check_mpi(MPI_Waitall(request_count, requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall gas halo");
    }
    for (std::size_t index = 0; index < receive_buffers_.size(); ++index) {
        const auto face = static_cast<Face>(index);
        if (decomposition_.neighbor_rank(face) != MPI_PROC_NULL) {
            unpack_face(face);
        }
    }
}

bool DistributedGasGrid::is_site_accessible(const Point3& site_position) const
{
    const auto containing_voxel = locate_voxel(site_position);
    if (!containing_voxel) {
        return false;
    }
    if (!owns(*containing_voxel)) {
        throw std::out_of_range(
            "site position is not owned by this MPI rank's gas grid");
    }
    if (gas_state(*containing_voxel) == GasState::OutsideAccessible) {
        return true;
    }

    const std::array<VoxelCoord, 6> candidates{{
        {containing_voxel->x - 1, containing_voxel->y, containing_voxel->z},
        {containing_voxel->x + 1, containing_voxel->y, containing_voxel->z},
        {containing_voxel->x, containing_voxel->y - 1, containing_voxel->z},
        {containing_voxel->x, containing_voxel->y + 1, containing_voxel->z},
        {containing_voxel->x, containing_voxel->y, containing_voxel->z - 1},
        {containing_voxel->x, containing_voxel->y, containing_voxel->z + 1}
    }};
    for (const auto& candidate : candidates) {
        const auto normalized = normalized_neighbor(candidate);
        if (normalized && *normalized != *containing_voxel
            && gas_state(*normalized) == GasState::OutsideAccessible) {
            return true;
        }
    }
    return false;
}

std::size_t DistributedGasGrid::face_index(Face face) noexcept
{
    return static_cast<std::size_t>(face);
}

Face DistributedGasGrid::opposite_face(Face face) noexcept
{
    switch (face) {
    case Face::XLow:
        return Face::XHigh;
    case Face::XHigh:
        return Face::XLow;
    case Face::YLow:
        return Face::YHigh;
    case Face::YHigh:
        return Face::YLow;
    case Face::ZLow:
        return Face::ZHigh;
    case Face::ZHigh:
        return Face::ZLow;
    }
    return Face::XLow;
}

std::size_t DistributedGasGrid::gas_state_index(GasState gas_state) noexcept
{
    return static_cast<std::size_t>(gas_state);
}

std::optional<DistributedGasGrid::LocalCoord> DistributedGasGrid::local_coord(
    const VoxelCoord& global_voxel_coord) const noexcept
{
    if (global_voxel_coord.x < 0
        || global_voxel_coord.y < 0
        || global_voxel_coord.z < 0
        || static_cast<std::uint64_t>(global_voxel_coord.x)
            >= global_grid_spec_.dimensions.x
        || static_cast<std::uint64_t>(global_voxel_coord.y)
            >= global_grid_spec_.dimensions.y
        || static_cast<std::uint64_t>(global_voxel_coord.z)
            >= global_grid_spec_.dimensions.z) {
        return std::nullopt;
    }

    const auto& range = owned_range();
    if (owns(global_voxel_coord)) {
        return LocalCoord{
            static_cast<std::uint64_t>(
                global_voxel_coord.x - range.begin.x + 1),
            static_cast<std::uint64_t>(
                global_voxel_coord.y - range.begin.y + 1),
            static_cast<std::uint64_t>(
                global_voxel_coord.z - range.begin.z + 1)
        };
    }

    const bool x_owned = global_voxel_coord.x >= range.begin.x
        && global_voxel_coord.x < range.end.x;
    const bool y_owned = global_voxel_coord.y >= range.begin.y
        && global_voxel_coord.y < range.end.y;
    const bool z_owned = global_voxel_coord.z >= range.begin.z
        && global_voxel_coord.z < range.end.z;

    const auto low_neighbor_coordinate = [](std::int64_t begin,
                                            std::uint64_t dimension,
                                            bool periodic)
        -> std::optional<std::int64_t> {
        if (begin > 0) {
            return begin - 1;
        }
        return periodic
            ? std::optional<std::int64_t>(
                static_cast<std::int64_t>(dimension - 1))
            : std::nullopt;
    };
    const auto high_neighbor_coordinate = [](std::int64_t end,
                                             std::uint64_t dimension,
                                             bool periodic)
        -> std::optional<std::int64_t> {
        if (static_cast<std::uint64_t>(end) < dimension) {
            return end;
        }
        return periodic ? std::optional<std::int64_t>(0) : std::nullopt;
    };

    if (y_owned && z_owned) {
        const auto low = low_neighbor_coordinate(
            range.begin.x,
            global_grid_spec_.dimensions.x,
            global_grid_spec_.periodic.x);
        const auto high = high_neighbor_coordinate(
            range.end.x,
            global_grid_spec_.dimensions.x,
            global_grid_spec_.periodic.x);
        if (low && global_voxel_coord.x == *low) {
            return LocalCoord{
                0,
                static_cast<std::uint64_t>(
                    global_voxel_coord.y - range.begin.y + 1),
                static_cast<std::uint64_t>(
                    global_voxel_coord.z - range.begin.z + 1)
            };
        }
        if (high && global_voxel_coord.x == *high) {
            return LocalCoord{
                range.dimensions.x + 1,
                static_cast<std::uint64_t>(
                    global_voxel_coord.y - range.begin.y + 1),
                static_cast<std::uint64_t>(
                    global_voxel_coord.z - range.begin.z + 1)
            };
        }
    }
    if (x_owned && z_owned) {
        const auto low = low_neighbor_coordinate(
            range.begin.y,
            global_grid_spec_.dimensions.y,
            global_grid_spec_.periodic.y);
        const auto high = high_neighbor_coordinate(
            range.end.y,
            global_grid_spec_.dimensions.y,
            global_grid_spec_.periodic.y);
        if (low && global_voxel_coord.y == *low) {
            return LocalCoord{
                static_cast<std::uint64_t>(
                    global_voxel_coord.x - range.begin.x + 1),
                0,
                static_cast<std::uint64_t>(
                    global_voxel_coord.z - range.begin.z + 1)
            };
        }
        if (high && global_voxel_coord.y == *high) {
            return LocalCoord{
                static_cast<std::uint64_t>(
                    global_voxel_coord.x - range.begin.x + 1),
                range.dimensions.y + 1,
                static_cast<std::uint64_t>(
                    global_voxel_coord.z - range.begin.z + 1)
            };
        }
    }
    if (x_owned && y_owned) {
        const auto low = low_neighbor_coordinate(
            range.begin.z,
            global_grid_spec_.dimensions.z,
            global_grid_spec_.periodic.z);
        const auto high = high_neighbor_coordinate(
            range.end.z,
            global_grid_spec_.dimensions.z,
            global_grid_spec_.periodic.z);
        if (low && global_voxel_coord.z == *low) {
            return LocalCoord{
                static_cast<std::uint64_t>(
                    global_voxel_coord.x - range.begin.x + 1),
                static_cast<std::uint64_t>(
                    global_voxel_coord.y - range.begin.y + 1),
                0
            };
        }
        if (high && global_voxel_coord.z == *high) {
            return LocalCoord{
                static_cast<std::uint64_t>(
                    global_voxel_coord.x - range.begin.x + 1),
                static_cast<std::uint64_t>(
                    global_voxel_coord.y - range.begin.y + 1),
                range.dimensions.z + 1
            };
        }
    }
    return std::nullopt;
}

std::optional<VoxelCoord> DistributedGasGrid::normalized_neighbor(
    const VoxelCoord& global_voxel_coord) const noexcept
{
    VoxelCoord normalized = global_voxel_coord;
    const auto normalize_axis = [](std::int64_t& coordinate,
                                   std::uint64_t dimension,
                                   bool periodic) {
        if (coordinate < 0) {
            if (!periodic) {
                return false;
            }
            coordinate = static_cast<std::int64_t>(dimension - 1);
        } else if (static_cast<std::uint64_t>(coordinate) >= dimension) {
            if (!periodic) {
                return false;
            }
            coordinate = 0;
        }
        return true;
    };
    if (!normalize_axis(
            normalized.x,
            global_grid_spec_.dimensions.x,
            global_grid_spec_.periodic.x)
        || !normalize_axis(
            normalized.y,
            global_grid_spec_.dimensions.y,
            global_grid_spec_.periodic.y)
        || !normalize_axis(
            normalized.z,
            global_grid_spec_.dimensions.z,
            global_grid_spec_.periodic.z)) {
        return std::nullopt;
    }
    return normalized;
}

std::size_t DistributedGasGrid::state_index(
    const LocalCoord& local_coordinate) const noexcept
{
    const auto index = (local_coordinate.z * storage_dimensions_.y
        + local_coordinate.y) * storage_dimensions_.x + local_coordinate.x;
    return static_cast<std::size_t>(index);
}

std::size_t DistributedGasGrid::face_element_count(Face face) const
{
    const auto& dimensions = owned_range().dimensions;
    std::uint64_t count = 0;
    switch (face) {
    case Face::XLow:
    case Face::XHigh:
        count = checked_multiply(dimensions.y, dimensions.z);
        break;
    case Face::YLow:
    case Face::YHigh:
        count = checked_multiply(dimensions.x, dimensions.z);
        break;
    case Face::ZLow:
    case Face::ZHigh:
        count = checked_multiply(dimensions.x, dimensions.y);
        break;
    }
    if (count > static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max())) {
        throw std::length_error("gas-state face exceeds size_t range");
    }
    return static_cast<std::size_t>(count);
}

void DistributedGasGrid::pack_face(Face face)
{
    auto& buffer = send_buffers_[face_index(face)];
    const auto& dimensions = owned_range().dimensions;
    std::size_t offset = 0;
    if (face == Face::XLow || face == Face::XHigh) {
        const auto x = face == Face::XLow ? std::uint64_t{1} : dimensions.x;
        for (std::uint64_t z = 1; z <= dimensions.z; ++z) {
            for (std::uint64_t y = 1; y <= dimensions.y; ++y) {
                buffer[offset++] = states_[state_index({x, y, z})];
            }
        }
    } else if (face == Face::YLow || face == Face::YHigh) {
        const auto y = face == Face::YLow ? std::uint64_t{1} : dimensions.y;
        for (std::uint64_t z = 1; z <= dimensions.z; ++z) {
            for (std::uint64_t x = 1; x <= dimensions.x; ++x) {
                buffer[offset++] = states_[state_index({x, y, z})];
            }
        }
    } else {
        const auto z = face == Face::ZLow ? std::uint64_t{1} : dimensions.z;
        for (std::uint64_t y = 1; y <= dimensions.y; ++y) {
            for (std::uint64_t x = 1; x <= dimensions.x; ++x) {
                buffer[offset++] = states_[state_index({x, y, z})];
            }
        }
    }
}

void DistributedGasGrid::unpack_face(Face face)
{
    const auto& buffer = receive_buffers_[face_index(face)];
    const auto& dimensions = owned_range().dimensions;
    std::size_t offset = 0;
    if (face == Face::XLow || face == Face::XHigh) {
        const auto x = face == Face::XLow ? std::uint64_t{0} : dimensions.x + 1;
        for (std::uint64_t z = 1; z <= dimensions.z; ++z) {
            for (std::uint64_t y = 1; y <= dimensions.y; ++y) {
                states_[state_index({x, y, z})] = buffer[offset++];
            }
        }
    } else if (face == Face::YLow || face == Face::YHigh) {
        const auto y = face == Face::YLow ? std::uint64_t{0} : dimensions.y + 1;
        for (std::uint64_t z = 1; z <= dimensions.z; ++z) {
            for (std::uint64_t x = 1; x <= dimensions.x; ++x) {
                states_[state_index({x, y, z})] = buffer[offset++];
            }
        }
    } else {
        const auto z = face == Face::ZLow ? std::uint64_t{0} : dimensions.z + 1;
        for (std::uint64_t y = 1; y <= dimensions.y; ++y) {
            for (std::uint64_t x = 1; x <= dimensions.x; ++x) {
                states_[state_index({x, y, z})] = buffer[offset++];
            }
        }
    }
}

void DistributedGasGrid::fill_ghost_face(Face face, GasState gas_state_value)
{
    std::fill(
        receive_buffers_[face_index(face)].begin(),
        receive_buffers_[face_index(face)].end(),
        gas_state_value);
    unpack_face(face);
}

DistributedGasAccessibilityQuery::DistributedGasAccessibilityQuery(
    const DistributedGasGrid& gas_grid) noexcept
    : gas_grid_(gas_grid)
{
}

bool DistributedGasAccessibilityQuery::is_site_accessible(
    const Point3& site_position) const
{
    return gas_grid_.is_site_accessible(site_position);
}

}  // namespace gasaccess
