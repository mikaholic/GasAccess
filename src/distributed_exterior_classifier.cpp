#include "gasaccess/distributed_exterior_classifier.hpp"

#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace gasaccess {
namespace {

constexpr std::array<Face, 6> faces{
    Face::XLow,
    Face::XHigh,
    Face::YLow,
    Face::YHigh,
    Face::ZLow,
    Face::ZHigh};

constexpr std::array<VoxelCoord, 6> neighbor_offsets{
    VoxelCoord{-1, 0, 0},
    VoxelCoord{1, 0, 0},
    VoxelCoord{0, -1, 0},
    VoxelCoord{0, 1, 0},
    VoxelCoord{0, 0, -1},
    VoxelCoord{0, 0, 1}};

constexpr int count_tag_base = 5410;
constexpr int payload_tag_base = 5420;

std::size_t face_index(Face face) noexcept
{
    return static_cast<std::size_t>(face);
}

Face opposite_face(Face face) noexcept
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

void check_mpi(int error_code, const char* operation)
{
    if (error_code == MPI_SUCCESS) {
        return;
    }

    std::array<char, MPI_MAX_ERROR_STRING> error_message{};
    int error_length = 0;
    MPI_Error_string(error_code, error_message.data(), &error_length);
    throw std::runtime_error(
        std::string(operation) + " failed: "
        + std::string(error_message.data(), static_cast<std::size_t>(error_length)));
}

std::int64_t dimension_as_int64(std::uint64_t dimension)
{
    if (dimension > static_cast<std::uint64_t>(
                        std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error("grid dimension exceeds signed coordinate range");
    }
    return static_cast<std::int64_t>(dimension);
}

std::optional<VoxelCoord> normalized_neighbor(
    const VoxelCoord& voxel_coord,
    const GridSpec& grid_spec)
{
    VoxelCoord normalized = voxel_coord;

    const auto normalize_axis = [](
                                    std::int64_t coordinate,
                                    std::uint64_t dimension,
                                    bool periodic)
        -> std::optional<std::int64_t> {
        const auto signed_dimension = dimension_as_int64(dimension);
        if (coordinate >= 0 && coordinate < signed_dimension) {
            return coordinate;
        }
        if (!periodic) {
            return std::nullopt;
        }
        if (coordinate < 0) {
            return signed_dimension - 1;
        }
        return 0;
    };

    const auto x = normalize_axis(
        normalized.x,
        grid_spec.dimensions.x,
        grid_spec.periodic.x);
    const auto y = normalize_axis(
        normalized.y,
        grid_spec.dimensions.y,
        grid_spec.periodic.y);
    const auto z = normalize_axis(
        normalized.z,
        grid_spec.dimensions.z,
        grid_spec.periodic.z);
    if (!x.has_value() || !y.has_value() || !z.has_value()) {
        return std::nullopt;
    }

    normalized.x = *x;
    normalized.y = *y;
    normalized.z = *z;
    return normalized;
}

bool is_boundary_source(
    const VoxelCoord& voxel_coord,
    const GridSpec& grid_spec)
{
    const auto x_high = dimension_as_int64(grid_spec.dimensions.x) - 1;
    const auto y_high = dimension_as_int64(grid_spec.dimensions.y) - 1;
    const auto z_high = dimension_as_int64(grid_spec.dimensions.z) - 1;
    const auto& sources = grid_spec.reservoir_faces;

    return (sources.x_low && voxel_coord.x == 0)
        || (sources.x_high && voxel_coord.x == x_high)
        || (sources.y_low && voxel_coord.y == 0)
        || (sources.y_high && voxel_coord.y == y_high)
        || (sources.z_low && voxel_coord.z == 0)
        || (sources.z_high && voxel_coord.z == z_high);
}

std::uint64_t checked_product(std::uint64_t lhs, std::uint64_t rhs)
{
    if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
        throw std::overflow_error("frontier face size overflows uint64_t");
    }
    return lhs * rhs;
}

std::uint64_t face_element_count(
    Face face,
    const OwnedVoxelRange& owned_range)
{
    switch (face) {
    case Face::XLow:
    case Face::XHigh:
        return checked_product(
            owned_range.dimensions.y,
            owned_range.dimensions.z);
    case Face::YLow:
    case Face::YHigh:
        return checked_product(
            owned_range.dimensions.x,
            owned_range.dimensions.z);
    case Face::ZLow:
    case Face::ZHigh:
        return checked_product(
            owned_range.dimensions.x,
            owned_range.dimensions.y);
    }

    return 0;
}

std::uint64_t coordinate_difference(std::int64_t value, std::int64_t begin)
{
    if (value < begin) {
        throw std::logic_error("frontier coordinate is outside the owned face");
    }
    return static_cast<std::uint64_t>(value - begin);
}

std::uint64_t encode_face_offset(
    Face face,
    const VoxelCoord& source_coord,
    const OwnedVoxelRange& owned_range)
{
    const auto local_x = coordinate_difference(
        source_coord.x,
        owned_range.begin.x);
    const auto local_y = coordinate_difference(
        source_coord.y,
        owned_range.begin.y);
    const auto local_z = coordinate_difference(
        source_coord.z,
        owned_range.begin.z);

    switch (face) {
    case Face::XLow:
    case Face::XHigh:
        return local_y + owned_range.dimensions.y * local_z;
    case Face::YLow:
    case Face::YHigh:
        return local_x + owned_range.dimensions.x * local_z;
    case Face::ZLow:
    case Face::ZHigh:
        return local_x + owned_range.dimensions.x * local_y;
    }

    return 0;
}

VoxelCoord decode_face_offset(
    Face face,
    std::uint64_t offset,
    const OwnedVoxelRange& owned_range)
{
    if (offset >= face_element_count(face, owned_range)) {
        throw std::runtime_error("received frontier offset exceeds owned face");
    }

    VoxelCoord voxel_coord{};
    switch (face) {
    case Face::XLow:
    case Face::XHigh:
        voxel_coord.x = face == Face::XLow
            ? owned_range.begin.x
            : owned_range.end.x - 1;
        voxel_coord.y = owned_range.begin.y
            + static_cast<std::int64_t>(offset % owned_range.dimensions.y);
        voxel_coord.z = owned_range.begin.z
            + static_cast<std::int64_t>(offset / owned_range.dimensions.y);
        break;
    case Face::YLow:
    case Face::YHigh:
        voxel_coord.x = owned_range.begin.x
            + static_cast<std::int64_t>(offset % owned_range.dimensions.x);
        voxel_coord.y = face == Face::YLow
            ? owned_range.begin.y
            : owned_range.end.y - 1;
        voxel_coord.z = owned_range.begin.z
            + static_cast<std::int64_t>(offset / owned_range.dimensions.x);
        break;
    case Face::ZLow:
    case Face::ZHigh:
        voxel_coord.x = owned_range.begin.x
            + static_cast<std::int64_t>(offset % owned_range.dimensions.x);
        voxel_coord.y = owned_range.begin.y
            + static_cast<std::int64_t>(offset / owned_range.dimensions.x);
        voxel_coord.z = face == Face::ZLow
            ? owned_range.begin.z
            : owned_range.end.z - 1;
        break;
    }

    return voxel_coord;
}

int checked_mpi_count(std::uint64_t count)
{
    if (count > static_cast<std::uint64_t>(INT_MAX)) {
        throw std::overflow_error("MPI frontier payload exceeds INT_MAX entries");
    }
    return static_cast<int>(count);
}

void mark_outside(
    DistributedGasGrid& gas_grid,
    const VoxelCoord& voxel_coord,
    std::vector<VoxelCoord>& frontier,
    DistributedClassificationSummary& summary)
{
    if (gas_grid.gas_state(voxel_coord) != GasState::ClosedVoid) {
        return;
    }

    gas_grid.set_owned_gas_state(voxel_coord, GasState::OutsideAccessible);
    frontier.push_back(voxel_coord);
    ++summary.local_outside_accessible_count;
    --summary.local_closed_void_count;
    ++summary.local_visited_voxel_count;
}

void exchange_frontiers(
    const MpiDecomposition& decomposition,
    std::array<std::vector<std::uint64_t>, 6>& send_buffers,
    std::array<std::vector<std::uint64_t>, 6>& receive_buffers,
    DistributedClassificationSummary& summary)
{
    std::array<std::uint64_t, 6> send_counts{};
    std::array<std::uint64_t, 6> receive_counts{};
    std::array<MPI_Request, 12> requests{};
    int request_count = 0;
    const auto communicator = decomposition.spec().communicator;

    for (const auto face : faces) {
        const auto index = face_index(face);
        send_counts[index] = static_cast<std::uint64_t>(send_buffers[index].size());
        receive_buffers[index].clear();

        const auto neighbor_rank = decomposition.neighbor_rank(face);
        if (neighbor_rank == MPI_PROC_NULL) {
            continue;
        }
        if (neighbor_rank == decomposition.rank()) {
            const auto opposite_index = face_index(opposite_face(face));
            receive_counts[index] = static_cast<std::uint64_t>(
                send_buffers[opposite_index].size());
            continue;
        }

        check_mpi(
            MPI_Irecv(
                &receive_counts[index],
                1,
                MPI_UINT64_T,
                neighbor_rank,
                count_tag_base + static_cast<int>(opposite_face(face)),
                communicator,
                &requests[static_cast<std::size_t>(request_count)]),
            "MPI_Irecv(frontier count)");
        ++request_count;
        check_mpi(
            MPI_Isend(
                &send_counts[index],
                1,
                MPI_UINT64_T,
                neighbor_rank,
                count_tag_base + static_cast<int>(face),
                communicator,
                &requests[static_cast<std::size_t>(request_count)]),
            "MPI_Isend(frontier count)");
        ++request_count;
    }

    if (request_count != 0) {
        check_mpi(
            MPI_Waitall(request_count, requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall(frontier counts)");
    }

    request_count = 0;
    for (const auto face : faces) {
        const auto index = face_index(face);
        const auto neighbor_rank = decomposition.neighbor_rank(face);
        if (neighbor_rank == MPI_PROC_NULL) {
            continue;
        }
        if (neighbor_rank == decomposition.rank()) {
            const auto opposite_index = face_index(opposite_face(face));
            receive_buffers[index].assign(
                send_buffers[opposite_index].begin(),
                send_buffers[opposite_index].end());
            continue;
        }

        receive_buffers[index].resize(
            static_cast<std::size_t>(receive_counts[index]));
        if (receive_counts[index] != 0) {
            check_mpi(
                MPI_Irecv(
                    receive_buffers[index].data(),
                    checked_mpi_count(receive_counts[index]),
                    MPI_UINT64_T,
                    neighbor_rank,
                    payload_tag_base + static_cast<int>(opposite_face(face)),
                    communicator,
                    &requests[static_cast<std::size_t>(request_count)]),
                "MPI_Irecv(frontier payload)");
            ++request_count;
        }
        if (send_counts[index] != 0) {
            check_mpi(
                MPI_Isend(
                    send_buffers[index].data(),
                    checked_mpi_count(send_counts[index]),
                    MPI_UINT64_T,
                    neighbor_rank,
                    payload_tag_base + static_cast<int>(face),
                    communicator,
                    &requests[static_cast<std::size_t>(request_count)]),
                "MPI_Isend(frontier payload)");
            ++request_count;
        }
    }

    if (request_count != 0) {
        check_mpi(
            MPI_Waitall(request_count, requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall(frontier payloads)");
    }

    for (const auto face : faces) {
        const auto index = face_index(face);
        summary.sent_frontier_entry_count += send_counts[index];
        summary.received_frontier_entry_count += receive_counts[index];
    }
}

}  // namespace

DistributedClassificationSummary DistributedExteriorClassifier::classify(
    DistributedGasGrid& gas_grid)
{
    DistributedClassificationSummary summary{};
    frontier_.clear();
    frontier_head_ = 0;
    for (auto& send_buffer : send_buffers_) {
        send_buffer.clear();
    }
    for (auto& receive_buffer : receive_buffers_) {
        receive_buffer.clear();
    }

    const auto& grid_spec = gas_grid.global_grid_spec();
    const auto& owned_range = gas_grid.owned_range();
    for (auto z = owned_range.begin.z; z < owned_range.end.z; ++z) {
        for (auto y = owned_range.begin.y; y < owned_range.end.y; ++y) {
            for (auto x = owned_range.begin.x; x < owned_range.end.x; ++x) {
                const VoxelCoord voxel_coord{x, y, z};
                if (gas_grid.gas_state(voxel_coord) == GasState::Solid) {
                    ++summary.local_solid_count;
                    continue;
                }

                gas_grid.set_owned_gas_state(voxel_coord, GasState::ClosedVoid);
                ++summary.local_closed_void_count;
                if (is_boundary_source(voxel_coord, grid_spec)) {
                    mark_outside(gas_grid, voxel_coord, frontier_, summary);
                }
            }
        }
    }

    for (const auto& source_coord : grid_spec.explicit_source_voxels) {
        if (gas_grid.owns(source_coord)) {
            mark_outside(gas_grid, source_coord, frontier_, summary);
        }
    }

    int local_active = frontier_.empty() ? 0 : 1;
    int global_active = 0;
    check_mpi(
        MPI_Allreduce(
            &local_active,
            &global_active,
            1,
            MPI_INT,
            MPI_MAX,
            gas_grid.decomposition().spec().communicator),
        "MPI_Allreduce(initial frontier)");

    while (global_active != 0) {
        for (auto& send_buffer : send_buffers_) {
            send_buffer.clear();
        }

        while (frontier_head_ < frontier_.size()) {
            const auto current_coord = frontier_[frontier_head_];
            ++frontier_head_;

            for (std::size_t neighbor_index = 0;
                 neighbor_index < neighbor_offsets.size();
                 ++neighbor_index) {
                const auto offset = neighbor_offsets[neighbor_index];
                const VoxelCoord candidate{
                    current_coord.x + offset.x,
                    current_coord.y + offset.y,
                    current_coord.z + offset.z};
                const auto neighbor_coord = normalized_neighbor(candidate, grid_spec);
                if (!neighbor_coord.has_value() || *neighbor_coord == current_coord) {
                    continue;
                }

                if (gas_grid.owns(*neighbor_coord)) {
                    mark_outside(gas_grid, *neighbor_coord, frontier_, summary);
                    continue;
                }

                const auto face = faces[neighbor_index];
                if (gas_grid.decomposition().neighbor_rank(face) != MPI_PROC_NULL) {
                    send_buffers_[face_index(face)].push_back(
                        encode_face_offset(face, current_coord, owned_range));
                }
            }
        }

        frontier_.clear();
        frontier_head_ = 0;
        exchange_frontiers(
            gas_grid.decomposition(),
            send_buffers_,
            receive_buffers_,
            summary);
        ++summary.communication_round_count;

        for (const auto face : faces) {
            for (const auto face_offset : receive_buffers_[face_index(face)]) {
                const auto voxel_coord = decode_face_offset(
                    face,
                    face_offset,
                    owned_range);
                mark_outside(gas_grid, voxel_coord, frontier_, summary);
            }
        }

        local_active = frontier_.empty() ? 0 : 1;
        check_mpi(
            MPI_Allreduce(
                &local_active,
                &global_active,
                1,
                MPI_INT,
                MPI_MAX,
                gas_grid.decomposition().spec().communicator),
            "MPI_Allreduce(frontier termination)");
    }

    gas_grid.exchange_ghost_states();
    return summary;
}

}  // namespace gasaccess
