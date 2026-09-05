#include "gasaccess/distributed_closing_region_repair.hpp"

#include <algorithm>
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

constexpr int count_tag_base = 5710;
constexpr int payload_tag_base = 5720;
constexpr std::size_t encoded_change_width = 4;

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
    static_cast<void>(
        MPI_Error_string(error_code, error_message.data(), &error_length));
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
        return coordinate < 0 ? signed_dimension - 1 : 0;
    };

    const auto x = normalize_axis(
        voxel_coord.x,
        grid_spec.dimensions.x,
        grid_spec.periodic.x);
    const auto y = normalize_axis(
        voxel_coord.y,
        grid_spec.dimensions.y,
        grid_spec.periodic.y);
    const auto z = normalize_axis(
        voxel_coord.z,
        grid_spec.dimensions.z,
        grid_spec.periodic.z);
    if (!x.has_value() || !y.has_value() || !z.has_value()) {
        return std::nullopt;
    }
    return VoxelCoord{*x, *y, *z};
}

bool is_reservoir_source(
    const VoxelCoord& voxel_coord,
    const GridSpec& grid_spec)
{
    const auto x_high = dimension_as_int64(grid_spec.dimensions.x) - 1;
    const auto y_high = dimension_as_int64(grid_spec.dimensions.y) - 1;
    const auto z_high = dimension_as_int64(grid_spec.dimensions.z) - 1;
    const auto& sources = grid_spec.reservoir_faces;
    if ((sources.x_low && voxel_coord.x == 0)
        || (sources.x_high && voxel_coord.x == x_high)
        || (sources.y_low && voxel_coord.y == 0)
        || (sources.y_high && voxel_coord.y == y_high)
        || (sources.z_low && voxel_coord.z == 0)
        || (sources.z_high && voxel_coord.z == z_high)) {
        return true;
    }
    return std::find(
        grid_spec.explicit_source_voxels.begin(),
        grid_spec.explicit_source_voxels.end(),
        voxel_coord) != grid_spec.explicit_source_voxels.end();
}

bool voxel_coord_less(const VoxelCoord& lhs, const VoxelCoord& rhs) noexcept
{
    if (lhs.z != rhs.z) {
        return lhs.z < rhs.z;
    }
    if (lhs.y != rhs.y) {
        return lhs.y < rhs.y;
    }
    return lhs.x < rhs.x;
}

void sort_unique(std::vector<VoxelCoord>& voxel_coords)
{
    std::sort(voxel_coords.begin(), voxel_coords.end(), voxel_coord_less);
    voxel_coords.erase(
        std::unique(voxel_coords.begin(), voxel_coords.end()),
        voxel_coords.end());
}

std::size_t owned_index(
    const VoxelCoord& voxel_coord,
    const OwnedVoxelRange& owned_range)
{
    if (voxel_coord.x < owned_range.begin.x
        || voxel_coord.x >= owned_range.end.x
        || voxel_coord.y < owned_range.begin.y
        || voxel_coord.y >= owned_range.end.y
        || voxel_coord.z < owned_range.begin.z
        || voxel_coord.z >= owned_range.end.z) {
        throw std::out_of_range("voxel coordinate is outside the owned range");
    }
    const auto x = static_cast<std::uint64_t>(
        voxel_coord.x - owned_range.begin.x);
    const auto y = static_cast<std::uint64_t>(
        voxel_coord.y - owned_range.begin.y);
    const auto z = static_cast<std::uint64_t>(
        voxel_coord.z - owned_range.begin.z);
    const auto index = x + owned_range.dimensions.x
        * (y + owned_range.dimensions.y * z);
    if (index > static_cast<std::uint64_t>(
                    std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("owned voxel index exceeds size_t");
    }
    return static_cast<std::size_t>(index);
}

std::uint64_t checked_product(std::uint64_t lhs, std::uint64_t rhs)
{
    if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
        throw std::overflow_error("closing-repair face size overflow");
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

std::uint64_t encode_face_offset(
    Face face,
    const VoxelCoord& source_coord,
    const OwnedVoxelRange& owned_range)
{
    const auto x = static_cast<std::uint64_t>(
        source_coord.x - owned_range.begin.x);
    const auto y = static_cast<std::uint64_t>(
        source_coord.y - owned_range.begin.y);
    const auto z = static_cast<std::uint64_t>(
        source_coord.z - owned_range.begin.z);
    switch (face) {
    case Face::XLow:
    case Face::XHigh:
        return y + owned_range.dimensions.y * z;
    case Face::YLow:
    case Face::YHigh:
        return x + owned_range.dimensions.x * z;
    case Face::ZLow:
    case Face::ZHigh:
        return x + owned_range.dimensions.x * y;
    }
    return 0;
}

VoxelCoord decode_face_offset(
    Face face,
    std::uint64_t offset,
    const OwnedVoxelRange& owned_range)
{
    if (offset >= face_element_count(face, owned_range)) {
        throw std::runtime_error(
            "received closing-repair offset exceeds owned face");
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
        throw std::overflow_error(
            "MPI closing-repair payload exceeds INT_MAX entries");
    }
    return static_cast<int>(count);
}

void exchange_frontiers(
    const MpiDecomposition& decomposition,
    std::array<std::vector<std::uint64_t>, 6>& send_buffers,
    std::array<std::vector<std::uint64_t>, 6>& receive_buffers,
    DistributedClosingRegionRepairResult& result)
{
    std::array<std::uint64_t, 6> send_counts{};
    std::array<std::uint64_t, 6> receive_counts{};
    std::array<MPI_Request, 12> requests{};
    int request_count = 0;
    const auto communicator = decomposition.spec().communicator;

    for (const auto face : faces) {
        const auto index = face_index(face);
        send_counts[index] = static_cast<std::uint64_t>(
            send_buffers[index].size());
        receive_buffers[index].clear();
        const auto neighbor_rank = decomposition.neighbor_rank(face);
        if (neighbor_rank == MPI_PROC_NULL) {
            continue;
        }
        if (neighbor_rank == decomposition.rank()) {
            receive_counts[index] = static_cast<std::uint64_t>(
                send_buffers[face_index(opposite_face(face))].size());
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
            "MPI_Irecv(closing-repair count)");
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
            "MPI_Isend(closing-repair count)");
        ++request_count;
    }
    if (request_count != 0) {
        check_mpi(
            MPI_Waitall(request_count, requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall(closing-repair counts)");
    }

    request_count = 0;
    for (const auto face : faces) {
        const auto index = face_index(face);
        const auto neighbor_rank = decomposition.neighbor_rank(face);
        if (neighbor_rank == MPI_PROC_NULL) {
            continue;
        }
        if (neighbor_rank == decomposition.rank()) {
            const auto& source = send_buffers[face_index(opposite_face(face))];
            receive_buffers[index].assign(source.begin(), source.end());
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
                "MPI_Irecv(closing-repair payload)");
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
                "MPI_Isend(closing-repair payload)");
            ++request_count;
        }
    }
    if (request_count != 0) {
        check_mpi(
            MPI_Waitall(request_count, requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall(closing-repair payloads)");
    }

    for (const auto face : faces) {
        const auto index = face_index(face);
        result.sent_frontier_entry_count += send_counts[index];
        result.received_frontier_entry_count += receive_counts[index];
    }
}

void reduce_search_state(
    MPI_Comm communicator,
    bool local_active,
    bool local_source_found,
    bool& global_active,
    bool& global_source_found)
{
    const std::array<int, 2> local_values{
        local_active ? 1 : 0,
        local_source_found ? 1 : 0};
    std::array<int, 2> global_values{};
    check_mpi(
        MPI_Allreduce(
            local_values.data(),
            global_values.data(),
            static_cast<int>(global_values.size()),
            MPI_INT,
            MPI_MAX,
            communicator),
        "MPI_Allreduce(closing-repair search state)");
    global_active = global_values[0] != 0;
    global_source_found = global_values[1] != 0;
}

void validate_collectively(
    const DistributedGasGrid& gas_grid,
    NewlySolidVoxelCoordView newly_solid_voxel_view)
{
    bool local_valid = newly_solid_voxel_view.count == 0
        || newly_solid_voxel_view.voxels != nullptr;
    if (local_valid) {
        for (std::size_t index = 0;
             index < newly_solid_voxel_view.count;
             ++index) {
            const auto& voxel = newly_solid_voxel_view.voxels[index];
            if (!gas_grid.owns(voxel.voxel_coord)
                || gas_grid.gas_state(voxel.voxel_coord) != GasState::Solid
                || (voxel.previous_state != GasState::OutsideAccessible
                    && voxel.previous_state != GasState::ClosedVoid)) {
                local_valid = false;
                break;
            }
        }
    }

    const int local_value = local_valid ? 1 : 0;
    int global_value = 0;
    check_mpi(
        MPI_Allreduce(
            &local_value,
            &global_value,
            1,
            MPI_INT,
            MPI_MIN,
            gas_grid.decomposition().spec().communicator),
        "MPI_Allreduce(closing-repair validation)");
    if (global_value == 0) {
        throw std::invalid_argument(
            "distributed closing repair received invalid rank-local input");
    }
}

}  // namespace

DistributedClosingRegionRepairResult DistributedClosingRegionRepair::repair(
    DistributedGasGrid& gas_grid,
    NewlySolidVoxelCoordView newly_solid_voxel_view)
{
    validate_collectively(gas_grid, newly_solid_voxel_view);
    gather_newly_solid_voxels(gas_grid, newly_solid_voxel_view);
    prepare_workspace(gas_grid);
    begin_repair_epoch();

    DistributedClosingRegionRepairResult result{};
    seed_voxels_.clear();
    const auto& grid_spec = gas_grid.global_grid_spec();
    const auto& owned_range = gas_grid.owned_range();
    const auto communicator = gas_grid.decomposition().spec().communicator;

    for (const auto& newly_solid : global_newly_solid_voxels_) {
        if (newly_solid.previous_state != GasState::OutsideAccessible) {
            continue;
        }
        for (const auto& offset : neighbor_offsets) {
            const VoxelCoord candidate{
                newly_solid.voxel_coord.x + offset.x,
                newly_solid.voxel_coord.y + offset.y,
                newly_solid.voxel_coord.z + offset.z};
            const auto neighbor = normalized_neighbor(candidate, grid_spec);
            if (neighbor.has_value()
                && *neighbor != newly_solid.voxel_coord) {
                seed_voxels_.push_back(*neighbor);
            }
        }
    }
    sort_unique(seed_voxels_);

    for (const auto& seed_voxel : seed_voxels_) {
        begin_search_epoch();
        frontier_.clear();
        visited_voxels_.clear();
        bool local_source_found = false;

        const auto visit_owned = [&](const VoxelCoord& voxel_coord) {
            if (gas_grid.gas_state(voxel_coord)
                != GasState::OutsideAccessible) {
                return;
            }
            const auto index = owned_index(voxel_coord, owned_range);
            if (confirmed_outside_epochs_[index] == repair_epoch_) {
                local_source_found = true;
                return;
            }
            if (search_epochs_[index] == search_epoch_) {
                return;
            }
            search_epochs_[index] = search_epoch_;
            frontier_.push_back(voxel_coord);
            visited_voxels_.push_back(voxel_coord);
            ++result.local_visited_voxel_count;
            if (is_reservoir_source(voxel_coord, grid_spec)) {
                local_source_found = true;
            }
        };

        if (gas_grid.owns(seed_voxel)) {
            visit_owned(seed_voxel);
        }
        bool global_active = false;
        bool global_source_found = false;
        reduce_search_state(
            communicator,
            !frontier_.empty(),
            local_source_found,
            global_active,
            global_source_found);
        if (!global_active && !global_source_found) {
            continue;
        }
        ++result.seed_search_count;

        while (global_active && !global_source_found) {
            for (auto& send_buffer : send_buffers_) {
                send_buffer.clear();
            }

            std::size_t frontier_head = 0;
            while (frontier_head < frontier_.size()) {
                const auto current = frontier_[frontier_head];
                ++frontier_head;
                for (std::size_t neighbor_index = 0;
                     neighbor_index < neighbor_offsets.size();
                     ++neighbor_index) {
                    const auto& offset = neighbor_offsets[neighbor_index];
                    const VoxelCoord candidate{
                        current.x + offset.x,
                        current.y + offset.y,
                        current.z + offset.z};
                    const auto neighbor = normalized_neighbor(candidate, grid_spec);
                    if (!neighbor.has_value() || *neighbor == current) {
                        continue;
                    }
                    if (gas_grid.owns(*neighbor)) {
                        visit_owned(*neighbor);
                        continue;
                    }
                    const auto face = faces[neighbor_index];
                    if (gas_grid.decomposition().neighbor_rank(face)
                        != MPI_PROC_NULL) {
                        send_buffers_[face_index(face)].push_back(
                            encode_face_offset(face, current, owned_range));
                    }
                }
            }
            frontier_.clear();

            exchange_frontiers(
                gas_grid.decomposition(),
                send_buffers_,
                receive_buffers_,
                result);
            ++result.communication_round_count;
            for (const auto face : faces) {
                for (const auto face_offset :
                     receive_buffers_[face_index(face)]) {
                    visit_owned(decode_face_offset(
                        face,
                        face_offset,
                        owned_range));
                }
            }

            reduce_search_state(
                communicator,
                !frontier_.empty(),
                local_source_found,
                global_active,
                global_source_found);
        }

        if (global_source_found) {
            for (const auto& voxel_coord : visited_voxels_) {
                confirmed_outside_epochs_[owned_index(
                    voxel_coord,
                    owned_range)] = repair_epoch_;
            }
            continue;
        }

        for (const auto& voxel_coord : visited_voxels_) {
            gas_grid.set_owned_gas_state(voxel_coord, GasState::ClosedVoid);
            result.newly_closed_owned_voxel_coords.push_back(voxel_coord);
            ++result.local_closed_voxel_count;
        }
    }

    const std::uint64_t local_participated =
        result.local_visited_voxel_count == 0 ? 0U : 1U;
    check_mpi(
        MPI_Allreduce(
            &local_participated,
            &result.participating_rank_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            communicator),
        "MPI_Allreduce(closing participating ranks)");
    sort_unique(result.newly_closed_owned_voxel_coords);
    return result;
}

void DistributedClosingRegionRepair::gather_newly_solid_voxels(
    const DistributedGasGrid& gas_grid,
    NewlySolidVoxelCoordView newly_solid_voxel_view)
{
    const auto count_limit = static_cast<std::uint64_t>(INT_MAX)
        / encoded_change_width;
    if (newly_solid_voxel_view.count
        > static_cast<std::size_t>(count_limit)) {
        throw std::overflow_error(
            "too many local newly solid voxels for MPI gather");
    }

    std::vector<std::int64_t> local_values;
    local_values.reserve(
        newly_solid_voxel_view.count * encoded_change_width);
    for (std::size_t index = 0;
         index < newly_solid_voxel_view.count;
         ++index) {
        const auto& voxel = newly_solid_voxel_view.voxels[index];
        local_values.push_back(voxel.voxel_coord.x);
        local_values.push_back(voxel.voxel_coord.y);
        local_values.push_back(voxel.voxel_coord.z);
        local_values.push_back(
            static_cast<std::int64_t>(voxel.previous_state));
    }

    const auto rank_count = static_cast<std::size_t>(
        gas_grid.decomposition().size());
    const int local_value_count = static_cast<int>(local_values.size());
    std::vector<int> receive_counts(rank_count, 0);
    check_mpi(
        MPI_Allgather(
            &local_value_count,
            1,
            MPI_INT,
            receive_counts.data(),
            1,
            MPI_INT,
            gas_grid.decomposition().spec().communicator),
        "MPI_Allgather(newly solid voxel counts)");

    std::vector<int> displacements(rank_count, 0);
    std::uint64_t total_value_count = 0;
    for (std::size_t rank_index = 0;
         rank_index < rank_count;
         ++rank_index) {
        if (receive_counts[rank_index] < 0
            || receive_counts[rank_index]
                % static_cast<int>(encoded_change_width) != 0
            || total_value_count > static_cast<std::uint64_t>(INT_MAX)) {
            throw std::overflow_error("invalid newly-solid gather count");
        }
        displacements[rank_index] = static_cast<int>(total_value_count);
        total_value_count += static_cast<std::uint64_t>(
            receive_counts[rank_index]);
    }
    if (total_value_count > static_cast<std::uint64_t>(INT_MAX)) {
        throw std::overflow_error("newly-solid gather exceeds MPI int count");
    }

    std::vector<std::int64_t> global_values(
        static_cast<std::size_t>(total_value_count));
    check_mpi(
        MPI_Allgatherv(
            local_values.data(),
            local_value_count,
            MPI_INT64_T,
            global_values.data(),
            receive_counts.data(),
            displacements.data(),
            MPI_INT64_T,
            gas_grid.decomposition().spec().communicator),
        "MPI_Allgatherv(newly solid voxels)");

    global_newly_solid_voxels_.clear();
    global_newly_solid_voxels_.reserve(
        global_values.size() / encoded_change_width);
    for (std::size_t index = 0;
         index < global_values.size();
         index += encoded_change_width) {
        const auto state = global_values[index + 3];
        if (state != static_cast<std::int64_t>(GasState::OutsideAccessible)
            && state != static_cast<std::int64_t>(GasState::ClosedVoid)) {
            throw std::runtime_error(
                "received invalid newly-solid previous state");
        }
        global_newly_solid_voxels_.push_back({
            {global_values[index],
             global_values[index + 1],
             global_values[index + 2]},
            static_cast<GasState>(state)});
    }
    std::sort(
        global_newly_solid_voxels_.begin(),
        global_newly_solid_voxels_.end(),
        [](const DistributedRemovedVoxel& lhs,
           const DistributedRemovedVoxel& rhs) {
            return voxel_coord_less(lhs.voxel_coord, rhs.voxel_coord);
        });
}

void DistributedClosingRegionRepair::prepare_workspace(
    const DistributedGasGrid& gas_grid)
{
    const auto voxel_count = static_cast<std::size_t>(
        gas_grid.owned_voxel_count());
    if (search_epochs_.size() == voxel_count
        && confirmed_outside_epochs_.size() == voxel_count) {
        return;
    }
    search_epochs_.assign(voxel_count, 0);
    confirmed_outside_epochs_.assign(voxel_count, 0);
    search_epoch_ = 0;
    repair_epoch_ = 0;
}

void DistributedClosingRegionRepair::begin_repair_epoch()
{
    ++repair_epoch_;
    if (repair_epoch_ != 0) {
        return;
    }
    std::fill(
        confirmed_outside_epochs_.begin(),
        confirmed_outside_epochs_.end(),
        0);
    repair_epoch_ = 1;
}

void DistributedClosingRegionRepair::begin_search_epoch()
{
    ++search_epoch_;
    if (search_epoch_ != 0) {
        return;
    }
    std::fill(search_epochs_.begin(), search_epochs_.end(), 0);
    search_epoch_ = 1;
}

}  // namespace gasaccess
