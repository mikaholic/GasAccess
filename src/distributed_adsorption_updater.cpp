#include "gasaccess/distributed_adsorption_updater.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
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

constexpr int count_tag_base = 5510;
constexpr int payload_tag_base = 5520;
constexpr std::size_t local_voxel_capacity = 27;
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
        return coordinate < 0 ? signed_dimension - 1 : 0;
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
    normalized = {*x, *y, *z};
    return normalized;
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

bool contains_coord(
    const VoxelCoord* voxel_coords,
    std::size_t count,
    const VoxelCoord& voxel_coord)
{
    return std::find(
        voxel_coords,
        voxel_coords + static_cast<std::ptrdiff_t>(count),
        voxel_coord) != voxel_coords + static_cast<std::ptrdiff_t>(count);
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

DistributedClassificationSummary current_classification(
    const DistributedGasGrid& gas_grid)
{
    if (gas_grid.owned_gas_state_count(GasState::Unclassified) != 0) {
        throw std::invalid_argument(
            "distributed adsorption requires a fully classified gas grid");
    }

    DistributedClassificationSummary summary{};
    summary.local_solid_count =
        gas_grid.owned_gas_state_count(GasState::Solid);
    summary.local_outside_accessible_count =
        gas_grid.owned_gas_state_count(GasState::OutsideAccessible);
    summary.local_closed_void_count =
        gas_grid.owned_gas_state_count(GasState::ClosedVoid);
    return summary;
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

TopologyCheckResult evaluate_local_topology(
    const DistributedGasGrid& gas_grid,
    const DistributedRemovedVoxel& removed_voxel)
{
    if (!gas_grid.owns(removed_voxel.voxel_coord)
        || gas_grid.gas_state(removed_voxel.voxel_coord) != GasState::Solid) {
        throw std::invalid_argument("removed voxel is not owned solid occupancy");
    }
    if (removed_voxel.previous_state != GasState::OutsideAccessible
        && removed_voxel.previous_state != GasState::ClosedVoid) {
        throw std::invalid_argument(
            "removed voxel was not classified as empty before adsorption");
    }
    if (removed_voxel.previous_state == GasState::ClosedVoid) {
        return {TopologyDecision::Safe, 0, 0};
    }

    const auto& grid_spec = gas_grid.global_grid_spec();
    if (is_reservoir_source(removed_voxel.voxel_coord, grid_spec)) {
        return {TopologyDecision::RequiresConnectivityRepair, 0, 0};
    }

    std::array<VoxelCoord, 6> accessible_neighbors{};
    std::size_t accessible_neighbor_count = 0;
    for (const auto& offset : neighbor_offsets) {
        const VoxelCoord candidate{
            removed_voxel.voxel_coord.x + offset.x,
            removed_voxel.voxel_coord.y + offset.y,
            removed_voxel.voxel_coord.z + offset.z};
        const auto neighbor = normalized_neighbor(candidate, grid_spec);
        if (!neighbor.has_value()
            || *neighbor == removed_voxel.voxel_coord
            || gas_grid.gas_state(*neighbor) != GasState::OutsideAccessible
            || contains_coord(
                accessible_neighbors.data(),
                accessible_neighbor_count,
                *neighbor)) {
            continue;
        }
        accessible_neighbors[accessible_neighbor_count] = *neighbor;
        ++accessible_neighbor_count;
    }
    if (accessible_neighbor_count <= 1) {
        return {
            TopologyDecision::Safe,
            accessible_neighbor_count,
            accessible_neighbor_count};
    }

    std::array<VoxelCoord, local_voxel_capacity> neighborhood{};
    std::size_t neighborhood_count = 0;
    for (std::int64_t z_offset = -1; z_offset <= 1; ++z_offset) {
        for (std::int64_t y_offset = -1; y_offset <= 1; ++y_offset) {
            for (std::int64_t x_offset = -1; x_offset <= 1; ++x_offset) {
                const VoxelCoord candidate{
                    removed_voxel.voxel_coord.x + x_offset,
                    removed_voxel.voxel_coord.y + y_offset,
                    removed_voxel.voxel_coord.z + z_offset};
                const auto normalized = normalized_neighbor(candidate, grid_spec);
                if (!normalized.has_value()
                    || contains_coord(
                        neighborhood.data(),
                        neighborhood_count,
                        *normalized)) {
                    continue;
                }
                if (!gas_grid.owns(*normalized)) {
                    return {
                        TopologyDecision::RequiresConnectivityRepair,
                        accessible_neighbor_count,
                        0};
                }
                neighborhood[neighborhood_count] = *normalized;
                ++neighborhood_count;
            }
        }
    }

    std::array<VoxelCoord, local_voxel_capacity> visited{};
    visited[0] = accessible_neighbors[0];
    std::size_t visited_count = 1;
    std::size_t frontier_index_value = 0;
    while (frontier_index_value < visited_count) {
        const auto current = visited[frontier_index_value];
        ++frontier_index_value;
        for (const auto& offset : neighbor_offsets) {
            const VoxelCoord candidate{
                current.x + offset.x,
                current.y + offset.y,
                current.z + offset.z};
            const auto neighbor = normalized_neighbor(candidate, grid_spec);
            if (!neighbor.has_value()
                || !contains_coord(
                    neighborhood.data(),
                    neighborhood_count,
                    *neighbor)
                || gas_grid.gas_state(*neighbor)
                    != GasState::OutsideAccessible
                || contains_coord(visited.data(), visited_count, *neighbor)) {
                continue;
            }
            if (visited_count >= visited.size()) {
                throw std::logic_error(
                    "distributed local topology search exceeded capacity");
            }
            visited[visited_count] = *neighbor;
            ++visited_count;
        }
    }

    for (std::size_t index = 1; index < accessible_neighbor_count; ++index) {
        if (!contains_coord(
            visited.data(),
            visited_count,
            accessible_neighbors[index])) {
            return {
                TopologyDecision::RequiresConnectivityRepair,
                accessible_neighbor_count,
                visited_count};
        }
    }
    return {
        TopologyDecision::Safe,
        accessible_neighbor_count,
        visited_count};
}

std::uint64_t face_element_count(
    Face face,
    const OwnedVoxelRange& owned_range)
{
    switch (face) {
    case Face::XLow:
    case Face::XHigh:
        return owned_range.dimensions.y * owned_range.dimensions.z;
    case Face::YLow:
    case Face::YHigh:
        return owned_range.dimensions.x * owned_range.dimensions.z;
    case Face::ZLow:
    case Face::ZHigh:
        return owned_range.dimensions.x * owned_range.dimensions.y;
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
        throw std::runtime_error("received repair offset exceeds owned face");
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
        throw std::overflow_error("MPI repair payload exceeds INT_MAX entries");
    }
    return static_cast<int>(count);
}

void exchange_frontiers(
    const MpiDecomposition& decomposition,
    std::array<std::vector<std::uint64_t>, 6>& send_buffers,
    std::array<std::vector<std::uint64_t>, 6>& receive_buffers,
    DistributedAdsorptionUpdateResult& result)
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
            "MPI_Irecv(repair count)");
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
            "MPI_Isend(repair count)");
        ++request_count;
    }
    if (request_count != 0) {
        check_mpi(
            MPI_Waitall(request_count, requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall(repair counts)");
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
                "MPI_Irecv(repair payload)");
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
                "MPI_Isend(repair payload)");
            ++request_count;
        }
    }
    if (request_count != 0) {
        check_mpi(
            MPI_Waitall(request_count, requests.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall(repair payloads)");
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
        "MPI_Allreduce(repair search state)");
    global_active = global_values[0] != 0;
    global_source_found = global_values[1] != 0;
}

void sort_changed_voxels(std::vector<VoxelCoord>& voxel_coords)
{
    std::sort(voxel_coords.begin(), voxel_coords.end(), voxel_coord_less);
    voxel_coords.erase(
        std::unique(voxel_coords.begin(), voxel_coords.end()),
        voxel_coords.end());
}

}  // namespace

bool DistributedAdsorptionUpdateResult::geometry_changed() const noexcept
{
    return global_newly_solid_count != 0;
}

bool DistributedAdsorptionUpdateResult::used_full_reclassification() const noexcept
{
    return full_reclassification_performed;
}

bool DistributedAdsorptionUpdateResult::used_distributed_repair() const noexcept
{
    return distributed_repair_performed;
}

DistributedAdsorptionUpdater::DistributedAdsorptionUpdater(
    double precursor_radius,
    ConnectivityRepairMode repair_mode)
    : precursor_radius_(precursor_radius),
      repair_mode_(repair_mode)
{
    if (!std::isfinite(precursor_radius_) || precursor_radius_ < 0.0) {
        throw std::invalid_argument(
            "precursor radius must be finite and nonnegative");
    }
    if (repair_mode_ != ConnectivityRepairMode::AffectedRegion
        && repair_mode_ != ConnectivityRepairMode::FullReclassification) {
        throw std::invalid_argument("invalid connectivity repair mode");
    }
}

double DistributedAdsorptionUpdater::precursor_radius() const noexcept
{
    return precursor_radius_;
}

ConnectivityRepairMode DistributedAdsorptionUpdater::repair_mode() const noexcept
{
    return repair_mode_;
}

DistributedAdsorptionUpdateResult
DistributedAdsorptionUpdater::apply_adsorption(
    DistributedGasGrid& gas_grid,
    AtomView adsorbed_atoms)
{
    DistributedAdsorptionUpdateResult result{};
    result.classification = current_classification(gas_grid);
    const auto& owned_range = gas_grid.owned_range();

    if (repair_mode_ == ConnectivityRepairMode::FullReclassification) {
        previous_states_.resize(
            static_cast<std::size_t>(gas_grid.owned_voxel_count()));
        for (auto z = owned_range.begin.z; z < owned_range.end.z; ++z) {
            for (auto y = owned_range.begin.y; y < owned_range.end.y; ++y) {
                for (auto x = owned_range.begin.x; x < owned_range.end.x; ++x) {
                    const VoxelCoord voxel_coord{x, y, z};
                    previous_states_[owned_index(voxel_coord, owned_range)] =
                        gas_grid.gas_state(voxel_coord);
                }
            }
        }
    }

    result.local_newly_solid_count = gas_grid.voxelize_owned_atoms(
        adsorbed_atoms,
        precursor_radius_,
        removed_voxels_);
    const auto communicator = gas_grid.decomposition().spec().communicator;
    check_mpi(
        MPI_Allreduce(
            &result.local_newly_solid_count,
            &result.global_newly_solid_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            communicator),
        "MPI_Allreduce(newly solid count)");
    if (!result.geometry_changed()) {
        return result;
    }

    result.changed_owned_voxel_coords.reserve(removed_voxels_.size());
    for (const auto& removed_voxel : removed_voxels_) {
        result.changed_owned_voxel_coords.push_back(removed_voxel.voxel_coord);
    }

    if (repair_mode_ == ConnectivityRepairMode::FullReclassification) {
        result.full_reclassification_performed = true;
        result.classification = exterior_classifier_.classify(gas_grid);
        result.changed_owned_voxel_coords.clear();
        for (auto z = owned_range.begin.z; z < owned_range.end.z; ++z) {
            for (auto y = owned_range.begin.y; y < owned_range.end.y; ++y) {
                for (auto x = owned_range.begin.x; x < owned_range.end.x; ++x) {
                    const VoxelCoord voxel_coord{x, y, z};
                    if (previous_states_[owned_index(voxel_coord, owned_range)]
                        != gas_grid.gas_state(voxel_coord)) {
                        result.changed_owned_voxel_coords.push_back(voxel_coord);
                    }
                }
            }
        }
        sort_changed_voxels(result.changed_owned_voxel_coords);
        return result;
    }

    std::uint64_t local_outside_removal_count = 0;
    for (const auto& removed_voxel : removed_voxels_) {
        if (removed_voxel.previous_state == GasState::OutsideAccessible) {
            ++local_outside_removal_count;
        }
    }
    std::uint64_t global_outside_removal_count = 0;
    check_mpi(
        MPI_Allreduce(
            &local_outside_removal_count,
            &global_outside_removal_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            communicator),
        "MPI_Allreduce(outside removal count)");

    bool requires_repair = global_outside_removal_count != 0;
    if (global_outside_removal_count != 0
        && result.global_newly_solid_count == 1) {
        result.topology_filter_performed = true;
        int local_requires_repair = 0;
        if (!removed_voxels_.empty()) {
            const auto topology_result = evaluate_local_topology(
                gas_grid,
                removed_voxels_.front());
            result.topology_visited_voxel_count = static_cast<std::uint64_t>(
                topology_result.visited_voxel_count);
            local_requires_repair = topology_result.is_safe() ? 0 : 1;
        }
        int global_requires_repair = 0;
        check_mpi(
            MPI_Allreduce(
                &local_requires_repair,
                &global_requires_repair,
                1,
                MPI_INT,
                MPI_MAX,
                communicator),
            "MPI_Allreduce(topology decision)");
        requires_repair = global_requires_repair != 0;
        result.topology_filter_declared_safe = !requires_repair;
    }

    if (requires_repair) {
        gather_removed_voxels(gas_grid);
        result.distributed_repair_performed = true;
        repair_affected_regions(gas_grid, result);
    }

    gas_grid.exchange_ghost_states();
    result.classification = current_classification(gas_grid);
    sort_changed_voxels(result.changed_owned_voxel_coords);
    return result;
}

void DistributedAdsorptionUpdater::gather_removed_voxels(
    const DistributedGasGrid& gas_grid)
{
    const auto global_count_limit = static_cast<std::uint64_t>(INT_MAX)
        / encoded_change_width;
    if (removed_voxels_.size() > static_cast<std::size_t>(global_count_limit)) {
        throw std::overflow_error("too many local removed voxels for MPI gather");
    }

    std::vector<std::int64_t> local_values;
    local_values.reserve(removed_voxels_.size() * encoded_change_width);
    for (const auto& removed_voxel : removed_voxels_) {
        local_values.push_back(removed_voxel.voxel_coord.x);
        local_values.push_back(removed_voxel.voxel_coord.y);
        local_values.push_back(removed_voxel.voxel_coord.z);
        local_values.push_back(
            static_cast<std::int64_t>(removed_voxel.previous_state));
    }

    const auto communicator = gas_grid.decomposition().spec().communicator;
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
            communicator),
        "MPI_Allgather(removed voxel counts)");

    std::vector<int> displacements(rank_count, 0);
    std::uint64_t total_value_count = 0;
    for (std::size_t rank_index = 0;
         rank_index < rank_count;
         ++rank_index) {
        if (receive_counts[rank_index] < 0
            || receive_counts[rank_index]
                % static_cast<int>(encoded_change_width) != 0
            || total_value_count > static_cast<std::uint64_t>(INT_MAX)) {
            throw std::overflow_error("invalid removed-voxel gather count");
        }
        displacements[rank_index] = static_cast<int>(total_value_count);
        total_value_count += static_cast<std::uint64_t>(
            receive_counts[rank_index]);
    }
    if (total_value_count > static_cast<std::uint64_t>(INT_MAX)) {
        throw std::overflow_error("removed-voxel gather exceeds MPI int count");
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
            communicator),
        "MPI_Allgatherv(removed voxels)");

    global_removed_voxels_.clear();
    global_removed_voxels_.reserve(
        global_values.size() / encoded_change_width);
    for (std::size_t value_index = 0;
         value_index < global_values.size();
         value_index += encoded_change_width) {
        const auto state_value = global_values[value_index + 3];
        if (state_value != static_cast<std::int64_t>(GasState::OutsideAccessible)
            && state_value != static_cast<std::int64_t>(GasState::ClosedVoid)) {
            throw std::runtime_error("received invalid previous gas state");
        }
        global_removed_voxels_.push_back({
            {global_values[value_index],
             global_values[value_index + 1],
             global_values[value_index + 2]},
            static_cast<GasState>(state_value)});
    }
    std::sort(
        global_removed_voxels_.begin(),
        global_removed_voxels_.end(),
        [](const DistributedRemovedVoxel& lhs,
           const DistributedRemovedVoxel& rhs) {
            return voxel_coord_less(lhs.voxel_coord, rhs.voxel_coord);
        });
}

void DistributedAdsorptionUpdater::prepare_repair_workspace(
    const DistributedGasGrid& gas_grid)
{
    const auto local_voxel_count = static_cast<std::size_t>(
        gas_grid.owned_voxel_count());
    if (search_epochs_.size() == local_voxel_count
        && confirmed_outside_epochs_.size() == local_voxel_count) {
        return;
    }
    search_epochs_.assign(local_voxel_count, 0);
    confirmed_outside_epochs_.assign(local_voxel_count, 0);
    search_epoch_ = 0;
    repair_epoch_ = 0;
}

void DistributedAdsorptionUpdater::begin_repair_epoch()
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

void DistributedAdsorptionUpdater::begin_search_epoch()
{
    ++search_epoch_;
    if (search_epoch_ != 0) {
        return;
    }
    std::fill(search_epochs_.begin(), search_epochs_.end(), 0);
    search_epoch_ = 1;
}

void DistributedAdsorptionUpdater::repair_affected_regions(
    DistributedGasGrid& gas_grid,
    DistributedAdsorptionUpdateResult& result)
{
    prepare_repair_workspace(gas_grid);
    begin_repair_epoch();
    seed_voxels_.clear();
    const auto& grid_spec = gas_grid.global_grid_spec();
    const auto& owned_range = gas_grid.owned_range();
    const auto communicator = gas_grid.decomposition().spec().communicator;

    for (const auto& removed_voxel : global_removed_voxels_) {
        if (removed_voxel.previous_state != GasState::OutsideAccessible) {
            continue;
        }
        for (const auto& offset : neighbor_offsets) {
            const VoxelCoord candidate{
                removed_voxel.voxel_coord.x + offset.x,
                removed_voxel.voxel_coord.y + offset.y,
                removed_voxel.voxel_coord.z + offset.z};
            const auto neighbor = normalized_neighbor(candidate, grid_spec);
            if (neighbor.has_value()
                && *neighbor != removed_voxel.voxel_coord) {
                seed_voxels_.push_back(*neighbor);
            }
        }
    }
    sort_changed_voxels(seed_voxels_);

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
            ++result.local_repair_visited_voxel_count;
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
        ++result.repair_seed_search_count;

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
                    const auto offset = neighbor_offsets[neighbor_index];
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
            ++result.repair_communication_round_count;
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
            result.changed_owned_voxel_coords.push_back(voxel_coord);
            ++result.local_repair_closed_voxel_count;
        }
    }
}

}  // namespace gasaccess
