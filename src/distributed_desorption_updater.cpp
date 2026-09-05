#include "gasaccess/distributed_desorption_updater.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>

namespace gasaccess {
namespace {

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

void sort_changed_voxels(std::vector<VoxelCoord>& voxel_coords)
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

DistributedClassificationSummary current_classification(
    const DistributedGasGrid& gas_grid)
{
    const int local_unclassified =
        gas_grid.owned_gas_state_count(GasState::Unclassified) == 0 ? 0 : 1;
    int global_unclassified = 0;
    check_mpi(
        MPI_Allreduce(
            &local_unclassified,
            &global_unclassified,
            1,
            MPI_INT,
            MPI_MAX,
            gas_grid.decomposition().spec().communicator),
        "MPI_Allreduce(desorption classification validation)");
    if (global_unclassified != 0) {
        throw std::invalid_argument(
            "distributed desorption requires a fully classified gas grid");
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

AtomChangeOccupancyResult apply_occupancy_collectively(
    DistributedGasGrid& gas_grid,
    AtomView removed_atoms,
    double precursor_radius,
    std::vector<DistributedVoxelOccupancyChange>& occupancy_changes)
{
    occupancy_changes.clear();
    AtomChangeOccupancyResult occupancy_result{};
    std::exception_ptr local_error;
    try {
        occupancy_result = gas_grid.apply_owned_atom_changes(
            {{nullptr, 0}, removed_atoms},
            precursor_radius,
            occupancy_changes);
        if (occupancy_result.newly_solid_count != 0) {
            throw std::logic_error(
                "removal-only occupancy update created a solid voxel");
        }
    } catch (...) {
        local_error = std::current_exception();
    }

    const int local_failed = local_error == nullptr ? 0 : 1;
    int global_failed = 0;
    check_mpi(
        MPI_Allreduce(
            &local_failed,
            &global_failed,
            1,
            MPI_INT,
            MPI_MAX,
            gas_grid.decomposition().spec().communicator),
        "MPI_Allreduce(desorption occupancy validation)");
    if (global_failed == 0) {
        return occupancy_result;
    }

    if (!occupancy_changes.empty()) {
        for (const auto& change : occupancy_changes) {
            gas_grid.set_owned_blocker_count(
                change.voxel_coord,
                change.previous_blocker_count);
            gas_grid.set_owned_gas_state(
                change.voxel_coord,
                change.previous_state);
        }
    }
    occupancy_changes.clear();
    throw std::invalid_argument(
        "distributed desorption occupancy update failed on at least one rank");
}

}  // namespace

bool DistributedDesorptionUpdateResult::geometry_changed() const noexcept
{
    return global_newly_gas_count != 0;
}

bool DistributedDesorptionUpdateResult::used_full_reclassification() const noexcept
{
    return full_reclassification_performed;
}

bool DistributedDesorptionUpdateResult::used_distributed_opening_repair()
    const noexcept
{
    return distributed_opening_repair_performed;
}

DistributedDesorptionUpdater::DistributedDesorptionUpdater(
    double precursor_radius,
    DesorptionRepairMode repair_mode)
    : precursor_radius_(precursor_radius),
      repair_mode_(repair_mode)
{
    if (precursor_radius_ < 0.0
        || !std::isfinite(precursor_radius_)) {
        throw std::invalid_argument(
            "precursor radius must be finite and nonnegative");
    }
    if (repair_mode_ != DesorptionRepairMode::OpeningRegion
        && repair_mode_ != DesorptionRepairMode::FullReclassification) {
        throw std::invalid_argument("invalid desorption repair mode");
    }
}

double DistributedDesorptionUpdater::precursor_radius() const noexcept
{
    return precursor_radius_;
}

DesorptionRepairMode DistributedDesorptionUpdater::repair_mode() const noexcept
{
    return repair_mode_;
}

DistributedDesorptionUpdateResult
DistributedDesorptionUpdater::apply_desorption(
    DistributedGasGrid& gas_grid,
    AtomView removed_atoms)
{
    DistributedDesorptionUpdateResult result{};
    result.classification = current_classification(gas_grid);
    const auto& owned_range = gas_grid.owned_range();

    if (repair_mode_ == DesorptionRepairMode::FullReclassification) {
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

    const auto occupancy_result = apply_occupancy_collectively(
        gas_grid,
        removed_atoms,
        precursor_radius_,
        occupancy_changes_);
    result.local_blocker_count_changed_voxel_count =
        occupancy_result.blocker_count_changed_voxel_count;
    result.local_newly_gas_count = occupancy_result.newly_gas_count;

    const std::array<std::uint64_t, 2> local_counts{
        result.local_blocker_count_changed_voxel_count,
        result.local_newly_gas_count};
    std::array<std::uint64_t, 2> global_counts{};
    check_mpi(
        MPI_Allreduce(
            local_counts.data(),
            global_counts.data(),
            static_cast<int>(global_counts.size()),
            MPI_UINT64_T,
            MPI_SUM,
            gas_grid.decomposition().spec().communicator),
        "MPI_Allreduce(desorption occupancy counts)");
    result.global_blocker_count_changed_voxel_count = global_counts[0];
    result.global_newly_gas_count = global_counts[1];
    if (!result.geometry_changed()) {
        return result;
    }

    newly_gas_voxel_coords_.clear();
    newly_gas_voxel_coords_.reserve(
        static_cast<std::size_t>(result.local_newly_gas_count));
    for (const auto& occupancy_change : occupancy_changes_) {
        if (occupancy_change.previous_blocker_count != 0
            && occupancy_change.blocker_count == 0) {
            newly_gas_voxel_coords_.push_back(occupancy_change.voxel_coord);
        }
    }
    if (newly_gas_voxel_coords_.size()
        != static_cast<std::size_t>(result.local_newly_gas_count)) {
        throw std::logic_error(
            "distributed voxelizer gas-count result is inconsistent");
    }

    if (repair_mode_ == DesorptionRepairMode::FullReclassification) {
        result.full_reclassification_performed = true;
        result.classification = exterior_classifier_.classify(gas_grid);
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

    result.changed_owned_voxel_coords = newly_gas_voxel_coords_;
    for (const auto& voxel_coord : newly_gas_voxel_coords_) {
        gas_grid.set_owned_gas_state(voxel_coord, GasState::ClosedVoid);
    }

    auto repair_result = opening_region_repair_.repair(
        gas_grid,
        {newly_gas_voxel_coords_.data(), newly_gas_voxel_coords_.size()});
    result.local_opening_visited_voxel_count =
        repair_result.local_visited_voxel_count;
    result.local_repair_opened_voxel_count =
        repair_result.local_opened_voxel_count;
    result.opening_participating_rank_count =
        repair_result.participating_rank_count;
    result.repair_communication_round_count =
        repair_result.communication_round_count;
    result.sent_frontier_entry_count =
        repair_result.sent_frontier_entry_count;
    result.received_frontier_entry_count =
        repair_result.received_frontier_entry_count;
    result.distributed_opening_repair_performed =
        result.opening_participating_rank_count != 0;
    result.changed_owned_voxel_coords.insert(
        result.changed_owned_voxel_coords.end(),
        repair_result.newly_opened_owned_voxel_coords.begin(),
        repair_result.newly_opened_owned_voxel_coords.end());
    sort_changed_voxels(result.changed_owned_voxel_coords);

    gas_grid.exchange_ghost_states();
    result.classification = current_classification(gas_grid);
    return result;
}

}  // namespace gasaccess
