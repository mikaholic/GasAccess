#include "gasaccess/distributed_atom_change_updater.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

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

DistributedClassificationSummary local_classification(
    const DistributedGasGrid& gas_grid)
{
    DistributedClassificationSummary summary{};
    summary.local_solid_count =
        gas_grid.owned_gas_state_count(GasState::Solid);
    summary.local_outside_accessible_count =
        gas_grid.owned_gas_state_count(GasState::OutsideAccessible);
    summary.local_closed_void_count =
        gas_grid.owned_gas_state_count(GasState::ClosedVoid);
    return summary;
}

DistributedClassificationSummary validated_classification(
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
        "MPI_Allreduce(atom-change classification validation)");
    if (global_unclassified != 0) {
        throw std::invalid_argument(
            "distributed atom-change update requires a fully classified grid");
    }
    return local_classification(gas_grid);
}

AtomChangeOccupancyResult apply_occupancy_collectively(
    DistributedGasGrid& gas_grid,
    const AtomChangeBatch& atom_changes,
    double precursor_radius,
    std::vector<DistributedVoxelOccupancyChange>& occupancy_changes)
{
    occupancy_changes.clear();
    AtomChangeOccupancyResult result{};
    std::exception_ptr local_error;
    try {
        result = gas_grid.apply_owned_atom_changes(
            atom_changes,
            precursor_radius,
            occupancy_changes);
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
        "MPI_Allreduce(atom-change occupancy validation)");
    if (global_failed == 0) {
        return result;
    }

    for (const auto& change : occupancy_changes) {
        gas_grid.set_owned_blocker_count(
            change.voxel_coord,
            change.previous_blocker_count);
        gas_grid.set_owned_gas_state(
            change.voxel_coord,
            change.previous_state);
    }
    occupancy_changes.clear();
    throw std::invalid_argument(
        "distributed atom-change occupancy failed on at least one rank");
}

AccessibilityRepairKind incremental_repair_kind(
    std::uint64_t newly_solid_count,
    std::uint64_t newly_gas_count) noexcept
{
    if (newly_solid_count != 0 && newly_gas_count != 0) {
        return AccessibilityRepairKind::Mixed;
    }
    if (newly_solid_count != 0) {
        return AccessibilityRepairKind::Closing;
    }
    if (newly_gas_count != 0) {
        return AccessibilityRepairKind::Opening;
    }
    return AccessibilityRepairKind::None;
}

}  // namespace

bool DistributedAtomChangeUpdateResult::geometry_changed() const noexcept
{
    return global_newly_solid_count != 0 || global_newly_gas_count != 0;
}

bool DistributedAtomChangeUpdateResult::used_full_reclassification()
    const noexcept
{
    return full_reclassification_performed;
}

bool DistributedAtomChangeUpdateResult::used_distributed_closing_repair()
    const noexcept
{
    return distributed_closing_repair_performed;
}

bool DistributedAtomChangeUpdateResult::used_distributed_opening_repair()
    const noexcept
{
    return distributed_opening_repair_performed;
}

DistributedAtomChangeUpdater::DistributedAtomChangeUpdater(
    double precursor_radius,
    AtomChangeRepairMode repair_mode)
    : precursor_radius_(precursor_radius),
      repair_mode_(repair_mode)
{
    if (!std::isfinite(precursor_radius_) || precursor_radius_ < 0.0) {
        throw std::invalid_argument(
            "precursor radius must be finite and nonnegative");
    }
    if (repair_mode_ != AtomChangeRepairMode::Incremental
        && repair_mode_ != AtomChangeRepairMode::FullReclassification) {
        throw std::invalid_argument("invalid atom-change repair mode");
    }
}

double DistributedAtomChangeUpdater::precursor_radius() const noexcept
{
    return precursor_radius_;
}

AtomChangeRepairMode DistributedAtomChangeUpdater::repair_mode() const noexcept
{
    return repair_mode_;
}

DistributedAtomChangeUpdateResult
DistributedAtomChangeUpdater::apply_atom_changes(
    DistributedGasGrid& gas_grid,
    const AtomChangeBatch& atom_changes)
{
    DistributedAtomChangeUpdateResult result{};
    result.classification = validated_classification(gas_grid);
    const auto& owned_range = gas_grid.owned_range();

    if (repair_mode_ == AtomChangeRepairMode::FullReclassification) {
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
        atom_changes,
        precursor_radius_,
        occupancy_changes_);
    result.local_blocker_count_changed_voxel_count =
        occupancy_result.blocker_count_changed_voxel_count;
    result.local_newly_solid_count = occupancy_result.newly_solid_count;
    result.local_newly_gas_count = occupancy_result.newly_gas_count;

    const std::array<std::uint64_t, 3> local_counts{
        result.local_blocker_count_changed_voxel_count,
        result.local_newly_solid_count,
        result.local_newly_gas_count};
    std::array<std::uint64_t, 3> global_counts{};
    check_mpi(
        MPI_Allreduce(
            local_counts.data(),
            global_counts.data(),
            static_cast<int>(global_counts.size()),
            MPI_UINT64_T,
            MPI_SUM,
            gas_grid.decomposition().spec().communicator),
        "MPI_Allreduce(atom-change occupancy counts)");
    result.global_blocker_count_changed_voxel_count = global_counts[0];
    result.global_newly_solid_count = global_counts[1];
    result.global_newly_gas_count = global_counts[2];
    if (!result.geometry_changed()) {
        return result;
    }

    if (repair_mode_ == AtomChangeRepairMode::FullReclassification) {
        result.repair_kind = AccessibilityRepairKind::FullReclassification;
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
        sort_unique(result.changed_owned_voxel_coords);
        return result;
    }

    result.repair_kind = incremental_repair_kind(
        result.global_newly_solid_count,
        result.global_newly_gas_count);
    newly_solid_voxels_.clear();
    newly_solid_voxels_.reserve(
        static_cast<std::size_t>(result.local_newly_solid_count));
    newly_gas_voxel_coords_.clear();
    newly_gas_voxel_coords_.reserve(
        static_cast<std::size_t>(result.local_newly_gas_count));
    for (const auto& change : occupancy_changes_) {
        if (change.previous_blocker_count == 0
            && change.blocker_count != 0) {
            newly_solid_voxels_.push_back({
                change.voxel_coord,
                change.previous_state});
            result.changed_owned_voxel_coords.push_back(change.voxel_coord);
        } else if (change.previous_blocker_count != 0
                   && change.blocker_count == 0) {
            newly_gas_voxel_coords_.push_back(change.voxel_coord);
            result.changed_owned_voxel_coords.push_back(change.voxel_coord);
        }
    }
    if (newly_solid_voxels_.size()
            != static_cast<std::size_t>(result.local_newly_solid_count)
        || newly_gas_voxel_coords_.size()
            != static_cast<std::size_t>(result.local_newly_gas_count)) {
        throw std::logic_error(
            "distributed atom-change transition counts are inconsistent");
    }

    if (result.global_newly_solid_count != 0) {
        auto closing_result = closing_repair_.repair(
            gas_grid,
            {newly_solid_voxels_.data(), newly_solid_voxels_.size()});
        result.local_closing_visited_voxel_count =
            closing_result.local_visited_voxel_count;
        result.local_repair_closed_voxel_count =
            closing_result.local_closed_voxel_count;
        result.closing_participating_rank_count =
            closing_result.participating_rank_count;
        result.closing_seed_search_count = closing_result.seed_search_count;
        result.closing_communication_round_count =
            closing_result.communication_round_count;
        result.closing_sent_frontier_entry_count =
            closing_result.sent_frontier_entry_count;
        result.closing_received_frontier_entry_count =
            closing_result.received_frontier_entry_count;
        result.distributed_closing_repair_performed =
            result.closing_participating_rank_count != 0;
        closing_changed_voxel_coords_ = std::move(
            closing_result.newly_closed_owned_voxel_coords);
    } else {
        closing_changed_voxel_coords_.clear();
    }

    if (result.global_newly_gas_count != 0) {
        for (const auto& voxel_coord : newly_gas_voxel_coords_) {
            gas_grid.set_owned_gas_state(voxel_coord, GasState::ClosedVoid);
        }
        const bool face_ghost_states_current =
            result.global_newly_solid_count == 0;
        const auto opening_result = opening_repair_.repair(
            gas_grid,
            {newly_gas_voxel_coords_.data(), newly_gas_voxel_coords_.size()},
            face_ghost_states_current);
        result.local_opening_visited_voxel_count =
            opening_result.local_visited_voxel_count;
        result.local_repair_opened_voxel_count =
            opening_result.local_opened_voxel_count;
        result.opening_participating_rank_count =
            opening_result.participating_rank_count;
        result.opening_communication_round_count =
            opening_result.communication_round_count;
        result.opening_sent_frontier_entry_count =
            opening_result.sent_frontier_entry_count;
        result.opening_received_frontier_entry_count =
            opening_result.received_frontier_entry_count;
        result.distributed_opening_repair_performed =
            result.opening_participating_rank_count != 0;
        for (const auto& voxel_coord :
             opening_result.newly_opened_owned_voxel_coords) {
            if (!std::binary_search(
                    closing_changed_voxel_coords_.begin(),
                    closing_changed_voxel_coords_.end(),
                    voxel_coord,
                    voxel_coord_less)) {
                result.changed_owned_voxel_coords.push_back(voxel_coord);
            }
        }
    }

    for (const auto& voxel_coord : closing_changed_voxel_coords_) {
        if (gas_grid.gas_state(voxel_coord) == GasState::ClosedVoid) {
            result.changed_owned_voxel_coords.push_back(voxel_coord);
        }
    }

    sort_unique(result.changed_owned_voxel_coords);
    gas_grid.exchange_ghost_states();
    result.classification = local_classification(gas_grid);
    return result;
}

DistributedAtomChangeUpdateResult
DistributedAtomChangeUpdater::apply_adsorption(
    DistributedGasGrid& gas_grid,
    AtomView added_atoms)
{
    return apply_atom_changes(gas_grid, {added_atoms, {nullptr, 0}});
}

DistributedAtomChangeUpdateResult
DistributedAtomChangeUpdater::apply_desorption(
    DistributedGasGrid& gas_grid,
    AtomView removed_atoms)
{
    return apply_atom_changes(gas_grid, {{nullptr, 0}, removed_atoms});
}

}  // namespace gasaccess
