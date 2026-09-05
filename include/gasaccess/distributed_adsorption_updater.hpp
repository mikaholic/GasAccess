#ifndef GASACCESS_DISTRIBUTED_ADSORPTION_UPDATER_HPP
#define GASACCESS_DISTRIBUTED_ADSORPTION_UPDATER_HPP

#include "gasaccess/adsorption_updater.hpp"
#include "gasaccess/distributed_exterior_classifier.hpp"
#include "gasaccess/mpi_gas_grid.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace gasaccess {

struct DistributedAdsorptionUpdateResult {
    std::uint64_t local_newly_solid_count = 0;
    std::uint64_t global_newly_solid_count = 0;
    std::vector<VoxelCoord> changed_owned_voxel_coords{};
    DistributedClassificationSummary classification{};
    bool topology_filter_performed = false;
    bool topology_filter_declared_safe = false;
    bool full_reclassification_performed = false;
    bool distributed_repair_performed = false;
    std::uint64_t topology_visited_voxel_count = 0;
    std::uint64_t local_repair_visited_voxel_count = 0;
    std::uint64_t local_repair_closed_voxel_count = 0;
    std::uint64_t repair_seed_search_count = 0;
    std::uint64_t repair_communication_round_count = 0;
    std::uint64_t sent_frontier_entry_count = 0;
    std::uint64_t received_frontier_entry_count = 0;

    bool geometry_changed() const noexcept;
    bool used_full_reclassification() const noexcept;
    bool used_distributed_repair() const noexcept;
};

class DistributedAdsorptionUpdater {
public:
    explicit DistributedAdsorptionUpdater(
        double precursor_radius,
        ConnectivityRepairMode repair_mode =
            ConnectivityRepairMode::AffectedRegion);

    double precursor_radius() const noexcept;
    ConnectivityRepairMode repair_mode() const noexcept;

    // Collective over the grid communicator. Adsorbed atoms must already be
    // synchronized far enough to cover every owned voxel within
    // R_atom + R_precursor. The grid must already be fully classified.
    DistributedAdsorptionUpdateResult apply_adsorption(
        DistributedGasGrid& gas_grid,
        AtomView adsorbed_atoms);

private:
    void gather_removed_voxels(const DistributedGasGrid& gas_grid);
    void prepare_repair_workspace(const DistributedGasGrid& gas_grid);
    void begin_repair_epoch();
    void begin_search_epoch();
    void repair_affected_regions(
        DistributedGasGrid& gas_grid,
        DistributedAdsorptionUpdateResult& result);

    double precursor_radius_ = 0.0;
    ConnectivityRepairMode repair_mode_ = ConnectivityRepairMode::AffectedRegion;
    DistributedExteriorClassifier exterior_classifier_{};
    std::vector<DistributedRemovedVoxel> removed_voxels_{};
    std::vector<DistributedRemovedVoxel> global_removed_voxels_{};
    std::vector<GasState> previous_states_{};
    std::vector<std::uint32_t> search_epochs_{};
    std::vector<std::uint32_t> confirmed_outside_epochs_{};
    std::vector<VoxelCoord> seed_voxels_{};
    std::vector<VoxelCoord> frontier_{};
    std::vector<VoxelCoord> visited_voxels_{};
    std::array<std::vector<std::uint64_t>, 6> send_buffers_{};
    std::array<std::vector<std::uint64_t>, 6> receive_buffers_{};
    std::uint32_t search_epoch_ = 0;
    std::uint32_t repair_epoch_ = 0;
};

}  // namespace gasaccess

#endif
