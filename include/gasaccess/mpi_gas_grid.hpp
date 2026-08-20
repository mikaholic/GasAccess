#ifndef GASACCESS_MPI_GAS_GRID_HPP
#define GASACCESS_MPI_GAS_GRID_HPP

#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/gas_grid.hpp"

#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX 1
#endif
#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace gasaccess {

struct AlignedGridGeometry {
    GridSpacing spacing{};
    GridDimensions dimensions{};
};

AlignedGridGeometry make_aligned_grid_geometry(
    const Point3& global_lower,
    const Point3& global_upper,
    const GridDimensions& process_grid,
    const GridSpacing& requested_spacing);

void validate_atom_ghost_coverage(
    double atom_ghost_distance,
    double maximum_excluded_radius);

struct MpiDecompositionSpec {
    MPI_Comm communicator = MPI_COMM_NULL;
    Point3 global_lower{};
    Point3 global_upper{};
    Point3 local_lower{};
    Point3 local_upper{};
    GridDimensions process_grid{};
    VoxelCoord process_location{};
    double atom_ghost_distance = 0.0;
    double maximum_excluded_radius = 0.0;
};

struct OwnedVoxelRange {
    VoxelCoord begin{};
    VoxelCoord end{};
    GridDimensions dimensions{};
};

enum class Face : std::uint8_t {
    XLow = 0,
    XHigh = 1,
    YLow = 2,
    YHigh = 3,
    ZLow = 4,
    ZHigh = 5
};

// Records an owned gas voxel removed by solidification and its state
// immediately before the occupancy change.
struct DistributedRemovedVoxel {
    VoxelCoord voxel_coord{};
    GasState previous_state = GasState::Unclassified;
};

class MpiDecomposition {
public:
    MpiDecomposition(
        const GridSpec& global_grid_spec,
        MpiDecompositionSpec decomposition_spec);

    const MpiDecompositionSpec& spec() const noexcept;
    const OwnedVoxelRange& owned_range() const noexcept;
    int rank() const noexcept;
    int size() const noexcept;
    int neighbor_rank(Face face) const noexcept;
    bool owns(const VoxelCoord& global_voxel_coord) const noexcept;

private:
    static std::size_t face_index(Face face) noexcept;

    MpiDecompositionSpec spec_{};
    OwnedVoxelRange owned_range_{};
    int rank_ = 0;
    int size_ = 0;
    std::array<int, 6> neighbor_ranks_{};
};

class DistributedGasGrid {
public:
    DistributedGasGrid(
        GridSpec global_grid_spec,
        MpiDecompositionSpec decomposition_spec);

    DistributedGasGrid(const DistributedGasGrid&) = delete;
    DistributedGasGrid& operator=(const DistributedGasGrid&) = delete;

    const GridSpec& global_grid_spec() const noexcept;
    const MpiDecomposition& decomposition() const noexcept;
    const OwnedVoxelRange& owned_range() const noexcept;
    std::uint64_t owned_voxel_count() const noexcept;

    bool owns(const VoxelCoord& global_voxel_coord) const noexcept;
    std::optional<VoxelCoord> locate_voxel(const Point3& point) const noexcept;
    Point3 voxel_center(const VoxelCoord& global_voxel_coord) const;

    GasState gas_state(const VoxelCoord& global_voxel_coord) const;
    std::uint64_t owned_gas_state_count(GasState gas_state) const;
    void set_owned_gas_state(
        const VoxelCoord& global_voxel_coord,
        GasState gas_state);
    void fill_owned_gas_state(GasState gas_state);

    std::uint64_t voxelize_owned_atoms(
        AtomView atom_view,
        double precursor_radius);
    std::uint64_t voxelize_owned_atoms(
        AtomView atom_view,
        double precursor_radius,
        std::vector<DistributedRemovedVoxel>& removed_voxels);
    void exchange_ghost_states();
    bool is_site_accessible(const Point3& site_position) const;

private:
    struct LocalCoord {
        std::uint64_t x = 0;
        std::uint64_t y = 0;
        std::uint64_t z = 0;
    };

    static std::size_t face_index(Face face) noexcept;
    static Face opposite_face(Face face) noexcept;
    static std::size_t gas_state_index(GasState gas_state) noexcept;

    std::optional<LocalCoord> local_coord(
        const VoxelCoord& global_voxel_coord) const noexcept;
    std::optional<VoxelCoord> normalized_neighbor(
        const VoxelCoord& global_voxel_coord) const noexcept;
    std::size_t state_index(const LocalCoord& local_coord) const noexcept;
    std::size_t face_element_count(Face face) const;
    void pack_face(Face face);
    void unpack_face(Face face);
    void fill_ghost_face(Face face, GasState gas_state);
    std::uint64_t voxelize_owned_atoms_impl(
        AtomView atom_view,
        double precursor_radius,
        std::vector<DistributedRemovedVoxel>* removed_voxels);

    GridSpec global_grid_spec_{};
    MpiDecomposition decomposition_;
    GridDimensions storage_dimensions_{};
    std::uint64_t owned_voxel_count_ = 0;
    std::vector<GasState> states_{};
    std::array<std::uint64_t, 4> owned_state_counts_{};
    std::array<std::vector<GasState>, 6> send_buffers_{};
    std::array<std::vector<GasState>, 6> receive_buffers_{};
};

class DistributedGasAccessibilityQuery {
public:
    explicit DistributedGasAccessibilityQuery(
        const DistributedGasGrid& gas_grid) noexcept;

    bool is_site_accessible(const Point3& site_position) const;

private:
    const DistributedGasGrid& gas_grid_;
};

}  // namespace gasaccess

#endif
