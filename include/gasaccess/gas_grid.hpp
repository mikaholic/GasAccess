#ifndef GASACCESS_GAS_GRID_HPP
#define GASACCESS_GAS_GRID_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace gasaccess {

using VoxelId = std::uint64_t;

struct VoxelCoord {
    std::int64_t x = 0;
    std::int64_t y = 0;
    std::int64_t z = 0;
};

constexpr bool operator==(const VoxelCoord& lhs, const VoxelCoord& rhs) noexcept
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

constexpr bool operator!=(const VoxelCoord& lhs, const VoxelCoord& rhs) noexcept
{
    return !(lhs == rhs);
}

struct Point3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct GridSpacing {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct GridDimensions {
    std::uint64_t x = 0;
    std::uint64_t y = 0;
    std::uint64_t z = 0;
};

struct PeriodicAxes {
    bool x = false;
    bool y = false;
    bool z = false;
};

struct ReservoirFaces {
    bool x_low = false;
    bool x_high = false;
    bool y_low = false;
    bool y_high = false;
    bool z_low = false;
    bool z_high = false;
};

enum class GasState : std::uint8_t {
    Unclassified = 0,
    Solid = 1,
    OutsideAccessible = 2,
    ClosedVoid = 3
};

static_assert(sizeof(GasState) == 1, "GasState must occupy one byte");
static_assert(sizeof(VoxelId) == 8, "VoxelId must occupy eight bytes");

struct GridSpec {
    Point3 origin{};
    GridSpacing spacing{};
    GridDimensions dimensions{};
    PeriodicAxes periodic{};
    ReservoirFaces reservoir_faces{};
    std::vector<VoxelCoord> explicit_source_voxels{};
};

struct NeighborList {
    std::array<VoxelId, 6> ids{};
    std::size_t count = 0;
};

class GasGrid {
public:
    explicit GasGrid(GridSpec grid_spec);

    const GridSpec& grid_spec() const noexcept;
    VoxelId voxel_count() const noexcept;

    bool contains(const VoxelCoord& voxel_coord) const noexcept;
    VoxelId voxel_id(const VoxelCoord& voxel_coord) const;
    VoxelCoord voxel_coord(VoxelId voxel_id) const;

    std::optional<VoxelCoord> locate_voxel(const Point3& point) const noexcept;
    Point3 voxel_center(const VoxelCoord& voxel_coord) const;
    Point3 voxel_center(VoxelId voxel_id) const;

    NeighborList neighbors(VoxelId voxel_id) const;

    bool is_reservoir_source(const VoxelCoord& voxel_coord) const;
    bool is_reservoir_source(VoxelId voxel_id) const;
    const std::vector<VoxelId>& explicit_source_ids() const noexcept;

    GasState gas_state(VoxelId voxel_id) const;
    GasState gas_state(const VoxelCoord& voxel_coord) const;
    VoxelId gas_state_count(GasState gas_state) const;
    void set_gas_state(VoxelId voxel_id, GasState gas_state);
    void fill_gas_state(GasState gas_state);

private:
    static std::size_t state_index(GasState gas_state);
    static VoxelId validate_and_count_voxels(const GridSpec& grid_spec);
    static void validate_boundary_conditions(const GridSpec& grid_spec);

    std::optional<std::int64_t> locate_axis(
        double position,
        double origin,
        double spacing,
        std::uint64_t dimension,
        bool periodic) const noexcept;
    bool is_boundary_source(const VoxelCoord& voxel_coord) const noexcept;
    void initialize_explicit_sources();
    void add_unique_neighbor(NeighborList& neighbor_list, const VoxelCoord& voxel_coord) const;

    GridSpec grid_spec_;
    VoxelId voxel_count_ = 0;
    std::vector<GasState> states_;
    std::array<VoxelId, 4> state_counts_{};
    std::vector<VoxelId> explicit_source_ids_;
};

}  // namespace gasaccess

#endif
