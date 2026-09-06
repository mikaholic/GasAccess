#ifndef GASACCESS_ACCESSIBILITY_QUERY_HPP
#define GASACCESS_ACCESSIBILITY_QUERY_HPP

#include "gasaccess/gas_grid.hpp"

#include <cstddef>

namespace gasaccess {

inline constexpr std::size_t max_site_voxel_count = 27;

struct VoxelIdView {
    const VoxelId* voxel_ids = nullptr;
    std::size_t count = 0;
};

class GasAccessibilityQuery {
public:
    explicit GasAccessibilityQuery(const GasGrid& gas_grid) noexcept;

    bool is_voxel_outside_accessible(VoxelId voxel_id) const;

    // Tests the voxel containing the site and its six face neighbors. A site
    // outside any non-periodic grid boundary is inaccessible.
    bool is_site_accessible(const Point3& site_position) const;

    // Tests a caller-supplied local stencil without allocation. At most
    // max_site_voxel_count identifiers are accepted.
    bool is_site_accessible(VoxelIdView voxel_id_view) const;

private:
    const GasGrid& gas_grid_;
};

}  // namespace gasaccess

#endif
