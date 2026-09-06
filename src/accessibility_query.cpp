#include "gasaccess/accessibility_query.hpp"

#include <stdexcept>

namespace gasaccess {

GasAccessibilityQuery::GasAccessibilityQuery(const GasGrid& gas_grid) noexcept
    : gas_grid_(gas_grid)
{
}

bool GasAccessibilityQuery::is_voxel_outside_accessible(VoxelId voxel_id) const
{
    return gas_grid_.gas_state(voxel_id) == GasState::OutsideAccessible;
}

bool GasAccessibilityQuery::is_site_accessible(const Point3& site_position) const
{
    const auto site_voxel_coord = gas_grid_.locate_voxel(site_position);
    if (!site_voxel_coord) {
        return false;
    }

    const auto site_voxel_id = gas_grid_.voxel_id(*site_voxel_coord);
    if (is_voxel_outside_accessible(site_voxel_id)) {
        return true;
    }

    const auto neighbor_list = gas_grid_.neighbors(site_voxel_id);
    for (std::size_t index = 0; index < neighbor_list.count; ++index) {
        if (is_voxel_outside_accessible(neighbor_list.ids[index])) {
            return true;
        }
    }
    return false;
}

bool GasAccessibilityQuery::is_site_accessible(VoxelIdView voxel_id_view) const
{
    if (voxel_id_view.count > max_site_voxel_count) {
        throw std::invalid_argument("site voxel stencil exceeds the maximum size");
    }
    if (voxel_id_view.count != 0 && voxel_id_view.voxel_ids == nullptr) {
        throw std::invalid_argument("site voxel stencil has a null pointer");
    }

    for (std::size_t index = 0; index < voxel_id_view.count; ++index) {
        if (is_voxel_outside_accessible(voxel_id_view.voxel_ids[index])) {
            return true;
        }
    }
    return false;
}

}  // namespace gasaccess
