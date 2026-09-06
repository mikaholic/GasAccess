#ifndef GASACCESS_EXTERIOR_CLASSIFIER_HPP
#define GASACCESS_EXTERIOR_CLASSIFIER_HPP

#include "gasaccess/gas_grid.hpp"

namespace gasaccess {

struct ClassificationSummary {
    VoxelId solid_count = 0;
    VoxelId outside_accessible_count = 0;
    VoxelId closed_void_count = 0;
};

class ExteriorClassifier {
public:
    // Performs a full serial reference classification. GasState::Solid is
    // treated as occupancy; every other cached state is recomputed.
    ClassificationSummary classify(GasGrid& gas_grid) const;
};

}  // namespace gasaccess

#endif
