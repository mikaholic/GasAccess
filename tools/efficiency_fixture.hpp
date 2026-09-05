#ifndef GASACCESS_TOOLS_EFFICIENCY_FIXTURE_HPP
#define GASACCESS_TOOLS_EFFICIENCY_FIXTURE_HPP

#include "gasaccess/gas_grid.hpp"

#include <cstdint>
#include <string>

namespace gasaccess::testing {

enum class EfficiencyRepairCase {
    DetectionBaseline,
    Best,
    Medium,
    Worst
};

struct EfficiencyRepairFixture {
    GridSpec grid_spec{};
    GridDimensions process_grid{};
    VoxelCoord cavity_begin{};
    VoxelCoord cavity_end{};
    VoxelCoord opening{};
    VoxelCoord cavity_probe{};
    VoxelCoord outside_probe{};
    std::uint64_t expected_initial_closed_count = 0;
    std::uint64_t expected_newly_closed_count = 0;
};

EfficiencyRepairCase parse_efficiency_repair_case(const std::string& value);
const char* efficiency_repair_case_name(EfficiencyRepairCase repair_case) noexcept;

EfficiencyRepairFixture make_efficiency_repair_fixture(
    EfficiencyRepairCase repair_case,
    const GridDimensions& dimensions,
    int process_count);

bool is_efficiency_fixture_solid(
    EfficiencyRepairCase repair_case,
    const EfficiencyRepairFixture& fixture,
    const VoxelCoord& coordinate) noexcept;

}  // namespace gasaccess::testing

#endif
