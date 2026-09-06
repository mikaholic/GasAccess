#ifndef GASACCESS_TOOLS_EFFICIENCY_FIXTURE_HPP
#define GASACCESS_TOOLS_EFFICIENCY_FIXTURE_HPP

#include "gasaccess/gas_grid.hpp"

#include <cstdint>
#include <string>

namespace gasaccess::testing {

enum class EfficiencyChangeKind {
    Adsorption,
    Desorption,
    Mixed
};

enum class EfficiencyRepairCase {
    DetectionBaseline,
    Best,
    Medium,
    Worst
};

struct EfficiencyRepairFixture {
    GridSpec grid_spec{};
    GridDimensions process_grid{};
    VoxelCoord closing_cavity_begin{};
    VoxelCoord closing_cavity_end{};
    VoxelCoord opening_cavity_begin{};
    VoxelCoord opening_cavity_end{};
    VoxelCoord added_voxel{};
    VoxelCoord removed_voxel{};
    VoxelCoord closing_probe{};
    VoxelCoord opening_probe{};
    VoxelCoord outside_probe{};
    bool has_added_atom = false;
    bool has_removed_atom = false;
    VoxelBlockerCount initial_added_blocker_count = 0;
    VoxelBlockerCount initial_removed_blocker_count = 0;
    std::uint64_t expected_initial_closed_count = 0;
    std::uint64_t expected_final_closed_count = 0;
    std::uint64_t expected_newly_solid_count = 0;
    std::uint64_t expected_newly_gas_count = 0;
    std::uint64_t expected_newly_closed_count = 0;
    std::uint64_t expected_newly_opened_count = 0;
    std::uint64_t expected_changed_count = 0;
};

EfficiencyChangeKind parse_efficiency_change_kind(const std::string& value);
const char* efficiency_change_kind_name(EfficiencyChangeKind change_kind) noexcept;
EfficiencyRepairCase parse_efficiency_repair_case(const std::string& value);
const char* efficiency_repair_case_name(EfficiencyRepairCase repair_case) noexcept;

EfficiencyRepairFixture make_efficiency_repair_fixture(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case,
    const GridDimensions& dimensions,
    int process_count);

bool is_efficiency_fixture_solid(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case,
    const EfficiencyRepairFixture& fixture,
    const VoxelCoord& coordinate) noexcept;

}  // namespace gasaccess::testing

#endif
