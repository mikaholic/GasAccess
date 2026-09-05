#include "efficiency_fixture.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace gasaccess::testing {
namespace {

std::int64_t as_int64(std::uint64_t value, const char* name)
{
    if (value > static_cast<std::uint64_t>(
                    std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error(std::string(name) + " exceeds int64");
    }
    return static_cast<std::int64_t>(value);
}

std::uint64_t checked_product(std::uint64_t lhs, std::uint64_t rhs)
{
    if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
        throw std::overflow_error("efficiency fixture volume overflows uint64");
    }
    return lhs * rhs;
}

std::uint64_t checked_volume(const GridDimensions& dimensions)
{
    return checked_product(
        checked_product(dimensions.x, dimensions.y),
        dimensions.z);
}

std::uint64_t region_volume(
    const VoxelCoord& begin,
    const VoxelCoord& end)
{
    if (end.x < begin.x || end.y < begin.y || end.z < begin.z) {
        throw std::logic_error("efficiency cavity bounds are invalid");
    }
    return checked_volume({
        static_cast<std::uint64_t>(end.x - begin.x),
        static_cast<std::uint64_t>(end.y - begin.y),
        static_cast<std::uint64_t>(end.z - begin.z)});
}

bool in_closed_range(std::int64_t value, std::int64_t low, std::int64_t high)
{
    return value >= low && value <= high;
}

bool on_cavity_shell(
    const VoxelCoord& begin,
    const VoxelCoord& end,
    const VoxelCoord& coordinate) noexcept
{
    const bool within_x = in_closed_range(
        coordinate.x, begin.x - 1, end.x);
    const bool within_y = in_closed_range(
        coordinate.y, begin.y - 1, end.y);
    const bool within_z = in_closed_range(
        coordinate.z, begin.z - 1, end.z);
    const bool on_x_shell = within_y && within_z
        && (coordinate.x == begin.x - 1 || coordinate.x == end.x);
    const bool on_y_shell = within_x && within_z
        && (coordinate.y == begin.y - 1 || coordinate.y == end.y);
    const bool on_z_shell = within_x && within_y
        && (coordinate.z == begin.z - 1 || coordinate.z == end.z);
    return on_x_shell || on_y_shell || on_z_shell;
}

void configure_single_cavity(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case,
    const GridDimensions& dimensions,
    std::uint64_t rank_count,
    EfficiencyRepairFixture& fixture)
{
    if (repair_case == EfficiencyRepairCase::Best) {
        const auto local_width = dimensions.x / rank_count;
        if (local_width < 12) {
            throw std::invalid_argument(
                "best-case fixture requires at least 12 x voxels per rank");
        }
        const auto cavity_x = std::min<std::uint64_t>(8, local_width - 4);
        const auto cavity_y = std::min<std::uint64_t>(8, dimensions.y - 4);
        const auto cavity_z = std::min<std::uint64_t>(8, dimensions.z - 3);
        fixture.closing_cavity_begin = {
            2,
            2,
            as_int64(dimensions.z - cavity_z - 2, "cavity z begin")};
        fixture.closing_cavity_end = {
            fixture.closing_cavity_begin.x
                + as_int64(cavity_x, "cavity x"),
            fixture.closing_cavity_begin.y
                + as_int64(cavity_y, "cavity y"),
            fixture.closing_cavity_begin.z
                + as_int64(cavity_z, "cavity z")};
    } else {
        const auto cavity_depth = repair_case == EfficiencyRepairCase::Worst
            ? (dimensions.z * 3U) / 4U
            : dimensions.z / 4U;
        fixture.closing_cavity_begin = {0, 0, 0};
        fixture.closing_cavity_end = {
            as_int64(dimensions.x, "x dimension"),
            as_int64(dimensions.y, "y dimension"),
            as_int64(cavity_depth, "cavity depth")};
    }
    fixture.opening_cavity_begin = fixture.closing_cavity_begin;
    fixture.opening_cavity_end = fixture.closing_cavity_end;
    const auto cavity_volume = region_volume(
        fixture.closing_cavity_begin,
        fixture.closing_cavity_end);
    const auto cavity_center = VoxelCoord{
        (fixture.closing_cavity_begin.x + fixture.closing_cavity_end.x) / 2,
        (fixture.closing_cavity_begin.y + fixture.closing_cavity_end.y) / 2,
        fixture.closing_cavity_end.z};
    fixture.added_voxel = cavity_center;
    fixture.removed_voxel = cavity_center;
    fixture.closing_probe = fixture.closing_cavity_begin;
    fixture.opening_probe = fixture.opening_cavity_begin;

    const bool baseline =
        repair_case == EfficiencyRepairCase::DetectionBaseline;
    if (change_kind == EfficiencyChangeKind::Adsorption) {
        fixture.has_added_atom = true;
        fixture.initial_added_blocker_count = baseline ? 1U : 0U;
        fixture.expected_initial_closed_count = baseline ? cavity_volume : 0U;
        fixture.expected_final_closed_count = cavity_volume;
        if (!baseline) {
            fixture.expected_newly_solid_count = 1;
            fixture.expected_newly_closed_count = cavity_volume;
            fixture.expected_changed_count = cavity_volume + 1U;
        }
        return;
    }

    fixture.has_removed_atom = true;
    fixture.initial_removed_blocker_count = baseline ? 2U : 1U;
    fixture.expected_initial_closed_count = cavity_volume;
    fixture.expected_final_closed_count = baseline ? cavity_volume : 0U;
    if (!baseline) {
        fixture.expected_newly_gas_count = 1;
        fixture.expected_newly_opened_count = cavity_volume + 1U;
        fixture.expected_changed_count = cavity_volume + 1U;
    }
}

void configure_mixed_fixture(
    EfficiencyRepairCase repair_case,
    const GridDimensions& dimensions,
    std::uint64_t rank_count,
    EfficiencyRepairFixture& fixture)
{
    fixture.has_added_atom = true;
    fixture.has_removed_atom = true;
    if (repair_case == EfficiencyRepairCase::DetectionBaseline) {
        fixture.added_voxel = {
            as_int64(dimensions.x / 2U, "baseline x"),
            as_int64(dimensions.y / 2U, "baseline y"),
            as_int64(dimensions.z / 2U, "baseline z")};
        fixture.removed_voxel = fixture.added_voxel;
        fixture.initial_added_blocker_count = 1;
        fixture.initial_removed_blocker_count = 1;
        fixture.closing_probe = {0, 0, 0};
        fixture.opening_probe = fixture.closing_probe;
        return;
    }

    if (repair_case == EfficiencyRepairCase::Best) {
        const auto local_width = dimensions.x / rank_count;
        if (local_width < 12 || dimensions.y < 23) {
            throw std::invalid_argument(
                "mixed best case requires 12 local x and 23 y voxels");
        }
        const std::int64_t cavity_extent = 8;
        const auto z_begin = as_int64(dimensions.z, "z dimension") - 10;
        fixture.closing_cavity_begin = {2, 2, z_begin};
        fixture.closing_cavity_end = {
            2 + cavity_extent,
            2 + cavity_extent,
            z_begin + cavity_extent};
        fixture.opening_cavity_begin = {
            2,
            as_int64(dimensions.y, "y dimension") - 10,
            z_begin};
        fixture.opening_cavity_end = {
            2 + cavity_extent,
            fixture.opening_cavity_begin.y + cavity_extent,
            z_begin + cavity_extent};
    } else {
        const auto cavity_depth = repair_case == EfficiencyRepairCase::Worst
            ? (dimensions.z * 3U) / 4U
            : dimensions.z / 4U;
        const auto separator_y = dimensions.y / 2U;
        fixture.closing_cavity_begin = {0, 0, 0};
        fixture.closing_cavity_end = {
            as_int64(dimensions.x, "x dimension"),
            as_int64(separator_y, "separator y"),
            as_int64(cavity_depth, "cavity depth")};
        fixture.opening_cavity_begin = {
            0,
            as_int64(separator_y + 1U, "opening cavity y"),
            0};
        fixture.opening_cavity_end = {
            as_int64(dimensions.x, "x dimension"),
            as_int64(dimensions.y, "y dimension"),
            as_int64(cavity_depth, "cavity depth")};
    }

    fixture.added_voxel = {
        (fixture.closing_cavity_begin.x + fixture.closing_cavity_end.x) / 2,
        (fixture.closing_cavity_begin.y + fixture.closing_cavity_end.y) / 2,
        fixture.closing_cavity_end.z};
    fixture.removed_voxel = {
        (fixture.opening_cavity_begin.x + fixture.opening_cavity_end.x) / 2,
        (fixture.opening_cavity_begin.y + fixture.opening_cavity_end.y) / 2,
        fixture.opening_cavity_end.z};
    fixture.closing_probe = fixture.closing_cavity_begin;
    fixture.opening_probe = fixture.opening_cavity_begin;
    fixture.initial_removed_blocker_count = 1;
    const auto closing_volume = region_volume(
        fixture.closing_cavity_begin,
        fixture.closing_cavity_end);
    const auto opening_volume = region_volume(
        fixture.opening_cavity_begin,
        fixture.opening_cavity_end);
    fixture.expected_initial_closed_count = opening_volume;
    fixture.expected_final_closed_count = closing_volume;
    fixture.expected_newly_solid_count = 1;
    fixture.expected_newly_gas_count = 1;
    fixture.expected_newly_closed_count = closing_volume;
    fixture.expected_newly_opened_count = opening_volume + 1U;
    fixture.expected_changed_count = closing_volume + opening_volume + 2U;
}

}  // namespace

EfficiencyChangeKind parse_efficiency_change_kind(const std::string& value)
{
    if (value == "adsorption") {
        return EfficiencyChangeKind::Adsorption;
    }
    if (value == "desorption") {
        return EfficiencyChangeKind::Desorption;
    }
    if (value == "mixed") {
        return EfficiencyChangeKind::Mixed;
    }
    throw std::invalid_argument(
        "change kind must be adsorption, desorption, or mixed");
}

const char* efficiency_change_kind_name(
    EfficiencyChangeKind change_kind) noexcept
{
    switch (change_kind) {
    case EfficiencyChangeKind::Adsorption:
        return "adsorption";
    case EfficiencyChangeKind::Desorption:
        return "desorption";
    case EfficiencyChangeKind::Mixed:
        return "mixed";
    }
    return "unknown";
}

EfficiencyRepairCase parse_efficiency_repair_case(const std::string& value)
{
    if (value == "baseline") {
        return EfficiencyRepairCase::DetectionBaseline;
    }
    if (value == "best") {
        return EfficiencyRepairCase::Best;
    }
    if (value == "medium") {
        return EfficiencyRepairCase::Medium;
    }
    if (value == "worst") {
        return EfficiencyRepairCase::Worst;
    }
    throw std::invalid_argument(
        "repair case must be baseline, best, medium, or worst");
}

const char* efficiency_repair_case_name(
    EfficiencyRepairCase repair_case) noexcept
{
    switch (repair_case) {
    case EfficiencyRepairCase::DetectionBaseline:
        return "baseline";
    case EfficiencyRepairCase::Best:
        return "best";
    case EfficiencyRepairCase::Medium:
        return "medium";
    case EfficiencyRepairCase::Worst:
        return "worst";
    }
    return "unknown";
}

EfficiencyRepairFixture make_efficiency_repair_fixture(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case,
    const GridDimensions& dimensions,
    int process_count)
{
    if (process_count <= 0) {
        throw std::invalid_argument("MPI process count must be positive");
    }
    if (dimensions.x < 16 || dimensions.y < 12 || dimensions.z < 12) {
        throw std::invalid_argument(
            "efficiency repair dimensions must be at least 16x12x12");
    }
    const auto rank_count = static_cast<std::uint64_t>(process_count);
    if (dimensions.x % rank_count != 0) {
        throw std::invalid_argument(
            "efficiency fixture x dimension must be divisible by rank count");
    }

    EfficiencyRepairFixture fixture{};
    fixture.grid_spec.spacing = {1.0, 1.0, 1.0};
    fixture.grid_spec.dimensions = dimensions;
    fixture.grid_spec.reservoir_faces.z_high = true;
    fixture.process_grid = {rank_count, 1, 1};
    fixture.outside_probe = {
        0,
        0,
        as_int64(dimensions.z, "z dimension") - 1};

    if (change_kind == EfficiencyChangeKind::Mixed) {
        configure_mixed_fixture(
            repair_case,
            dimensions,
            rank_count,
            fixture);
    } else {
        configure_single_cavity(
            change_kind,
            repair_case,
            dimensions,
            rank_count,
            fixture);
    }
    return fixture;
}

bool is_efficiency_fixture_solid(
    EfficiencyChangeKind change_kind,
    EfficiencyRepairCase repair_case,
    const EfficiencyRepairFixture& fixture,
    const VoxelCoord& coordinate) noexcept
{
    if (change_kind == EfficiencyChangeKind::Mixed) {
        if (repair_case == EfficiencyRepairCase::DetectionBaseline) {
            return coordinate == fixture.added_voxel;
        }
        if (repair_case == EfficiencyRepairCase::Best) {
            const bool closing_shell = on_cavity_shell(
                fixture.closing_cavity_begin,
                fixture.closing_cavity_end,
                coordinate);
            const bool opening_shell = on_cavity_shell(
                fixture.opening_cavity_begin,
                fixture.opening_cavity_end,
                coordinate);
            return (closing_shell && coordinate != fixture.added_voxel)
                || opening_shell;
        }
        const auto separator_y = fixture.closing_cavity_end.y;
        const auto roof_z = fixture.closing_cavity_end.z;
        if (coordinate.z == roof_z) {
            return coordinate != fixture.added_voxel;
        }
        return coordinate.y == separator_y && coordinate.z < roof_z;
    }

    const bool plug_present =
        change_kind == EfficiencyChangeKind::Desorption
        || repair_case == EfficiencyRepairCase::DetectionBaseline;
    if (repair_case != EfficiencyRepairCase::Best) {
        return coordinate.z == fixture.closing_cavity_end.z
            && (plug_present || coordinate != fixture.added_voxel);
    }
    return on_cavity_shell(
               fixture.closing_cavity_begin,
               fixture.closing_cavity_end,
               coordinate)
        && (plug_present || coordinate != fixture.added_voxel);
}

}  // namespace gasaccess::testing
