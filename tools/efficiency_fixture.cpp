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

std::uint64_t checked_volume(const GridDimensions& dimensions)
{
    if (dimensions.x != 0
        && dimensions.y > std::numeric_limits<std::uint64_t>::max()
                / dimensions.x) {
        throw std::overflow_error("efficiency fixture volume overflows uint64");
    }
    const auto xy = dimensions.x * dimensions.y;
    if (dimensions.z != 0
        && xy > std::numeric_limits<std::uint64_t>::max() / dimensions.z) {
        throw std::overflow_error("efficiency fixture volume overflows uint64");
    }
    return xy * dimensions.z;
}

bool in_closed_range(std::int64_t value, std::int64_t low, std::int64_t high)
{
    return value >= low && value <= high;
}

}  // namespace

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

    if (repair_case == EfficiencyRepairCase::Best) {
        const auto local_width = dimensions.x / rank_count;
        if (local_width < 12) {
            throw std::invalid_argument(
                "best-case fixture requires at least 12 x voxels per rank");
        }
        const auto cavity_x = std::min<std::uint64_t>(8, local_width - 4);
        const auto cavity_y = std::min<std::uint64_t>(8, dimensions.y - 4);
        const auto cavity_z = std::min<std::uint64_t>(8, dimensions.z - 3);
        fixture.cavity_begin = {
            2,
            2,
            as_int64(dimensions.z - cavity_z - 2, "cavity z begin")};
        fixture.cavity_end = {
            fixture.cavity_begin.x + as_int64(cavity_x, "cavity x"),
            fixture.cavity_begin.y + as_int64(cavity_y, "cavity y"),
            fixture.cavity_begin.z + as_int64(cavity_z, "cavity z")};
        fixture.opening = {
            fixture.cavity_begin.x + as_int64(cavity_x / 2, "opening x"),
            fixture.cavity_begin.y + as_int64(cavity_y / 2, "opening y"),
            fixture.cavity_end.z};
        fixture.cavity_probe = fixture.cavity_begin;
        fixture.expected_newly_closed_count = checked_volume(
            {cavity_x, cavity_y, cavity_z});
        return fixture;
    }

    const auto cavity_depth = repair_case == EfficiencyRepairCase::Worst
        ? (dimensions.z * 3U) / 4U
        : dimensions.z / 4U;
    fixture.cavity_begin = {0, 0, 0};
    fixture.cavity_end = {
        as_int64(dimensions.x, "x dimension"),
        as_int64(dimensions.y, "y dimension"),
        as_int64(cavity_depth, "cavity depth")};
    fixture.opening = {
        as_int64(dimensions.x / 2U, "opening x"),
        as_int64(dimensions.y / 2U, "opening y"),
        as_int64(cavity_depth, "opening z")};
    fixture.cavity_probe = {0, 0, 0};
    const auto cavity_volume = checked_volume(
        {dimensions.x, dimensions.y, cavity_depth});
    if (repair_case == EfficiencyRepairCase::DetectionBaseline) {
        fixture.expected_initial_closed_count = cavity_volume;
    } else {
        fixture.expected_newly_closed_count = cavity_volume;
    }
    return fixture;
}

bool is_efficiency_fixture_solid(
    EfficiencyRepairCase repair_case,
    const EfficiencyRepairFixture& fixture,
    const VoxelCoord& coordinate) noexcept
{
    if (repair_case != EfficiencyRepairCase::Best) {
        if (coordinate.z != fixture.opening.z) {
            return false;
        }
        return repair_case == EfficiencyRepairCase::DetectionBaseline
            || coordinate != fixture.opening;
    }

    const auto& begin = fixture.cavity_begin;
    const auto& end = fixture.cavity_end;
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
    return (on_x_shell || on_y_shell || on_z_shell)
        && coordinate != fixture.opening;
}

}  // namespace gasaccess::testing
