#include "mock_spparks_fixture.hpp"

#include <limits>
#include <stdexcept>

namespace gasaccess::testing {
namespace {

Point3 voxel_center(
    const MockFixtureSpec& fixture_spec,
    const VoxelCoord& coordinate) noexcept
{
    return {
        fixture_spec.global_lower.x
            + (static_cast<double>(coordinate.x) + 0.5)
                * fixture_spec.requested_spacing.x,
        fixture_spec.global_lower.y
            + (static_cast<double>(coordinate.y) + 0.5)
                * fixture_spec.requested_spacing.y,
        fixture_spec.global_lower.z
            + (static_cast<double>(coordinate.z) + 0.5)
                * fixture_spec.requested_spacing.z};
}

std::int64_t atom_id(
    const GridDimensions& dimensions,
    const VoxelCoord& coordinate)
{
    const auto id = static_cast<std::uint64_t>(coordinate.x)
        + dimensions.x
            * (static_cast<std::uint64_t>(coordinate.y)
               + dimensions.y * static_cast<std::uint64_t>(coordinate.z));
    if (id >= static_cast<std::uint64_t>(
                  std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error("mock SPPARKS atom ID exceeds int64");
    }
    return static_cast<std::int64_t>(id + 1);
}

}  // namespace

MockScenario parse_mock_scenario(const std::string& value)
{
    if (value == "open-trench") {
        return MockScenario::OpenTrench;
    }
    if (value == "sealed-trench") {
        return MockScenario::SealedTrench;
    }
    if (value == "million-slab") {
        return MockScenario::MillionSlab;
    }
    throw std::invalid_argument(
        "scenario must be open-trench, sealed-trench, or million-slab");
}

const char* mock_scenario_name(MockScenario scenario) noexcept
{
    switch (scenario) {
    case MockScenario::OpenTrench:
        return "open-trench";
    case MockScenario::SealedTrench:
        return "sealed-trench";
    case MockScenario::MillionSlab:
        return "million-slab";
    }
    return "unknown";
}

MockFixtureSpec make_mock_fixture_spec(MockScenario scenario)
{
    MockFixtureSpec fixture_spec{};
    fixture_spec.dimensions = scenario == MockScenario::MillionSlab
        ? GridDimensions{128, 128, 128}
        : GridDimensions{16, 16, 16};
    fixture_spec.global_lower = {0.0, 0.0, 0.0};
    fixture_spec.global_upper = {
        static_cast<double>(fixture_spec.dimensions.x),
        static_cast<double>(fixture_spec.dimensions.y),
        static_cast<double>(fixture_spec.dimensions.z)};
    fixture_spec.requested_spacing = {1.0, 1.0, 1.0};

    // This intentionally models the production distinction discussed with the
    // user: SPPARKS may be periodic in all axes while GasAccess opens z at the
    // top reservoir.
    fixture_spec.spparks_periodic = {true, true, true};
    fixture_spec.gas_periodic = {true, true, false};
    fixture_spec.reservoir_faces.z_high = true;
    fixture_spec.atom_radius = 0.30;
    fixture_spec.precursor_radius = 0.20;
    fixture_spec.ghost_distance = 1.0;
    return fixture_spec;
}

bool is_mock_solid(
    MockScenario scenario,
    const VoxelCoord& coordinate) noexcept
{
    if (scenario == MockScenario::MillionSlab) {
        return coordinate.z < 64;
    }

    const bool substrate = coordinate.z <= 2;
    const bool sidewalls = coordinate.z >= 3 && coordinate.z <= 11
        && ((coordinate.x >= 2 && coordinate.x <= 4)
            || (coordinate.x >= 11 && coordinate.x <= 13));
    const bool cap = scenario == MockScenario::SealedTrench
        && coordinate.z >= 9 && coordinate.z <= 11
        && coordinate.x >= 5 && coordinate.x <= 10;
    return substrate || sidewalls || cap;
}

bool is_inner_trench_probe(
    MockScenario scenario,
    const VoxelCoord& coordinate) noexcept
{
    if (scenario == MockScenario::MillionSlab) {
        return false;
    }
    return coordinate.z == 5
        && (coordinate.x == 4 || coordinate.x == 11);
}

std::vector<MockSpparksAtom> make_mock_atoms(
    MockScenario scenario,
    const MockFixtureSpec& fixture_spec,
    const MockSpparksDomain* owner)
{
    std::vector<MockSpparksAtom> atoms;
    if (scenario == MockScenario::MillionSlab && owner == nullptr) {
        atoms.reserve(128U * 128U * 64U);
    }
    for (std::uint64_t z = 0; z < fixture_spec.dimensions.z; ++z) {
        for (std::uint64_t y = 0; y < fixture_spec.dimensions.y; ++y) {
            for (std::uint64_t x = 0; x < fixture_spec.dimensions.x; ++x) {
                const VoxelCoord coordinate{
                    static_cast<std::int64_t>(x),
                    static_cast<std::int64_t>(y),
                    static_cast<std::int64_t>(z)};
                if (!is_mock_solid(scenario, coordinate)) {
                    continue;
                }
                const auto position = voxel_center(fixture_spec, coordinate);
                if (owner != nullptr && !owner->owns(position)) {
                    continue;
                }
                atoms.push_back({
                    atom_id(fixture_spec.dimensions, coordinate),
                    position,
                    fixture_spec.atom_radius});
            }
        }
    }
    return atoms;
}

GridDimensions make_mock_process_grid(int process_count)
{
    if (process_count <= 0) {
        throw std::invalid_argument("MPI process count must be positive");
    }
    int dimensions[3]{0, 0, 0};
    if (MPI_Dims_create(process_count, 3, dimensions) != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Dims_create failed");
    }
    return {
        static_cast<std::uint64_t>(dimensions[0]),
        static_cast<std::uint64_t>(dimensions[1]),
        static_cast<std::uint64_t>(dimensions[2])};
}

}  // namespace gasaccess::testing
