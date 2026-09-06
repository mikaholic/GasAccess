#ifndef GASACCESS_TOOLS_MOCK_SPPARKS_FIXTURE_HPP
#define GASACCESS_TOOLS_MOCK_SPPARKS_FIXTURE_HPP

#include "mock_spparks.hpp"

#include <string>
#include <vector>

namespace gasaccess::testing {

enum class MockScenario {
    OpenTrench,
    SealedTrench,
    MillionSlab
};

struct MockFixtureSpec {
    GridDimensions dimensions{};
    Point3 global_lower{};
    Point3 global_upper{};
    GridSpacing requested_spacing{};
    PeriodicAxes spparks_periodic{};
    PeriodicAxes gas_periodic{};
    ReservoirFaces reservoir_faces{};
    double atom_radius = 0.0;
    double precursor_radius = 0.0;
    double ghost_distance = 0.0;
};

MockScenario parse_mock_scenario(const std::string& value);
const char* mock_scenario_name(MockScenario scenario) noexcept;
MockFixtureSpec make_mock_fixture_spec(MockScenario scenario);
bool is_mock_solid(MockScenario scenario, const VoxelCoord& coordinate) noexcept;
bool is_inner_trench_probe(
    MockScenario scenario,
    const VoxelCoord& coordinate) noexcept;
std::vector<MockSpparksAtom> make_mock_atoms(
    MockScenario scenario,
    const MockFixtureSpec& fixture_spec,
    const MockSpparksDomain* owner = nullptr);
GridDimensions make_mock_process_grid(int process_count);

}  // namespace gasaccess::testing

#endif
