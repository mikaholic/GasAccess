#include "gasaccess/gas_grid.hpp"

#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridSpec;
using gasaccess::NeighborList;
using gasaccess::Point3;
using gasaccess::VoxelCoord;
using gasaccess::VoxelId;

class TestFailure : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

void require_condition(
    bool condition,
    const char* expression,
    const char* file,
    int line)
{
    if (condition) {
        return;
    }

    std::ostringstream message;
    message << file << ':' << line << ": check failed: " << expression;
    throw TestFailure(message.str());
}

void require_near(
    double actual,
    double expected,
    double tolerance,
    const char* file,
    int line)
{
    if (std::abs(actual - expected) <= tolerance) {
        return;
    }

    std::ostringstream message;
    message << file << ':' << line << ": expected " << expected
            << ", received " << actual;
    throw TestFailure(message.str());
}

template <typename Exception, typename Function>
void require_throws(Function&& function, const char* expression, const char* file, int line)
{
    try {
        std::forward<Function>(function)();
    } catch (const Exception&) {
        return;
    } catch (const std::exception& error) {
        std::ostringstream message;
        message << file << ':' << line << ": " << expression
                << " threw a different exception: " << error.what();
        throw TestFailure(message.str());
    }

    std::ostringstream message;
    message << file << ':' << line << ": " << expression << " did not throw";
    throw TestFailure(message.str());
}

#define REQUIRE(condition) \
    require_condition((condition), #condition, __FILE__, __LINE__)
#define REQUIRE_NEAR(actual, expected, tolerance) \
    require_near((actual), (expected), (tolerance), __FILE__, __LINE__)
#define REQUIRE_THROWS_AS(expression, exception_type) \
    require_throws<exception_type>([&]() { expression; }, #expression, __FILE__, __LINE__)

GridSpec make_grid_spec(
    std::uint64_t x = 3,
    std::uint64_t y = 3,
    std::uint64_t z = 3)
{
    GridSpec grid_spec{};
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {x, y, z};
    return grid_spec;
}

std::set<VoxelId> neighbor_ids(const NeighborList& neighbor_list)
{
    std::set<VoxelId> result;
    for (std::size_t index = 0; index < neighbor_list.count; ++index) {
        result.insert(neighbor_list.ids[index]);
    }
    return result;
}

std::set<VoxelId> ids_for(
    const GasGrid& grid,
    const std::vector<VoxelCoord>& voxel_coords)
{
    std::set<VoxelId> result;
    for (const auto& voxel_coord : voxel_coords) {
        result.insert(grid.voxel_id(voxel_coord));
    }
    return result;
}

void test_grid_indexing_and_state_storage()
{
    auto grid_spec = make_grid_spec(4, 3, 2);
    grid_spec.origin = {-1.0, 2.0, 10.0};
    grid_spec.spacing = {0.5, 0.5, 0.5};
    GasGrid grid(grid_spec);

    REQUIRE(grid.voxel_count() == 24);
    REQUIRE(grid.gas_state_count(GasState::Unclassified) == 24);
    REQUIRE(grid.gas_state_count(GasState::Solid) == 0);
    REQUIRE(grid.gas_state_count(GasState::OutsideAccessible) == 0);
    REQUIRE(grid.gas_state_count(GasState::ClosedVoid) == 0);

    for (std::int64_t z = 0; z < 2; ++z) {
        for (std::int64_t y = 0; y < 3; ++y) {
            for (std::int64_t x = 0; x < 4; ++x) {
                const VoxelCoord voxel_coord{x, y, z};
                const auto expected = static_cast<VoxelId>((z * 3 + y) * 4 + x);
                const auto id = grid.voxel_id(voxel_coord);
                REQUIRE(id == expected);
                REQUIRE(grid.voxel_coord(id) == voxel_coord);
                REQUIRE(grid.gas_state(id) == GasState::Unclassified);
            }
        }
    }

    REQUIRE(!grid.contains({-1, 0, 0}));
    REQUIRE(!grid.contains({4, 0, 0}));
    REQUIRE_THROWS_AS(grid.voxel_id({4, 0, 0}), std::out_of_range);
    REQUIRE_THROWS_AS(grid.voxel_coord(grid.voxel_count()), std::out_of_range);

    const auto center = grid.voxel_center({2, 1, 0});
    REQUIRE_NEAR(center.x, 0.25, 1.0e-12);
    REQUIRE_NEAR(center.y, 2.75, 1.0e-12);
    REQUIRE_NEAR(center.z, 10.25, 1.0e-12);

    const auto selected_id = grid.voxel_id({2, 1, 0});
    grid.set_gas_state(selected_id, GasState::Solid);
    REQUIRE(grid.gas_state(selected_id) == GasState::Solid);
    REQUIRE(grid.gas_state_count(GasState::Unclassified) == 23);
    REQUIRE(grid.gas_state_count(GasState::Solid) == 1);
    grid.set_gas_state(selected_id, GasState::Solid);
    REQUIRE(grid.gas_state_count(GasState::Unclassified) == 23);
    REQUIRE(grid.gas_state_count(GasState::Solid) == 1);

    grid.fill_gas_state(GasState::ClosedVoid);
    for (VoxelId id = 0; id < grid.voxel_count(); ++id) {
        REQUIRE(grid.gas_state(id) == GasState::ClosedVoid);
    }
    REQUIRE(grid.gas_state_count(GasState::Unclassified) == 0);
    REQUIRE(grid.gas_state_count(GasState::Solid) == 0);
    REQUIRE(grid.gas_state_count(GasState::OutsideAccessible) == 0);
    REQUIRE(grid.gas_state_count(GasState::ClosedVoid) == 24);

    REQUIRE_THROWS_AS(
        grid.set_gas_state(grid.voxel_count(), GasState::Solid),
        std::out_of_range);
    const auto invalid_state = static_cast<GasState>(255);
    REQUIRE_THROWS_AS(
        grid.set_gas_state(selected_id, invalid_state),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        grid.fill_gas_state(invalid_state),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        grid.gas_state_count(invalid_state),
        std::invalid_argument);
    REQUIRE(grid.gas_state_count(GasState::ClosedVoid) == 24);
}

void test_blocker_count_state_invariant()
{
    GasGrid grid(make_grid_spec(3, 2, 1));
    const auto id = grid.voxel_id({1, 1, 0});
    for (VoxelId voxel_id = 0; voxel_id < grid.voxel_count(); ++voxel_id) {
        REQUIRE(grid.blocker_count(voxel_id) == 0);
    }

    grid.set_gas_state(id, GasState::Solid);
    REQUIRE(grid.blocker_count(id) == 1);
    grid.set_gas_state(id, GasState::OutsideAccessible);
    REQUIRE(grid.blocker_count(id) == 0);

    grid.set_blocker_count(id, 2);
    REQUIRE(grid.blocker_count(id) == 2);
    REQUIRE(grid.gas_state(id) == GasState::Solid);
    REQUIRE(grid.gas_state_count(GasState::Solid) == 1);
    grid.set_blocker_count(id, 1);
    REQUIRE(grid.gas_state_count(GasState::Solid) == 1);
    grid.set_blocker_count(id, 0);
    REQUIRE(grid.gas_state(id) == GasState::Unclassified);
    REQUIRE(grid.gas_state_count(GasState::Solid) == 0);

    grid.set_gas_state(id, GasState::ClosedVoid);
    grid.set_blocker_count(id, 0);
    REQUIRE(grid.gas_state(id) == GasState::ClosedVoid);

    grid.fill_gas_state(GasState::Solid);
    for (VoxelId voxel_id = 0; voxel_id < grid.voxel_count(); ++voxel_id) {
        REQUIRE(grid.blocker_count(voxel_id) == 1);
    }
    grid.fill_gas_state(GasState::OutsideAccessible);
    for (VoxelId voxel_id = 0; voxel_id < grid.voxel_count(); ++voxel_id) {
        REQUIRE(grid.blocker_count(voxel_id) == 0);
    }

    REQUIRE_THROWS_AS(grid.blocker_count(grid.voxel_count()), std::out_of_range);
    REQUIRE_THROWS_AS(
        grid.set_blocker_count(grid.voxel_count(), 1),
        std::out_of_range);
}

void test_invalid_grid_specs()
{
    auto grid_spec = make_grid_spec();
    grid_spec.spacing.x = 0.0;
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::invalid_argument);

    grid_spec = make_grid_spec();
    grid_spec.spacing.y = -1.0;
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::invalid_argument);

    grid_spec = make_grid_spec();
    grid_spec.spacing.z = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::invalid_argument);

    grid_spec = make_grid_spec();
    grid_spec.origin.x = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::invalid_argument);

    grid_spec = make_grid_spec();
    grid_spec.origin.x = 1.0e16;
    grid_spec.spacing = {0.1, 0.1, 0.1};
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::invalid_argument);

    grid_spec = make_grid_spec(0, 3, 3);
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::invalid_argument);

    grid_spec = make_grid_spec(
        std::numeric_limits<std::uint64_t>::max(),
        1,
        1);
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::overflow_error);

    grid_spec = make_grid_spec(
        static_cast<std::uint64_t>(1) << 32U,
        static_cast<std::uint64_t>(1) << 32U,
        1);
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::overflow_error);

    grid_spec = make_grid_spec();
    grid_spec.periodic.x = true;
    grid_spec.reservoir_faces.x_low = true;
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::invalid_argument);

    grid_spec = make_grid_spec();
    grid_spec.periodic.y = true;
    grid_spec.reservoir_faces.y_high = true;
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::invalid_argument);

    grid_spec = make_grid_spec();
    grid_spec.periodic.z = true;
    grid_spec.reservoir_faces.z_low = true;
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::invalid_argument);

    grid_spec = make_grid_spec();
    grid_spec.explicit_source_voxels.push_back({3, 0, 0});
    REQUIRE_THROWS_AS(GasGrid grid(grid_spec), std::invalid_argument);
}

void test_nonperiodic_neighbors()
{
    GasGrid grid(make_grid_spec());

    const auto center_neighbors = grid.neighbors(grid.voxel_id({1, 1, 1}));
    REQUIRE(center_neighbors.count == 6);
    REQUIRE(neighbor_ids(center_neighbors) == ids_for(grid, {
        {0, 1, 1}, {2, 1, 1},
        {1, 0, 1}, {1, 2, 1},
        {1, 1, 0}, {1, 1, 2}
    }));

    const auto corner_neighbors = grid.neighbors(grid.voxel_id({0, 0, 0}));
    REQUIRE(corner_neighbors.count == 3);
    REQUIRE(neighbor_ids(corner_neighbors) == ids_for(grid, {
        {1, 0, 0}, {0, 1, 0}, {0, 0, 1}
    }));
}

void test_all_periodic_axis_combinations()
{
    for (unsigned int mask = 0; mask < 8U; ++mask) {
        auto grid_spec = make_grid_spec();
        grid_spec.periodic = {
            (mask & 1U) != 0U,
            (mask & 2U) != 0U,
            (mask & 4U) != 0U
        };
        GasGrid grid(grid_spec);

        std::vector<VoxelCoord> expected_coords{
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        };
        if (grid_spec.periodic.x) {
            expected_coords.push_back({2, 0, 0});
        }
        if (grid_spec.periodic.y) {
            expected_coords.push_back({0, 2, 0});
        }
        if (grid_spec.periodic.z) {
            expected_coords.push_back({0, 0, 2});
        }

        const auto corner_neighbors = grid.neighbors(grid.voxel_id({0, 0, 0}));
        const auto expected_ids = ids_for(grid, expected_coords);
        REQUIRE(corner_neighbors.count == expected_ids.size());
        REQUIRE(neighbor_ids(corner_neighbors) == expected_ids);

        const auto center_neighbors = grid.neighbors(grid.voxel_id({1, 1, 1}));
        REQUIRE(center_neighbors.count == 6);
        REQUIRE(neighbor_ids(center_neighbors).size() == center_neighbors.count);
    }
}

void test_small_periodic_dimensions_have_unique_neighbors()
{
    auto grid_spec = make_grid_spec(1, 1, 1);
    grid_spec.periodic = {true, true, true};
    GasGrid one_voxel_grid(grid_spec);
    REQUIRE(one_voxel_grid.neighbors(0).count == 0);

    grid_spec = make_grid_spec(2, 1, 1);
    grid_spec.periodic = {true, false, false};
    GasGrid two_voxel_grid(grid_spec);
    const auto line_neighbors = two_voxel_grid.neighbors(0);
    REQUIRE(line_neighbors.count == 1);
    REQUIRE(line_neighbors.ids[0] == 1);

    grid_spec = make_grid_spec(2, 2, 2);
    grid_spec.periodic = {true, true, true};
    GasGrid small_periodic_grid(grid_spec);
    const auto neighbors = small_periodic_grid.neighbors(0);
    REQUIRE(neighbors.count == 3);
    REQUIRE(neighbor_ids(neighbors) == ids_for(small_periodic_grid, {
        {1, 0, 0}, {0, 1, 0}, {0, 0, 1}
    }));
}

void test_reservoir_sources()
{
    auto grid_spec = make_grid_spec();
    grid_spec.reservoir_faces.x_low = true;
    grid_spec.reservoir_faces.z_high = true;
    grid_spec.explicit_source_voxels = {{1, 1, 1}, {1, 1, 1}};
    GasGrid grid(grid_spec);

    REQUIRE(grid.explicit_source_ids().size() == 1);
    REQUIRE(grid.is_reservoir_source({0, 2, 1}));
    REQUIRE(grid.is_reservoir_source({2, 1, 2}));
    REQUIRE(grid.is_reservoir_source({1, 1, 1}));
    REQUIRE(!grid.is_reservoir_source({2, 1, 0}));
    REQUIRE_THROWS_AS(grid.is_reservoir_source({3, 0, 0}), std::out_of_range);

    grid_spec = make_grid_spec();
    grid_spec.periodic = {true, true, true};
    grid_spec.explicit_source_voxels = {{1, 1, 1}};
    GasGrid fully_periodic_grid(grid_spec);
    REQUIRE(fully_periodic_grid.is_reservoir_source({1, 1, 1}));
    REQUIRE(!fully_periodic_grid.is_reservoir_source({0, 0, 0}));
}

void test_world_coordinate_mapping()
{
    auto grid_spec = make_grid_spec(4, 3, 2);
    grid_spec.origin = {1.0, 2.0, 3.0};
    grid_spec.spacing = {0.5, 0.5, 0.5};
    grid_spec.periodic = {true, false, true};
    GasGrid grid(grid_spec);

    REQUIRE((grid.locate_voxel({1.0, 2.0, 3.0}) == VoxelCoord{0, 0, 0}));
    REQUIRE((grid.locate_voxel({2.99, 3.49, 3.99}) == VoxelCoord{3, 2, 1}));
    REQUIRE((grid.locate_voxel({3.0, 2.0, 4.0}) == VoxelCoord{0, 0, 0}));
    REQUIRE((grid.locate_voxel({0.9, 2.0, 2.9}) == VoxelCoord{3, 0, 1}));
    REQUIRE(!grid.locate_voxel({1.0, 1.99, 3.0}));
    REQUIRE(!grid.locate_voxel({1.0, 3.5, 3.0}));
    REQUIRE(!grid.locate_voxel({
        std::numeric_limits<double>::quiet_NaN(),
        2.0,
        3.0
    }));
    REQUIRE(!grid.locate_voxel({
        std::numeric_limits<double>::infinity(),
        2.0,
        3.0
    }));

    for (VoxelId id = 0; id < grid.voxel_count(); ++id) {
        const auto coordinate = grid.voxel_coord(id);
        const auto mapped_coordinate = grid.locate_voxel(grid.voxel_center(id));
        REQUIRE(mapped_coordinate.has_value());
        REQUIRE(*mapped_coordinate == coordinate);
    }
}

void test_non_cubic_world_geometry()
{
    auto grid_spec = make_grid_spec(4, 3, 2);
    grid_spec.origin = {1.0, -2.0, 10.0};
    grid_spec.spacing = {0.5, 1.25, 2.0};
    grid_spec.periodic = {true, false, true};
    GasGrid grid(grid_spec);

    const auto center = grid.voxel_center({2, 1, 0});
    REQUIRE_NEAR(center.x, 2.25, 1.0e-12);
    REQUIRE_NEAR(center.y, -0.125, 1.0e-12);
    REQUIRE_NEAR(center.z, 11.0, 1.0e-12);

    REQUIRE((grid.locate_voxel({2.25, -0.125, 11.0}) == VoxelCoord{2, 1, 0}));
    REQUIRE((grid.locate_voxel({3.0, -2.0, 14.0}) == VoxelCoord{0, 0, 0}));
    REQUIRE((grid.locate_voxel({0.75, -0.75, 9.0}) == VoxelCoord{3, 1, 1}));
    REQUIRE(!grid.locate_voxel({1.0, 1.75, 10.0}));

    for (VoxelId id = 0; id < grid.voxel_count(); ++id) {
        const auto mapped_coordinate = grid.locate_voxel(grid.voxel_center(id));
        REQUIRE(mapped_coordinate.has_value());
        REQUIRE(*mapped_coordinate == grid.voxel_coord(id));
    }
}

}  // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"grid indexing and state storage", test_grid_indexing_and_state_storage},
        {"blocker-count state invariant", test_blocker_count_state_invariant},
        {"invalid grid specifications", test_invalid_grid_specs},
        {"nonperiodic neighbors", test_nonperiodic_neighbors},
        {"all periodic axis combinations", test_all_periodic_axis_combinations},
        {"small periodic dimensions", test_small_periodic_dimensions_have_unique_neighbors},
        {"reservoir sources", test_reservoir_sources},
        {"world coordinate mapping", test_world_coordinate_mapping},
        {"non-cubic world geometry", test_non_cubic_world_geometry}
    };

    std::size_t failure_count = 0;
    for (const auto& test : tests) {
        try {
            test.second();
            std::cout << "[PASS] " << test.first << '\n';
        } catch (const std::exception& error) {
            ++failure_count;
            std::cerr << "[FAIL] " << test.first << ": " << error.what() << '\n';
        }
    }

    if (failure_count != 0) {
        std::cerr << failure_count << " test group(s) failed\n";
        return 1;
    }

    std::cout << tests.size() << " test groups passed\n";
    return 0;
}
