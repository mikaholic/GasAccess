#include "gasaccess/accessibility_query.hpp"

#include <cstddef>
#include <exception>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace allocation_tracking {

bool enabled = false;
std::size_t allocation_count = 0;

}  // namespace allocation_tracking

#if defined(GASACCESS_TEST_WRAP_ALLOCATIONS)
extern "C" void* __real__Znwm(std::size_t size);
extern "C" void* __real__Znam(std::size_t size);

extern "C" void* __wrap__Znwm(std::size_t size)
{
    if (allocation_tracking::enabled) {
        ++allocation_tracking::allocation_count;
    }
    return __real__Znwm(size);
}

extern "C" void* __wrap__Znam(std::size_t size)
{
    if (allocation_tracking::enabled) {
        ++allocation_tracking::allocation_count;
    }
    return __real__Znam(size);
}
#endif

namespace {

using gasaccess::GasAccessibilityQuery;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridSpec;
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
#define REQUIRE_THROWS_AS(expression, exception_type) \
    require_throws<exception_type>([&]() { expression; }, #expression, __FILE__, __LINE__)

GridSpec make_grid_spec(
    std::uint64_t x = 4,
    std::uint64_t y = 3,
    std::uint64_t z = 3)
{
    GridSpec grid_spec{};
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {x, y, z};
    return grid_spec;
}

void test_voxel_and_coordinate_state_queries()
{
    GasGrid gas_grid(make_grid_spec());
    const auto outside_id = gas_grid.voxel_id({2, 1, 1});
    gas_grid.set_gas_state(outside_id, GasState::OutsideAccessible);
    gas_grid.set_gas_state(
        gas_grid.voxel_id({1, 1, 1}),
        GasState::ClosedVoid);

    const GasAccessibilityQuery query(gas_grid);
    REQUIRE(query.is_voxel_outside_accessible(outside_id));
    REQUIRE(!query.is_voxel_outside_accessible(gas_grid.voxel_id({1, 1, 1})));
    REQUIRE(gas_grid.gas_state({2, 1, 1}) == GasState::OutsideAccessible);
    REQUIRE_THROWS_AS(
        query.is_voxel_outside_accessible(gas_grid.voxel_count()),
        std::out_of_range);
    REQUIRE_THROWS_AS(gas_grid.gas_state({4, 0, 0}), std::out_of_range);
}

void test_default_seven_voxel_site_stencil()
{
    GasGrid gas_grid(make_grid_spec());
    gas_grid.fill_gas_state(GasState::Solid);
    const GasAccessibilityQuery query(gas_grid);

    gas_grid.set_gas_state(
        gas_grid.voxel_id({2, 1, 1}),
        GasState::OutsideAccessible);
    REQUIRE(query.is_site_accessible({1.5, 1.5, 1.5}));

    gas_grid.set_gas_state(
        gas_grid.voxel_id({2, 1, 1}),
        GasState::ClosedVoid);
    gas_grid.set_gas_state(
        gas_grid.voxel_id({2, 2, 2}),
        GasState::OutsideAccessible);
    REQUIRE(!query.is_site_accessible({1.5, 1.5, 1.5}));

    gas_grid.set_gas_state(
        gas_grid.voxel_id({1, 1, 1}),
        GasState::OutsideAccessible);
    REQUIRE(query.is_site_accessible({1.5, 1.5, 1.5}));

    REQUIRE(!query.is_site_accessible({-0.1, 1.5, 1.5}));
    REQUIRE(!query.is_site_accessible({4.0, 1.5, 1.5}));
}

void test_periodic_site_stencil()
{
    auto periodic_spec = make_grid_spec();
    periodic_spec.periodic.x = true;
    GasGrid periodic_grid(periodic_spec);
    periodic_grid.fill_gas_state(GasState::Solid);
    periodic_grid.set_gas_state(
        periodic_grid.voxel_id({3, 1, 1}),
        GasState::OutsideAccessible);
    const GasAccessibilityQuery periodic_query(periodic_grid);
    REQUIRE(periodic_query.is_site_accessible({0.5, 1.5, 1.5}));
    REQUIRE(periodic_query.is_site_accessible({4.5, 1.5, 1.5}));

    GasGrid nonperiodic_grid(make_grid_spec());
    nonperiodic_grid.fill_gas_state(GasState::Solid);
    nonperiodic_grid.set_gas_state(
        nonperiodic_grid.voxel_id({3, 1, 1}),
        GasState::OutsideAccessible);
    const GasAccessibilityQuery nonperiodic_query(nonperiodic_grid);
    REQUIRE(!nonperiodic_query.is_site_accessible({0.5, 1.5, 1.5}));
}

void test_explicit_site_stencil()
{
    GasGrid gas_grid(make_grid_spec());
    gas_grid.fill_gas_state(GasState::ClosedVoid);
    const auto outside_id = gas_grid.voxel_id({3, 2, 2});
    gas_grid.set_gas_state(outside_id, GasState::OutsideAccessible);
    const GasAccessibilityQuery query(gas_grid);

    const VoxelId inaccessible_ids[]{
        gas_grid.voxel_id({0, 0, 0}),
        gas_grid.voxel_id({1, 1, 1})
    };
    REQUIRE(!query.is_site_accessible({inaccessible_ids, 2}));

    const VoxelId accessible_ids[]{inaccessible_ids[0], outside_id};
    REQUIRE(query.is_site_accessible({accessible_ids, 2}));
    REQUIRE(!query.is_site_accessible({nullptr, 0}));

    REQUIRE_THROWS_AS(
        query.is_site_accessible({nullptr, 1}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        query.is_site_accessible({nullptr, gasaccess::max_site_voxel_count + 1}),
        std::invalid_argument);

    const VoxelId invalid_id = gas_grid.voxel_count();
    REQUIRE_THROWS_AS(
        query.is_site_accessible({&invalid_id, 1}),
        std::out_of_range);
}

void test_queries_do_not_allocate_or_modify_state()
{
    GasGrid gas_grid(make_grid_spec());
    gas_grid.fill_gas_state(GasState::ClosedVoid);
    const auto outside_id = gas_grid.voxel_id({2, 1, 1});
    gas_grid.set_gas_state(outside_id, GasState::OutsideAccessible);
    const VoxelId stencil[]{gas_grid.voxel_id({1, 1, 1}), outside_id};
    const GasAccessibilityQuery query(gas_grid);

    std::vector<GasState> states_before;
    states_before.reserve(static_cast<std::size_t>(gas_grid.voxel_count()));
    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        states_before.push_back(gas_grid.gas_state(voxel_id));
    }

    allocation_tracking::allocation_count = 0;
    allocation_tracking::enabled = true;
    bool all_accessible = true;
    for (int iteration = 0; iteration < 1000; ++iteration) {
        all_accessible = all_accessible
            && query.is_voxel_outside_accessible(outside_id)
            && query.is_site_accessible({1.5, 1.5, 1.5})
            && query.is_site_accessible({stencil, 2});
    }
    allocation_tracking::enabled = false;

    REQUIRE(all_accessible);
#if defined(GASACCESS_TEST_WRAP_ALLOCATIONS)
    REQUIRE(allocation_tracking::allocation_count == 0);
#endif
    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        REQUIRE(gas_grid.gas_state(voxel_id)
            == states_before[static_cast<std::size_t>(voxel_id)]);
    }
}

}  // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"voxel and coordinate state queries", test_voxel_and_coordinate_state_queries},
        {"default seven-voxel site stencil", test_default_seven_voxel_site_stencil},
        {"periodic site stencil", test_periodic_site_stencil},
        {"explicit site stencil", test_explicit_site_stencil},
        {"queries do not allocate or modify state",
         test_queries_do_not_allocate_or_modify_state}
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
