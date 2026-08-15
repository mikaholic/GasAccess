#include "gasaccess/deposition_updater.hpp"
#include "gasaccess/exterior_classifier.hpp"
#include "gasaccess/local_topology_filter.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using gasaccess::Atom;
using gasaccess::AtomVoxelizer;
using gasaccess::DepositionUpdater;
using gasaccess::ExteriorClassifier;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridSpec;
using gasaccess::LocalTopologyFilter;
using gasaccess::RemovedVoxel;
using gasaccess::TopologyDecision;
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
    std::uint64_t x,
    std::uint64_t y,
    std::uint64_t z)
{
    GridSpec grid_spec{};
    grid_spec.spacing = 1.0;
    grid_spec.dimensions = {x, y, z};
    return grid_spec;
}

std::vector<GasState> copy_states(const GasGrid& gas_grid)
{
    std::vector<GasState> states;
    states.reserve(static_cast<std::size_t>(gas_grid.voxel_count()));
    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        states.push_back(gas_grid.gas_state(voxel_id));
    }
    return states;
}

void require_same_states(const GasGrid& lhs, const GasGrid& rhs)
{
    REQUIRE(copy_states(lhs) == copy_states(rhs));
}

void test_input_validation_and_multi_voxel_fallback()
{
    auto grid_spec = make_grid_spec(3, 3, 3);
    grid_spec.reservoir_faces.z_high = true;
    GasGrid gas_grid(grid_spec);
    ExteriorClassifier{}.classify(gas_grid);
    const LocalTopologyFilter topology_filter;

    REQUIRE(topology_filter.evaluate(gas_grid, {nullptr, 0}).is_safe());
    REQUIRE_THROWS_AS(
        topology_filter.evaluate(gas_grid, {nullptr, 1}),
        std::invalid_argument);

    const auto first_id = gas_grid.voxel_id({1, 1, 1});
    RemovedVoxel first_removed{first_id, GasState::OutsideAccessible};
    REQUIRE_THROWS_AS(
        topology_filter.evaluate(gas_grid, {&first_removed, 1}),
        std::invalid_argument);

    gas_grid.set_gas_state(first_id, GasState::Solid);
    first_removed.previous_state = GasState::Solid;
    REQUIRE_THROWS_AS(
        topology_filter.evaluate(gas_grid, {&first_removed, 1}),
        std::invalid_argument);

    first_removed.previous_state = GasState::OutsideAccessible;
    const auto second_id = gas_grid.voxel_id({1, 1, 2});
    gas_grid.set_gas_state(second_id, GasState::Solid);
    const std::vector<RemovedVoxel> removed_voxels{
        first_removed,
        {second_id, GasState::OutsideAccessible}
    };
    const auto result = topology_filter.evaluate(
        gas_grid,
        {removed_voxels.data(), removed_voxels.size()});
    REQUIRE(result.decision == TopologyDecision::RequiresConnectivityRepair);
}

void test_open_region_is_proven_safe()
{
    auto grid_spec = make_grid_spec(5, 5, 5);
    grid_spec.reservoir_faces.z_high = true;
    GasGrid gas_grid(grid_spec);
    ExteriorClassifier{}.classify(gas_grid);
    const auto removed_id = gas_grid.voxel_id({2, 2, 2});
    gas_grid.set_gas_state(removed_id, GasState::Solid);
    const RemovedVoxel removed_voxel{removed_id, GasState::OutsideAccessible};

    const auto result = LocalTopologyFilter{}.evaluate(
        gas_grid,
        {&removed_voxel, 1});
    REQUIRE(result.is_safe());
    REQUIRE(result.accessible_neighbor_count == 6);
    REQUIRE(result.visited_voxel_count > result.accessible_neighbor_count);
}

void test_bridge_and_source_are_inconclusive()
{
    auto bridge_spec = make_grid_spec(5, 1, 1);
    bridge_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid bridge_grid(bridge_spec);
    ExteriorClassifier{}.classify(bridge_grid);
    const auto bridge_id = bridge_grid.voxel_id({2, 0, 0});
    bridge_grid.set_gas_state(bridge_id, GasState::Solid);
    const RemovedVoxel bridge_voxel{bridge_id, GasState::OutsideAccessible};
    const auto bridge_result = LocalTopologyFilter{}.evaluate(
        bridge_grid,
        {&bridge_voxel, 1});
    REQUIRE(bridge_result.decision == TopologyDecision::RequiresConnectivityRepair);
    REQUIRE(bridge_result.accessible_neighbor_count == 2);

    auto source_spec = make_grid_spec(3, 3, 3);
    source_spec.explicit_source_voxels = {{1, 1, 1}};
    GasGrid source_grid(source_spec);
    ExteriorClassifier{}.classify(source_grid);
    const auto source_id = source_grid.voxel_id({1, 1, 1});
    source_grid.set_gas_state(source_id, GasState::Solid);
    const RemovedVoxel source_voxel{source_id, GasState::OutsideAccessible};
    const auto source_result = LocalTopologyFilter{}.evaluate(
        source_grid,
        {&source_voxel, 1});
    REQUIRE(source_result.decision == TopologyDecision::RequiresConnectivityRepair);
}

void test_closed_void_removal_is_safe()
{
    GasGrid gas_grid(make_grid_spec(3, 3, 3));
    ExteriorClassifier{}.classify(gas_grid);
    const auto removed_id = gas_grid.voxel_id({1, 1, 1});
    REQUIRE(gas_grid.gas_state(removed_id) == GasState::ClosedVoid);
    gas_grid.set_gas_state(removed_id, GasState::Solid);
    const RemovedVoxel removed_voxel{removed_id, GasState::ClosedVoid};

    const auto result = LocalTopologyFilter{}.evaluate(
        gas_grid,
        {&removed_voxel, 1});
    REQUIRE(result.is_safe());
    REQUIRE(result.accessible_neighbor_count == 0);
}

void test_periodic_safe_path_and_periodic_bridge()
{
    auto safe_spec = make_grid_spec(3, 3, 1);
    safe_spec.periodic = {true, true, false};
    safe_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid safe_grid(safe_spec);
    ExteriorClassifier{}.classify(safe_grid);
    const auto safe_id = safe_grid.voxel_id({0, 1, 0});
    safe_grid.set_gas_state(safe_id, GasState::Solid);
    const RemovedVoxel safe_voxel{safe_id, GasState::OutsideAccessible};
    const auto safe_result = LocalTopologyFilter{}.evaluate(
        safe_grid,
        {&safe_voxel, 1});
    REQUIRE(safe_result.is_safe());
    REQUIRE(safe_result.accessible_neighbor_count == 4);

    auto bridge_spec = make_grid_spec(5, 1, 1);
    bridge_spec.periodic.x = true;
    bridge_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid bridge_grid(bridge_spec);
    bridge_grid.set_gas_state(
        bridge_grid.voxel_id({1, 0, 0}),
        GasState::Solid);
    ExteriorClassifier{}.classify(bridge_grid);
    const auto bridge_id = bridge_grid.voxel_id({4, 0, 0});
    bridge_grid.set_gas_state(bridge_id, GasState::Solid);
    const RemovedVoxel bridge_voxel{bridge_id, GasState::OutsideAccessible};
    const auto bridge_result = LocalTopologyFilter{}.evaluate(
        bridge_grid,
        {&bridge_voxel, 1});
    REQUIRE(bridge_result.decision == TopologyDecision::RequiresConnectivityRepair);
    REQUIRE(bridge_result.accessible_neighbor_count == 2);
}

void test_exhaustive_planar_patterns_match_full_reference()
{
    std::size_t safe_count = 0;
    std::size_t repair_count = 0;
    constexpr std::uint32_t pattern_count = 128;
    const std::array<VoxelCoord, 7> pattern_voxels{{
        {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {2, 1, 0},
        {0, 2, 0}, {1, 2, 0}, {2, 2, 0}
    }};

    for (std::uint32_t pattern = 0; pattern < pattern_count; ++pattern) {
        auto grid_spec = make_grid_spec(3, 3, 1);
        grid_spec.explicit_source_voxels = {{0, 0, 0}};
        GasGrid actual_grid(grid_spec);
        GasGrid reference_grid(grid_spec);
        for (std::size_t bit = 0; bit < pattern_voxels.size(); ++bit) {
            if ((pattern & (UINT32_C(1) << bit)) == 0U) {
                continue;
            }
            actual_grid.set_gas_state(
                actual_grid.voxel_id(pattern_voxels[bit]),
                GasState::Solid);
            reference_grid.set_gas_state(
                reference_grid.voxel_id(pattern_voxels[bit]),
                GasState::Solid);
        }
        ExteriorClassifier{}.classify(actual_grid);
        ExteriorClassifier{}.classify(reference_grid);

        const Atom atom{{1.5, 1.5, 0.5}, 0.0};
        const auto actual_result = DepositionUpdater(0.0).apply_deposition(
            actual_grid,
            {&atom, 1});
        AtomVoxelizer(0.0).voxelize(reference_grid, {&atom, 1});
        const auto reference_summary = ExteriorClassifier{}.classify(reference_grid);

        require_same_states(actual_grid, reference_grid);
        REQUIRE(actual_result.classification.solid_count
            == reference_summary.solid_count);
        REQUIRE(actual_result.classification.outside_accessible_count
            == reference_summary.outside_accessible_count);
        REQUIRE(actual_result.classification.closed_void_count
            == reference_summary.closed_void_count);
        if (actual_result.used_affected_region_repair()) {
            ++repair_count;
        } else {
            ++safe_count;
        }
    }
    REQUIRE(safe_count != 0);
    REQUIRE(repair_count != 0);
}

}  // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"input validation and multi-voxel fallback",
         test_input_validation_and_multi_voxel_fallback},
        {"open region is locally safe", test_open_region_is_proven_safe},
        {"bridge and source are inconclusive", test_bridge_and_source_are_inconclusive},
        {"closed-void removal is safe", test_closed_void_removal_is_safe},
        {"periodic safe path and bridge", test_periodic_safe_path_and_periodic_bridge},
        {"exhaustive planar differential patterns",
         test_exhaustive_planar_patterns_match_full_reference}
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
