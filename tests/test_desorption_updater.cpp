#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/deposition_updater.hpp"
#include "gasaccess/desorption_updater.hpp"
#include "gasaccess/opening_region_repair.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using gasaccess::Atom;
using gasaccess::AtomView;
using gasaccess::AtomVoxelizer;
using gasaccess::ClassificationSummary;
using gasaccess::DepositionUpdater;
using gasaccess::DesorptionRepairMode;
using gasaccess::DesorptionUpdateResult;
using gasaccess::DesorptionUpdater;
using gasaccess::ExteriorClassifier;
using gasaccess::GasAccessibilityQuery;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridSpec;
using gasaccess::OpeningRegionRepair;
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
void require_throws(
    Function&& function,
    const char* expression,
    const char* file,
    int line)
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
    require_throws<exception_type>( \
        [&]() { expression; }, #expression, __FILE__, __LINE__)

GridSpec make_grid_spec(
    std::uint64_t x,
    std::uint64_t y,
    std::uint64_t z)
{
    GridSpec grid_spec{};
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {x, y, z};
    return grid_spec;
}

void set_solid(GasGrid& gas_grid, const VoxelCoord& coordinate)
{
    gas_grid.set_gas_state(gas_grid.voxel_id(coordinate), GasState::Solid);
}

std::vector<GasState> copy_states(const GasGrid& gas_grid)
{
    std::vector<GasState> states;
    states.reserve(static_cast<std::size_t>(gas_grid.voxel_count()));
    for (VoxelId id = 0; id < gas_grid.voxel_count(); ++id) {
        states.push_back(gas_grid.gas_state(id));
    }
    return states;
}

std::vector<VoxelId> changed_state_ids(
    const std::vector<GasState>& previous_states,
    const GasGrid& gas_grid)
{
    std::vector<VoxelId> result;
    for (VoxelId id = 0; id < gas_grid.voxel_count(); ++id) {
        if (previous_states[static_cast<std::size_t>(id)]
            != gas_grid.gas_state(id)) {
            result.push_back(id);
        }
    }
    return result;
}

bool summaries_equal(
    const ClassificationSummary& lhs,
    const ClassificationSummary& rhs)
{
    return lhs.solid_count == rhs.solid_count
        && lhs.outside_accessible_count == rhs.outside_accessible_count
        && lhs.closed_void_count == rhs.closed_void_count;
}

void require_valid_changed_ids(const std::vector<VoxelId>& changed_voxel_ids)
{
    REQUIRE(std::is_sorted(changed_voxel_ids.begin(), changed_voxel_ids.end()));
    REQUIRE(std::adjacent_find(
        changed_voxel_ids.begin(),
        changed_voxel_ids.end()) == changed_voxel_ids.end());
}

void require_same_grid(const GasGrid& actual, const GasGrid& expected)
{
    REQUIRE(actual.voxel_count() == expected.voxel_count());
    for (VoxelId id = 0; id < actual.voxel_count(); ++id) {
        REQUIRE(actual.blocker_count(id) == expected.blocker_count(id));
        REQUIRE(actual.gas_state(id) == expected.gas_state(id));
    }
    for (const auto state : {
             GasState::Unclassified,
             GasState::Solid,
             GasState::OutsideAccessible,
             GasState::ClosedVoid}) {
        REQUIRE(actual.gas_state_count(state) == expected.gas_state_count(state));
    }
}

void require_matches_full_reclassification(
    GasGrid& actual_grid,
    GasGrid& reference_grid,
    AtomView removed_atoms,
    DesorptionUpdateResult& actual_result)
{
    const DesorptionUpdater incremental_updater(0.0);
    const DesorptionUpdater reference_updater(
        0.0,
        DesorptionRepairMode::FullReclassification);
    actual_result = incremental_updater.apply_desorption(
        actual_grid,
        removed_atoms);
    const auto reference_result = reference_updater.apply_desorption(
        reference_grid,
        removed_atoms);

    REQUIRE(actual_result.newly_gas_count == reference_result.newly_gas_count);
    REQUIRE(actual_result.blocker_count_changed_voxel_count
        == reference_result.blocker_count_changed_voxel_count);
    REQUIRE(actual_result.changed_voxel_ids == reference_result.changed_voxel_ids);
    REQUIRE(summaries_equal(
        actual_result.classification,
        reference_result.classification));
    require_valid_changed_ids(actual_result.changed_voxel_ids);
    require_same_grid(actual_grid, reference_grid);
}

void test_validation_and_empty_update()
{
    auto grid_spec = make_grid_spec(3, 1, 1);
    grid_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid unclassified_grid(grid_spec);
    const Atom atom{{1.5, 0.5, 0.5}, 0.0};
    AtomVoxelizer(0.0).voxelize(unclassified_grid, {&atom, 1});
    REQUIRE_THROWS_AS(
        DesorptionUpdater(0.0).apply_desorption(
            unclassified_grid,
            {&atom, 1}),
        std::invalid_argument);
    REQUIRE(unclassified_grid.blocker_count({1, 0, 0}) == 1);

    GasGrid gas_grid(grid_spec);
    ExteriorClassifier{}.classify(gas_grid);
    const auto states_before = copy_states(gas_grid);
    const auto empty_result = DesorptionUpdater(0.0).apply_desorption(
        gas_grid,
        {nullptr, 0});
    REQUIRE(!empty_result.geometry_changed());
    REQUIRE(empty_result.blocker_count_changed_voxel_count == 0);
    REQUIRE(empty_result.changed_voxel_ids.empty());
    REQUIRE(copy_states(gas_grid) == states_before);

    REQUIRE_THROWS_AS(
        DesorptionUpdater(0.0).apply_desorption(gas_grid, {nullptr, 1}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        DesorptionUpdater(0.0).apply_desorption(gas_grid, {&atom, 1}),
        std::underflow_error);
    REQUIRE(copy_states(gas_grid) == states_before);

    OpeningRegionRepair repair;
    const auto repair_result = repair.repair(gas_grid, {nullptr, 0});
    REQUIRE(repair_result.visited_voxel_count == 0);
    REQUIRE(repair_result.newly_opened_voxel_ids.empty());
    REQUIRE_THROWS_AS(repair.repair(gas_grid, {nullptr, 1}), std::invalid_argument);
    const auto outside_id = gas_grid.voxel_id({0, 0, 0});
    REQUIRE_THROWS_AS(repair.repair(gas_grid, {&outside_id, 1}), std::invalid_argument);

    REQUIRE_THROWS_AS(
        DesorptionUpdater(0.0, static_cast<DesorptionRepairMode>(99)),
        std::invalid_argument);
}

void test_overlapping_blocker_fast_path_and_final_removal()
{
    auto grid_spec = make_grid_spec(5, 1, 1);
    grid_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid gas_grid(grid_spec);
    const Atom atom{{2.5, 0.5, 0.5}, 0.0};
    const AtomVoxelizer voxelizer(0.0);
    voxelizer.voxelize(gas_grid, {&atom, 1});
    voxelizer.voxelize(gas_grid, {&atom, 1});
    ExteriorClassifier{}.classify(gas_grid);

    const DesorptionUpdater updater(0.0);
    auto result = updater.apply_desorption(gas_grid, {&atom, 1});
    REQUIRE(result.blocker_count_changed_voxel_count == 1);
    REQUIRE(result.newly_gas_count == 0);
    REQUIRE(!result.geometry_changed());
    REQUIRE(!result.used_opening_region_repair());
    REQUIRE(result.changed_voxel_ids.empty());
    REQUIRE(gas_grid.blocker_count({2, 0, 0}) == 1);
    REQUIRE(gas_grid.gas_state({2, 0, 0}) == GasState::Solid);

    result = updater.apply_desorption(gas_grid, {&atom, 1});
    REQUIRE(result.newly_gas_count == 1);
    REQUIRE(result.geometry_changed());
    REQUIRE(result.used_opening_region_repair());
    REQUIRE(result.opening_visited_voxel_count == 3);
    REQUIRE(result.repair_opened_voxel_count == 3);
    REQUIRE(result.changed_voxel_ids == std::vector<VoxelId>({
        gas_grid.voxel_id({2, 0, 0}),
        gas_grid.voxel_id({3, 0, 0}),
        gas_grid.voxel_id({4, 0, 0})
    }));
    REQUIRE(gas_grid.gas_state({2, 0, 0}) == GasState::OutsideAccessible);

    const auto states_before = copy_states(gas_grid);
    REQUIRE_THROWS_AS(
        updater.apply_desorption(gas_grid, {&atom, 1}),
        std::underflow_error);
    REQUIRE(copy_states(gas_grid) == states_before);
}

void test_new_void_remains_closed()
{
    auto grid_spec = make_grid_spec(3, 3, 3);
    GasGrid gas_grid(grid_spec);
    gas_grid.fill_gas_state(GasState::Solid);
    ExteriorClassifier{}.classify(gas_grid);
    const Atom atom{gas_grid.voxel_center({1, 1, 1}), 0.0};

    const auto result = DesorptionUpdater(0.0).apply_desorption(
        gas_grid,
        {&atom, 1});
    REQUIRE(result.newly_gas_count == 1);
    REQUIRE(!result.used_opening_region_repair());
    REQUIRE(result.opening_visited_voxel_count == 0);
    REQUIRE(result.repair_opened_voxel_count == 0);
    REQUIRE(result.changed_voxel_ids == std::vector<VoxelId>({
        gas_grid.voxel_id({1, 1, 1})
    }));
    REQUIRE(gas_grid.gas_state({1, 1, 1}) == GasState::ClosedVoid);
    REQUIRE(result.classification.solid_count == 26);
    REQUIRE(result.classification.outside_accessible_count == 0);
    REQUIRE(result.classification.closed_void_count == 1);
}

void test_local_opening_visits_only_newly_gas_voxel()
{
    auto grid_spec = make_grid_spec(3, 3, 1);
    grid_spec.reservoir_faces.x_low = true;
    GasGrid gas_grid(grid_spec);
    set_solid(gas_grid, {1, 1, 0});
    ExteriorClassifier{}.classify(gas_grid);
    REQUIRE(gas_grid.gas_state({2, 1, 0}) == GasState::OutsideAccessible);
    const Atom atom{{1.5, 1.5, 0.5}, 0.0};

    const auto result = DesorptionUpdater(0.0).apply_desorption(
        gas_grid,
        {&atom, 1});
    REQUIRE(result.newly_gas_count == 1);
    REQUIRE(result.opening_visited_voxel_count == 1);
    REQUIRE(result.repair_opened_voxel_count == 1);
    REQUIRE(result.changed_voxel_ids == std::vector<VoxelId>({
        gas_grid.voxel_id({1, 1, 0})
    }));
}

void test_one_voxel_opening_exposes_cavity()
{
    auto grid_spec = make_grid_spec(5, 1, 1);
    grid_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid actual_grid(grid_spec);
    GasGrid reference_grid(grid_spec);
    set_solid(actual_grid, {1, 0, 0});
    set_solid(reference_grid, {1, 0, 0});
    ExteriorClassifier{}.classify(actual_grid);
    ExteriorClassifier{}.classify(reference_grid);
    const Atom barrier{{1.5, 0.5, 0.5}, 0.0};

    DesorptionUpdateResult result{};
    require_matches_full_reclassification(
        actual_grid,
        reference_grid,
        {&barrier, 1},
        result);
    REQUIRE(result.used_opening_region_repair());
    REQUIRE(result.opening_visited_voxel_count == 4);
    REQUIRE(result.repair_opened_voxel_count == 4);
    REQUIRE(result.changed_voxel_ids == std::vector<VoxelId>({1, 2, 3, 4}));
}

void test_small_opening_does_not_scan_large_grid()
{
    auto grid_spec = make_grid_spec(10000, 1, 1);
    grid_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid actual_grid(grid_spec);
    GasGrid reference_grid(grid_spec);
    set_solid(actual_grid, {10, 0, 0});
    set_solid(actual_grid, {21, 0, 0});
    set_solid(reference_grid, {10, 0, 0});
    set_solid(reference_grid, {21, 0, 0});
    ExteriorClassifier{}.classify(actual_grid);
    ExteriorClassifier{}.classify(reference_grid);
    const Atom opening{{10.5, 0.5, 0.5}, 0.0};

    DesorptionUpdateResult result{};
    require_matches_full_reclassification(
        actual_grid,
        reference_grid,
        {&opening, 1},
        result);
    REQUIRE(result.opening_visited_voxel_count == 11);
    REQUIRE(result.opening_visited_voxel_count < actual_grid.voxel_count() / 100);
    REQUIRE(actual_grid.gas_state({20, 0, 0})
        == GasState::OutsideAccessible);
    REQUIRE(actual_grid.gas_state({22, 0, 0}) == GasState::ClosedVoid);
}

void test_batched_openings_expose_disconnected_cavities()
{
    auto grid_spec = make_grid_spec(5, 3, 1);
    grid_spec.explicit_source_voxels = {{0, 0, 0}, {0, 2, 0}};
    GasGrid actual_grid(grid_spec);
    GasGrid reference_grid(grid_spec);
    for (std::int64_t x = 0; x < 5; ++x) {
        set_solid(actual_grid, {x, 1, 0});
        set_solid(reference_grid, {x, 1, 0});
    }
    set_solid(actual_grid, {1, 0, 0});
    set_solid(actual_grid, {1, 2, 0});
    set_solid(reference_grid, {1, 0, 0});
    set_solid(reference_grid, {1, 2, 0});
    ExteriorClassifier{}.classify(actual_grid);
    ExteriorClassifier{}.classify(reference_grid);
    const std::vector<Atom> barriers{
        {{1.5, 0.5, 0.5}, 0.0},
        {{1.5, 2.5, 0.5}, 0.0}
    };

    DesorptionUpdateResult result{};
    require_matches_full_reclassification(
        actual_grid,
        reference_grid,
        {barriers.data(), barriers.size()},
        result);
    REQUIRE(result.newly_gas_count == 2);
    REQUIRE(result.opening_visited_voxel_count == 8);
    REQUIRE(result.repair_opened_voxel_count == 8);
    REQUIRE(result.changed_voxel_ids.size() == 8);
}

void test_desorption_exposes_explicit_source()
{
    auto grid_spec = make_grid_spec(3, 1, 1);
    grid_spec.explicit_source_voxels = {{1, 0, 0}};
    GasGrid gas_grid(grid_spec);
    gas_grid.fill_gas_state(GasState::Solid);
    ExteriorClassifier{}.classify(gas_grid);
    const Atom source_atom{{1.5, 0.5, 0.5}, 0.0};

    const auto result = DesorptionUpdater(0.0).apply_desorption(
        gas_grid,
        {&source_atom, 1});
    REQUIRE(result.newly_gas_count == 1);
    REQUIRE(result.opening_visited_voxel_count == 1);
    REQUIRE(result.repair_opened_voxel_count == 1);
    REQUIRE(gas_grid.gas_state({1, 0, 0}) == GasState::OutsideAccessible);
}

VoxelCoord axis_coordinate(unsigned int axis, std::int64_t value)
{
    if (axis == 0) {
        return {value, 0, 0};
    }
    if (axis == 1) {
        return {0, value, 0};
    }
    return {0, 0, value};
}

void test_opening_crosses_each_periodic_seam()
{
    for (unsigned int axis = 0; axis < 3; ++axis) {
        auto dimensions = gasaccess::GridDimensions{1, 1, 1};
        if (axis == 0) {
            dimensions.x = 7;
        } else if (axis == 1) {
            dimensions.y = 7;
        } else {
            dimensions.z = 7;
        }
        auto grid_spec = make_grid_spec(
            dimensions.x,
            dimensions.y,
            dimensions.z);
        grid_spec.periodic = {
            axis == 0,
            axis == 1,
            axis == 2
        };
        grid_spec.explicit_source_voxels = {axis_coordinate(axis, 0)};
        GasGrid actual_grid(grid_spec);
        GasGrid reference_grid(grid_spec);
        set_solid(actual_grid, axis_coordinate(axis, 1));
        set_solid(actual_grid, axis_coordinate(axis, 6));
        set_solid(reference_grid, axis_coordinate(axis, 1));
        set_solid(reference_grid, axis_coordinate(axis, 6));
        ExteriorClassifier{}.classify(actual_grid);
        ExteriorClassifier{}.classify(reference_grid);
        const Atom seam_barrier{
            actual_grid.voxel_center(axis_coordinate(axis, 6)),
            0.0
        };

        DesorptionUpdateResult result{};
        require_matches_full_reclassification(
            actual_grid,
            reference_grid,
            {&seam_barrier, 1},
            result);
        REQUIRE(result.opening_visited_voxel_count == 5);
        for (std::int64_t index = 2; index < 7; ++index) {
            REQUIRE(actual_grid.gas_state(axis_coordinate(axis, index))
                == GasState::OutsideAccessible);
        }
    }
}

void test_random_deposition_desorption_sequences_match_rebuild()
{
    std::mt19937_64 random_engine(0xbb67ae8584caa73bULL);
    std::uniform_int_distribution<int> x_distribution(0, 4);
    std::uniform_int_distribution<int> y_distribution(0, 3);
    std::uniform_int_distribution<int> z_distribution(0, 2);
    std::uniform_int_distribution<std::size_t> choice_distribution;
    std::bernoulli_distribution add_distribution(0.55);

    auto grid_spec = make_grid_spec(5, 4, 3);
    grid_spec.periodic = {true, true, false};
    grid_spec.reservoir_faces.z_high = true;
    GasGrid actual_grid(grid_spec);
    ExteriorClassifier{}.classify(actual_grid);
    const DepositionUpdater deposition_updater(0.0);
    const DesorptionUpdater desorption_updater(0.0);
    const AtomVoxelizer voxelizer(0.0);
    std::vector<Atom> active_atoms;

    for (int event = 0; event < 160; ++event) {
        const auto previous_states = copy_states(actual_grid);
        std::vector<VoxelId> actual_changed_ids;
        ClassificationSummary actual_summary{};
        if (active_atoms.empty() || add_distribution(random_engine)) {
            const VoxelCoord coordinate{
                x_distribution(random_engine),
                y_distribution(random_engine),
                z_distribution(random_engine)
            };
            const Atom atom{actual_grid.voxel_center(coordinate), 0.0};
            const auto result = deposition_updater.apply_deposition(
                actual_grid,
                {&atom, 1});
            active_atoms.push_back(atom);
            actual_changed_ids = result.changed_voxel_ids;
            actual_summary = result.classification;
        } else {
            choice_distribution = std::uniform_int_distribution<std::size_t>(
                0,
                active_atoms.size() - 1);
            const auto atom_index = choice_distribution(random_engine);
            const auto atom = active_atoms[atom_index];
            const auto result = desorption_updater.apply_desorption(
                actual_grid,
                {&atom, 1});
            active_atoms.erase(
                active_atoms.begin() + static_cast<std::ptrdiff_t>(atom_index));
            actual_changed_ids = result.changed_voxel_ids;
            actual_summary = result.classification;
        }

        GasGrid rebuilt_grid(grid_spec);
        voxelizer.voxelize(
            rebuilt_grid,
            {active_atoms.data(), active_atoms.size()});
        const auto rebuilt_summary = ExteriorClassifier{}.classify(rebuilt_grid);
        const auto expected_changed_ids = changed_state_ids(
            previous_states,
            rebuilt_grid);

        REQUIRE(actual_changed_ids == expected_changed_ids);
        REQUIRE(summaries_equal(actual_summary, rebuilt_summary));
        require_valid_changed_ids(actual_changed_ids);
        require_same_grid(actual_grid, rebuilt_grid);
        const GasAccessibilityQuery actual_query(actual_grid);
        const GasAccessibilityQuery rebuilt_query(rebuilt_grid);
        for (VoxelId id = 0; id < actual_grid.voxel_count(); ++id) {
            const auto point = actual_grid.voxel_center(id);
            REQUIRE(actual_query.is_site_accessible(point)
                == rebuilt_query.is_site_accessible(point));
        }
    }
}

}  // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"validation and empty update", test_validation_and_empty_update},
        {"overlap fast path and final removal",
         test_overlapping_blocker_fast_path_and_final_removal},
        {"new void remains closed", test_new_void_remains_closed},
        {"local opening visits only newly gas voxel",
         test_local_opening_visits_only_newly_gas_voxel},
        {"one-voxel opening exposes cavity",
         test_one_voxel_opening_exposes_cavity},
        {"small opening avoids full-grid scan",
         test_small_opening_does_not_scan_large_grid},
        {"batched openings expose disconnected cavities",
         test_batched_openings_expose_disconnected_cavities},
        {"desorption exposes explicit source",
         test_desorption_exposes_explicit_source},
        {"opening crosses periodic seams",
         test_opening_crosses_each_periodic_seam},
        {"random deposition/desorption sequences match rebuild",
         test_random_deposition_desorption_sequences_match_rebuild}
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
