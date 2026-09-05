#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/atom_change_event_buffer.hpp"
#include "gasaccess/atom_change_updater.hpp"
#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/exterior_classifier.hpp"

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

using gasaccess::AccessibilityRepairKind;
using gasaccess::Atom;
using gasaccess::AtomChangeBatch;
using gasaccess::AtomChangeEventBuffer;
using gasaccess::AtomChangeRepairMode;
using gasaccess::AtomChangeUpdateResult;
using gasaccess::AtomChangeUpdater;
using gasaccess::AtomView;
using gasaccess::AtomVoxelizer;
using gasaccess::ExteriorClassifier;
using gasaccess::GasAccessibilityQuery;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridSpec;
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
    std::uint64_t x = 12,
    std::uint64_t y = 6,
    std::uint64_t z = 8)
{
    GridSpec grid_spec{};
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {x, y, z};
    grid_spec.reservoir_faces.z_high = true;
    return grid_spec;
}

std::vector<Atom> atoms_for_voxels(
    const GasGrid& gas_grid,
    const std::vector<VoxelCoord>& voxel_coords,
    double radius = 0.0)
{
    std::vector<Atom> atoms;
    atoms.reserve(voxel_coords.size());
    for (const auto& voxel_coord : voxel_coords) {
        atoms.push_back({gas_grid.voxel_center(voxel_coord), radius});
    }
    return atoms;
}

std::vector<GasState> copy_states(const GasGrid& gas_grid)
{
    std::vector<GasState> states;
    states.reserve(static_cast<std::size_t>(gas_grid.voxel_count()));
    for (VoxelId voxel_id = 0;
         voxel_id < gas_grid.voxel_count();
         ++voxel_id) {
        states.push_back(gas_grid.gas_state(voxel_id));
    }
    return states;
}

std::vector<VoxelId> changed_state_ids(
    const std::vector<GasState>& previous_states,
    const GasGrid& gas_grid)
{
    std::vector<VoxelId> changed;
    for (VoxelId voxel_id = 0;
         voxel_id < gas_grid.voxel_count();
         ++voxel_id) {
        if (previous_states[static_cast<std::size_t>(voxel_id)]
            != gas_grid.gas_state(voxel_id)) {
            changed.push_back(voxel_id);
        }
    }
    return changed;
}

void require_same_grid(const GasGrid& lhs, const GasGrid& rhs)
{
    REQUIRE(lhs.voxel_count() == rhs.voxel_count());
    const GasAccessibilityQuery lhs_query(lhs);
    const GasAccessibilityQuery rhs_query(rhs);
    for (VoxelId voxel_id = 0; voxel_id < lhs.voxel_count(); ++voxel_id) {
        REQUIRE(lhs.blocker_count(voxel_id) == rhs.blocker_count(voxel_id));
        REQUIRE(lhs.gas_state(voxel_id) == rhs.gas_state(voxel_id));
        const auto center = lhs.voxel_center(voxel_id);
        REQUIRE(lhs_query.is_site_accessible(center)
            == rhs_query.is_site_accessible(center));
    }
    for (const auto state : {
             GasState::Unclassified,
             GasState::Solid,
             GasState::OutsideAccessible,
             GasState::ClosedVoid}) {
        REQUIRE(lhs.gas_state_count(state) == rhs.gas_state_count(state));
    }
}

void require_sorted_unique(const std::vector<VoxelId>& voxel_ids)
{
    REQUIRE(std::is_sorted(voxel_ids.begin(), voxel_ids.end()));
    REQUIRE(std::adjacent_find(voxel_ids.begin(), voxel_ids.end())
        == voxel_ids.end());
}

struct UpdatePair {
    AtomChangeUpdateResult incremental{};
    AtomChangeUpdateResult full{};
};

UpdatePair apply_and_compare(
    GasGrid& incremental_grid,
    GasGrid& full_grid,
    const AtomChangeBatch& changes,
    AccessibilityRepairKind expected_kind)
{
    const auto previous_states = copy_states(incremental_grid);
    UpdatePair results{};
    results.incremental = AtomChangeUpdater(0.0).apply_atom_changes(
        incremental_grid,
        changes);
    results.full = AtomChangeUpdater(
        0.0,
        AtomChangeRepairMode::FullReclassification).apply_atom_changes(
            full_grid,
            changes);

    REQUIRE(results.incremental.blocker_count_changed_voxel_count
        == results.full.blocker_count_changed_voxel_count);
    REQUIRE(results.incremental.newly_solid_count
        == results.full.newly_solid_count);
    REQUIRE(results.incremental.newly_gas_count
        == results.full.newly_gas_count);
    REQUIRE(results.incremental.repair_kind == expected_kind);
    REQUIRE(results.full.repair_kind
        == (results.full.geometry_changed()
            ? AccessibilityRepairKind::FullReclassification
            : AccessibilityRepairKind::None));
    REQUIRE(results.full.used_full_reclassification()
        == results.full.geometry_changed());
    const auto expected_changed = changed_state_ids(
        previous_states,
        incremental_grid);
    REQUIRE(results.incremental.changed_voxel_ids == expected_changed);
    REQUIRE(results.full.changed_voxel_ids == expected_changed);
    require_sorted_unique(results.incremental.changed_voxel_ids);
    require_sorted_unique(results.full.changed_voxel_ids);
    require_same_grid(incremental_grid, full_grid);
    return results;
}

void initialize_grids(
    GasGrid& incremental_grid,
    GasGrid& full_grid,
    const std::vector<Atom>& atoms)
{
    const AtomView atom_view{atoms.data(), atoms.size()};
    AtomVoxelizer(0.0).voxelize(incremental_grid, atom_view);
    AtomVoxelizer(0.0).voxelize(full_grid, atom_view);
    ExteriorClassifier{}.classify(incremental_grid);
    ExteriorClassifier{}.classify(full_grid);
}

void test_validation_empty_and_cancellation()
{
    const auto grid_spec = make_grid_spec();
    GasGrid incremental_grid(grid_spec);
    GasGrid full_grid(grid_spec);
    ExteriorClassifier{}.classify(incremental_grid);
    ExteriorClassifier{}.classify(full_grid);

    const auto empty_result = apply_and_compare(
        incremental_grid,
        full_grid,
        {{nullptr, 0}, {nullptr, 0}},
        AccessibilityRepairKind::None);
    REQUIRE(!empty_result.incremental.geometry_changed());

    const auto atom = atoms_for_voxels(incremental_grid, {{3, 3, 3}});
    const AtomView atom_view{atom.data(), atom.size()};
    const auto cancellation_result = apply_and_compare(
        incremental_grid,
        full_grid,
        {atom_view, atom_view},
        AccessibilityRepairKind::None);
    REQUIRE(cancellation_result.incremental.blocker_count_changed_voxel_count
        == 0);

    REQUIRE_THROWS_AS(
        AtomChangeUpdater(0.0).apply_atom_changes(
            incremental_grid,
            {{nullptr, 0}, atom_view}),
        std::underflow_error);

    GasGrid unclassified_grid(grid_spec);
    REQUIRE_THROWS_AS(
        AtomChangeUpdater(0.0).apply_atom_changes(
            unclassified_grid,
            {{nullptr, 0}, {nullptr, 0}}),
        std::invalid_argument);
}

void test_atom_move_is_atomic()
{
    const auto grid_spec = make_grid_spec();
    GasGrid incremental_grid(grid_spec);
    GasGrid full_grid(grid_spec);
    const auto old_atom = atoms_for_voxels(incremental_grid, {{3, 3, 3}});
    initialize_grids(incremental_grid, full_grid, old_atom);
    const auto new_atom = atoms_for_voxels(incremental_grid, {{8, 3, 3}});

    const auto results = apply_and_compare(
        incremental_grid,
        full_grid,
        {{new_atom.data(), new_atom.size()},
         {old_atom.data(), old_atom.size()}},
        AccessibilityRepairKind::Mixed);
    REQUIRE(results.incremental.newly_solid_count == 1);
    REQUIRE(results.incremental.newly_gas_count == 1);
    REQUIRE(incremental_grid.blocker_count({3, 3, 3}) == 0);
    REQUIRE(incremental_grid.gas_state({3, 3, 3})
        == GasState::OutsideAccessible);
    REQUIRE(incremental_grid.blocker_count({8, 3, 3}) == 1);
    REQUIRE(incremental_grid.gas_state({8, 3, 3}) == GasState::Solid);
}

std::vector<VoxelCoord> make_two_cavity_solids(const GridSpec& grid_spec)
{
    const auto is_open = [](const VoxelCoord& voxel_coord) {
        const bool cavity_a = voxel_coord.x >= 1 && voxel_coord.x <= 3
            && voxel_coord.y >= 1 && voxel_coord.y <= 3
            && voxel_coord.z >= 2 && voxel_coord.z <= 4;
        const bool channel_a = voxel_coord.x == 2 && voxel_coord.y == 2
            && voxel_coord.z >= 5;
        const bool cavity_b = voxel_coord.x >= 8 && voxel_coord.x <= 10
            && voxel_coord.y >= 1 && voxel_coord.y <= 3
            && voxel_coord.z >= 2 && voxel_coord.z <= 4;
        const bool channel_b_without_plug = voxel_coord.x == 9
            && voxel_coord.y == 2 && voxel_coord.z >= 6;
        return cavity_a || channel_a || cavity_b || channel_b_without_plug;
    };

    std::vector<VoxelCoord> solids;
    for (std::int64_t z = 0;
         z < static_cast<std::int64_t>(grid_spec.dimensions.z);
         ++z) {
        for (std::int64_t y = 0;
             y < static_cast<std::int64_t>(grid_spec.dimensions.y);
             ++y) {
            for (std::int64_t x = 0;
                 x < static_cast<std::int64_t>(grid_spec.dimensions.x);
                 ++x) {
                const VoxelCoord voxel_coord{x, y, z};
                if (!is_open(voxel_coord)) {
                    solids.push_back(voxel_coord);
                }
            }
        }
    }
    return solids;
}

void test_close_one_cavity_and_open_another()
{
    const auto grid_spec = make_grid_spec(12, 5, 8);
    GasGrid incremental_grid(grid_spec);
    GasGrid full_grid(grid_spec);
    const auto initial_atoms = atoms_for_voxels(
        incremental_grid,
        make_two_cavity_solids(grid_spec));
    initialize_grids(incremental_grid, full_grid, initial_atoms);

    const auto added = atoms_for_voxels(incremental_grid, {{2, 2, 5}});
    const auto removed = atoms_for_voxels(incremental_grid, {{9, 2, 5}});
    const auto results = apply_and_compare(
        incremental_grid,
        full_grid,
        {{added.data(), added.size()}, {removed.data(), removed.size()}},
        AccessibilityRepairKind::Mixed);

    REQUIRE(results.incremental.repair_closed_voxel_count > 0);
    REQUIRE(results.incremental.repair_opened_voxel_count > 0);
    REQUIRE(results.incremental.used_closing_repair());
    REQUIRE(results.incremental.used_opening_repair());
    REQUIRE(incremental_grid.gas_state({2, 2, 3}) == GasState::ClosedVoid);
    REQUIRE(incremental_grid.gas_state({9, 2, 3})
        == GasState::OutsideAccessible);
}

void test_partially_overlapping_atom_move()
{
    auto grid_spec = make_grid_spec(12, 8, 8);
    grid_spec.periodic.x = true;
    GasGrid incremental_grid(grid_spec);
    GasGrid full_grid(grid_spec);
    const auto old_atom = atoms_for_voxels(
        incremental_grid,
        {{5, 4, 4}},
        1.25);
    initialize_grids(incremental_grid, full_grid, old_atom);
    const auto new_atom = atoms_for_voxels(
        incremental_grid,
        {{6, 4, 4}},
        1.25);

    const auto results = apply_and_compare(
        incremental_grid,
        full_grid,
        {{new_atom.data(), new_atom.size()},
         {old_atom.data(), old_atom.size()}},
        AccessibilityRepairKind::Mixed);
    REQUIRE(results.incremental.blocker_count_changed_voxel_count == 10);
    REQUIRE(results.incremental.newly_solid_count == 5);
    REQUIRE(results.incremental.newly_gas_count == 5);
}

void test_channel_swap_excludes_transient_changes()
{
    const auto grid_spec = make_grid_spec();
    GasGrid incremental_grid(grid_spec);
    GasGrid full_grid(grid_spec);
    const VoxelCoord old_channel{5, 3, 4};
    const VoxelCoord new_channel{6, 3, 4};
    std::vector<VoxelCoord> initial_solids;
    for (std::int64_t y = 0; y < 6; ++y) {
        for (std::int64_t x = 0; x < 12; ++x) {
            const VoxelCoord voxel_coord{x, y, 4};
            if (voxel_coord != old_channel) {
                initial_solids.push_back(voxel_coord);
            }
        }
    }
    initialize_grids(
        incremental_grid,
        full_grid,
        atoms_for_voxels(incremental_grid, initial_solids));
    const auto added = atoms_for_voxels(incremental_grid, {old_channel});
    const auto removed = atoms_for_voxels(incremental_grid, {new_channel});

    const auto results = apply_and_compare(
        incremental_grid,
        full_grid,
        {{added.data(), added.size()}, {removed.data(), removed.size()}},
        AccessibilityRepairKind::Mixed);
    REQUIRE(results.incremental.repair_closed_voxel_count > 0);
    REQUIRE(results.incremental.repair_opened_voxel_count > 0);
    REQUIRE(results.incremental.changed_voxel_ids.size() == 2);
    REQUIRE(incremental_grid.gas_state({4, 3, 2})
        == GasState::OutsideAccessible);
}

void test_pure_direction_paths()
{
    const auto grid_spec = make_grid_spec();
    GasGrid incremental_grid(grid_spec);
    GasGrid full_grid(grid_spec);
    ExteriorClassifier{}.classify(incremental_grid);
    ExteriorClassifier{}.classify(full_grid);
    const auto atom = atoms_for_voxels(incremental_grid, {{4, 3, 3}});

    const auto adsorption = apply_and_compare(
        incremental_grid,
        full_grid,
        {{atom.data(), atom.size()}, {nullptr, 0}},
        AccessibilityRepairKind::Closing);
    REQUIRE(adsorption.incremental.newly_solid_count == 1);
    REQUIRE(adsorption.incremental.newly_gas_count == 0);

    const auto desorption = apply_and_compare(
        incremental_grid,
        full_grid,
        {{nullptr, 0}, {atom.data(), atom.size()}},
        AccessibilityRepairKind::Opening);
    REQUIRE(desorption.incremental.newly_solid_count == 0);
    REQUIRE(desorption.incremental.newly_gas_count == 1);
}

void test_owning_event_buffer_and_convenience_wrappers()
{
    auto grid_spec = make_grid_spec(5, 1, 1);
    grid_spec.reservoir_faces = {};
    grid_spec.reservoir_faces.x_low = true;
    GasGrid gas_grid(grid_spec);
    Atom old_atom{gas_grid.voxel_center({2, 0, 0}), 0.0};
    Atom new_atom{gas_grid.voxel_center({1, 0, 0}), 0.0};
    AtomVoxelizer(0.0).voxelize(gas_grid, {&old_atom, 1});
    ExteriorClassifier{}.classify(gas_grid);

    AtomChangeEventBuffer event_buffer;
    event_buffer.reserve(2, 2);
    event_buffer.record_move(old_atom, new_atom);
    old_atom.position = gas_grid.voxel_center({4, 0, 0});
    new_atom.position = gas_grid.voxel_center({4, 0, 0});
    REQUIRE(event_buffer.added_atom_count() == 1);
    REQUIRE(event_buffer.removed_atom_count() == 1);
    REQUIRE(event_buffer.added_atoms().atoms[0].position.x == 1.5);
    REQUIRE(event_buffer.removed_atoms().atoms[0].position.x == 2.5);

    AtomChangeUpdater updater(0.0);
    const auto move_result = updater.apply_atom_changes(
        gas_grid,
        event_buffer.atom_changes());
    REQUIRE(move_result.repair_kind == AccessibilityRepairKind::Mixed);
    REQUIRE(gas_grid.gas_state({1, 0, 0}) == GasState::Solid);
    REQUIRE(gas_grid.gas_state({2, 0, 0}) == GasState::ClosedVoid);

    const Atom moved_atom{gas_grid.voxel_center({1, 0, 0}), 0.0};
    event_buffer.clear();
    REQUIRE(event_buffer.empty());
    event_buffer.record_desorption(moved_atom);
    const auto desorption_result = updater.apply_desorption(
        gas_grid,
        event_buffer.removed_atoms());
    REQUIRE(desorption_result.repair_kind == AccessibilityRepairKind::Opening);
    REQUIRE(gas_grid.gas_state({4, 0, 0})
        == GasState::OutsideAccessible);

    const Atom adsorbed_atom{gas_grid.voxel_center({3, 0, 0}), 0.0};
    event_buffer.clear();
    event_buffer.record_adsorption(adsorbed_atom);
    const auto adsorption_result = updater.apply_adsorption(
        gas_grid,
        event_buffer.added_atoms());
    REQUIRE(adsorption_result.repair_kind == AccessibilityRepairKind::Closing);
    REQUIRE(gas_grid.gas_state({3, 0, 0}) == GasState::Solid);

    REQUIRE_THROWS_AS(
        event_buffer.append({{nullptr, 1}, {nullptr, 0}}),
        std::invalid_argument);
}

void test_deterministic_random_mixed_sequence()
{
    const auto grid_spec = make_grid_spec(8, 6, 6);
    GasGrid incremental_grid(grid_spec);
    GasGrid full_grid(grid_spec);
    std::vector<bool> occupied(
        static_cast<std::size_t>(incremental_grid.voxel_count()),
        false);
    std::mt19937_64 random_generator(418923U);
    for (VoxelId voxel_id = 0;
         voxel_id < incremental_grid.voxel_count();
         ++voxel_id) {
        occupied[static_cast<std::size_t>(voxel_id)] =
            random_generator() % 5U == 0;
    }

    std::vector<VoxelCoord> initial_voxels;
    for (VoxelId voxel_id = 0;
         voxel_id < incremental_grid.voxel_count();
         ++voxel_id) {
        if (occupied[static_cast<std::size_t>(voxel_id)]) {
            initial_voxels.push_back(incremental_grid.voxel_coord(voxel_id));
        }
    }
    initialize_grids(
        incremental_grid,
        full_grid,
        atoms_for_voxels(incremental_grid, initial_voxels));

    AtomChangeUpdater incremental_updater(0.0);
    AtomChangeUpdater full_updater(
        0.0,
        AtomChangeRepairMode::FullReclassification);
    for (std::size_t event_index = 0; event_index < 120; ++event_index) {
        std::vector<VoxelCoord> occupied_voxels;
        std::vector<VoxelCoord> empty_voxels;
        for (VoxelId voxel_id = 0;
             voxel_id < incremental_grid.voxel_count();
             ++voxel_id) {
            auto& destination = occupied[static_cast<std::size_t>(voxel_id)]
                ? occupied_voxels
                : empty_voxels;
            destination.push_back(incremental_grid.voxel_coord(voxel_id));
        }
        std::shuffle(
            occupied_voxels.begin(),
            occupied_voxels.end(),
            random_generator);
        std::shuffle(
            empty_voxels.begin(),
            empty_voxels.end(),
            random_generator);
        const std::size_t removal_count = std::min<std::size_t>(
            2,
            occupied_voxels.size());
        const std::size_t addition_count = std::min<std::size_t>(
            2,
            empty_voxels.size());
        occupied_voxels.resize(removal_count);
        empty_voxels.resize(addition_count);

        const auto added_atoms = atoms_for_voxels(
            incremental_grid,
            empty_voxels);
        const auto removed_atoms = atoms_for_voxels(
            incremental_grid,
            occupied_voxels);
        const AtomChangeBatch changes{
            {added_atoms.data(), added_atoms.size()},
            {removed_atoms.data(), removed_atoms.size()}};
        const auto previous_states = copy_states(incremental_grid);
        const auto incremental_result = incremental_updater.apply_atom_changes(
            incremental_grid,
            changes);
        const auto full_result = full_updater.apply_atom_changes(
            full_grid,
            changes);

        for (const auto& voxel_coord : occupied_voxels) {
            occupied[static_cast<std::size_t>(
                incremental_grid.voxel_id(voxel_coord))] = false;
        }
        for (const auto& voxel_coord : empty_voxels) {
            occupied[static_cast<std::size_t>(
                incremental_grid.voxel_id(voxel_coord))] = true;
        }

        GasGrid reconstructed_grid(grid_spec);
        std::vector<VoxelCoord> final_voxels;
        for (VoxelId voxel_id = 0;
             voxel_id < reconstructed_grid.voxel_count();
             ++voxel_id) {
            if (occupied[static_cast<std::size_t>(voxel_id)]) {
                final_voxels.push_back(reconstructed_grid.voxel_coord(voxel_id));
            }
        }
        const auto final_atoms = atoms_for_voxels(
            reconstructed_grid,
            final_voxels);
        AtomVoxelizer(0.0).voxelize(
            reconstructed_grid,
            {final_atoms.data(), final_atoms.size()});
        ExteriorClassifier{}.classify(reconstructed_grid);

        REQUIRE(incremental_result.newly_solid_count == addition_count);
        REQUIRE(incremental_result.newly_gas_count == removal_count);
        REQUIRE(incremental_result.changed_voxel_ids
            == changed_state_ids(previous_states, reconstructed_grid));
        REQUIRE(full_result.changed_voxel_ids
            == incremental_result.changed_voxel_ids);
        require_same_grid(incremental_grid, reconstructed_grid);
        require_same_grid(full_grid, reconstructed_grid);
    }
}

}  // namespace

int main()
{
    int failure_count = 0;
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"validation, empty, and cancellation", []() {
             test_validation_empty_and_cancellation();
         }},
        {"atom move is atomic", []() { test_atom_move_is_atomic(); }},
        {"close one cavity and open another", []() {
             test_close_one_cavity_and_open_another();
         }},
        {"partially overlapping atom move", []() {
             test_partially_overlapping_atom_move();
         }},
        {"channel swap excludes transient changes", []() {
             test_channel_swap_excludes_transient_changes();
         }},
        {"pure direction paths", []() { test_pure_direction_paths(); }},
        {"owning event buffer and convenience wrappers", []() {
             test_owning_event_buffer_and_convenience_wrappers();
         }},
        {"deterministic random mixed sequence", []() {
             test_deterministic_random_mixed_sequence();
         }}
    };

    for (const auto& test : tests) {
        try {
            test.second();
            std::cout << "[PASS] " << test.first << '\n';
        } catch (const std::exception& error) {
            ++failure_count;
            std::cerr << "[FAIL] " << test.first << ": " << error.what()
                      << '\n';
        }
    }
    if (failure_count == 0) {
        std::cout << tests.size() << " atom-change update test groups passed\n";
    }
    return failure_count == 0 ? 0 : 1;
}
