#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/deposition_updater.hpp"
#include "gasaccess/exterior_classifier.hpp"

#include <cstddef>
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
using gasaccess::GasAccessibilityQuery;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridSpec;
using gasaccess::Point3;

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

#define REQUIRE(condition) \
    require_condition((condition), #condition, __FILE__, __LINE__)

struct KmcSite {
    Point3 atom_position{};
    bool gas_reaction_enabled = false;
};

std::size_t evaluate_gas_reactions(
    const GasAccessibilityQuery& query,
    std::vector<KmcSite>& sites)
{
    std::size_t enabled_count = 0;
    for (auto& site : sites) {
        // This is the only GasAccess operation required in the KMC site loop.
        const bool accessible = query.is_site_accessible(site.atom_position);
        site.gas_reaction_enabled = accessible;
        if (accessible) {
            ++enabled_count;
        }
    }
    return enabled_count;
}

GridSpec make_trench_spec()
{
    GridSpec grid_spec{};
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {5, 1, 5};
    grid_spec.reservoir_faces.z_high = true;
    return grid_spec;
}

void test_single_call_contract_observes_pinch_off()
{
    GasGrid gas_grid(make_trench_spec());
    std::vector<Atom> initial_atoms;
    initial_atoms.reserve(9);
    for (std::int64_t z = 0; z < 4; ++z) {
        initial_atoms.push_back({{1.5, 0.5, static_cast<double>(z) + 0.5}, 0.0});
        initial_atoms.push_back({{3.5, 0.5, static_cast<double>(z) + 0.5}, 0.0});
    }
    initial_atoms.push_back({{2.5, 0.5, 0.5}, 0.0});

    REQUIRE(AtomVoxelizer(0.0).voxelize(
        gas_grid,
        {initial_atoms.data(), initial_atoms.size()}) == initial_atoms.size());
    ExteriorClassifier{}.classify(gas_grid);

    // The floor atom touches the open trench before closure. The second site
    // remains in the exterior gas and acts as a control.
    std::vector<KmcSite> sites{
        {{2.5, 0.5, 0.5}, false},
        {{0.5, 0.5, 3.5}, false}
    };
    const GasAccessibilityQuery query(gas_grid);

    REQUIRE(evaluate_gas_reactions(query, sites) == 2);
    REQUIRE(sites[0].gas_reaction_enabled);
    REQUIRE(sites[1].gas_reaction_enabled);

    const Atom roof_atom{{2.5, 0.5, 3.5}, 0.0};
    const auto update_result = DepositionUpdater(0.0).apply_deposition(
        gas_grid,
        {&roof_atom, 1});
    REQUIRE(update_result.geometry_changed());
    REQUIRE(update_result.repair_closed_voxel_count == 2);

    // Reuse the same query. It references the grid rather than a state copy.
    REQUIRE(evaluate_gas_reactions(query, sites) == 1);
    REQUIRE(!sites[0].gas_reaction_enabled);
    REQUIRE(sites[1].gas_reaction_enabled);
}

void test_query_uses_current_position_without_registration()
{
    GridSpec grid_spec{};
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {5, 3, 3};
    GasGrid gas_grid(grid_spec);
    gas_grid.fill_gas_state(GasState::Solid);
    gas_grid.set_gas_state(
        gas_grid.voxel_id({3, 1, 1}),
        GasState::OutsideAccessible);

    const GasAccessibilityQuery query(gas_grid);
    std::vector<KmcSite> sites{{{2.5, 1.5, 1.5}, false}};
    REQUIRE(evaluate_gas_reactions(query, sites) == 1);

    // Simulate KMC supplying a new position after MD relaxation. GasAccess has
    // no site registration or position cache to update.
    sites[0].atom_position = {0.5, 1.5, 1.5};
    REQUIRE(evaluate_gas_reactions(query, sites) == 0);
    REQUIRE(!sites[0].gas_reaction_enabled);
}

void test_single_call_contract_obeys_periodicity()
{
    GridSpec periodic_spec{};
    periodic_spec.spacing = {1.0, 1.0, 1.0};
    periodic_spec.dimensions = {4, 3, 3};
    periodic_spec.periodic.x = true;
    GasGrid periodic_grid(periodic_spec);
    periodic_grid.fill_gas_state(GasState::Solid);
    periodic_grid.set_gas_state(
        periodic_grid.voxel_id({3, 1, 1}),
        GasState::OutsideAccessible);

    std::vector<KmcSite> sites{{{0.5, 1.5, 1.5}, false}};
    const GasAccessibilityQuery periodic_query(periodic_grid);
    REQUIRE(evaluate_gas_reactions(periodic_query, sites) == 1);

    sites[0].atom_position = {4.5, 1.5, 1.5};
    REQUIRE(evaluate_gas_reactions(periodic_query, sites) == 1);

    periodic_spec.periodic.x = false;
    GasGrid nonperiodic_grid(periodic_spec);
    nonperiodic_grid.fill_gas_state(GasState::Solid);
    nonperiodic_grid.set_gas_state(
        nonperiodic_grid.voxel_id({3, 1, 1}),
        GasState::OutsideAccessible);
    const GasAccessibilityQuery nonperiodic_query(nonperiodic_grid);
    REQUIRE(evaluate_gas_reactions(nonperiodic_query, sites) == 0);
}

}  // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"single call observes pinch-off",
         test_single_call_contract_observes_pinch_off},
        {"query uses current position",
         test_query_uses_current_position_without_registration},
        {"single call obeys periodicity",
         test_single_call_contract_obeys_periodicity}
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
