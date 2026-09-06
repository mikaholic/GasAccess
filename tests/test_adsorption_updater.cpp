#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/adsorption_updater.hpp"

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
using gasaccess::ConnectivityRepairMode;
using gasaccess::AdsorptionUpdateResult;
using gasaccess::AdsorptionUpdater;
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
    std::uint64_t x = 5,
    std::uint64_t y = 3,
    std::uint64_t z = 5)
{
    GridSpec grid_spec{};
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {x, y, z};
    return grid_spec;
}

void set_solid(GasGrid& gas_grid, const std::vector<VoxelCoord>& voxel_coords)
{
    for (const auto& voxel_coord : voxel_coords) {
        gas_grid.set_gas_state(gas_grid.voxel_id(voxel_coord), GasState::Solid);
    }
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
    REQUIRE(lhs.voxel_count() == rhs.voxel_count());
    for (VoxelId voxel_id = 0; voxel_id < lhs.voxel_count(); ++voxel_id) {
        REQUIRE(lhs.gas_state(voxel_id) == rhs.gas_state(voxel_id));
    }
}

bool summaries_equal(
    const ClassificationSummary& lhs,
    const ClassificationSummary& rhs)
{
    return lhs.solid_count == rhs.solid_count
        && lhs.outside_accessible_count == rhs.outside_accessible_count
        && lhs.closed_void_count == rhs.closed_void_count;
}

AdsorptionUpdateResult apply_reference_update(
    GasGrid& gas_grid,
    const AtomVoxelizer& atom_voxelizer,
    AtomView atom_view)
{
    const auto previous_states = copy_states(gas_grid);
    AdsorptionUpdateResult result{};
    result.newly_solid_count = atom_voxelizer.voxelize(gas_grid, atom_view);
    if (result.newly_solid_count == 0) {
        for (const auto gas_state : previous_states) {
            if (gas_state == GasState::Solid) {
                ++result.classification.solid_count;
            } else if (gas_state == GasState::OutsideAccessible) {
                ++result.classification.outside_accessible_count;
            } else if (gas_state == GasState::ClosedVoid) {
                ++result.classification.closed_void_count;
            }
        }
        return result;
    }

    result.classification = ExteriorClassifier{}.classify(gas_grid);
    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        if (previous_states[static_cast<std::size_t>(voxel_id)]
            != gas_grid.gas_state(voxel_id)) {
            result.changed_voxel_ids.push_back(voxel_id);
        }
    }
    return result;
}

void require_valid_changed_ids(const AdsorptionUpdateResult& result)
{
    REQUIRE(std::is_sorted(
        result.changed_voxel_ids.begin(),
        result.changed_voxel_ids.end()));
    REQUIRE(std::adjacent_find(
        result.changed_voxel_ids.begin(),
        result.changed_voxel_ids.end()) == result.changed_voxel_ids.end());
}

void test_requires_classified_grid_and_valid_batch()
{
    GasGrid gas_grid(make_grid_spec());
    const AdsorptionUpdater updater(0.0);
    const Atom atom{{2.5, 1.5, 2.5}, 0.0};

    REQUIRE_THROWS_AS(
        updater.apply_adsorption(gas_grid, {&atom, 1}),
        std::invalid_argument);
    REQUIRE(gas_grid.gas_state(gas_grid.voxel_id({2, 1, 2}))
        == GasState::Unclassified);

    auto grid_spec = make_grid_spec();
    grid_spec.reservoir_faces.z_high = true;
    GasGrid classified_grid(grid_spec);
    ExteriorClassifier{}.classify(classified_grid);
    const auto states_before = copy_states(classified_grid);
    const std::vector<Atom> invalid_batch{
        atom,
        {{2.5, 1.5, 2.5}, -0.1}
    };
    REQUIRE_THROWS_AS(
        updater.apply_adsorption(
            classified_grid,
            {invalid_batch.data(), invalid_batch.size()}),
        std::invalid_argument);
    REQUIRE(copy_states(classified_grid) == states_before);
}

void test_no_op_and_overlapping_adsorption()
{
    auto grid_spec = make_grid_spec();
    grid_spec.reservoir_faces.z_high = true;
    GasGrid gas_grid(grid_spec);
    const Atom atom{{2.5, 1.5, 2.5}, 0.4};
    const AtomVoxelizer atom_voxelizer(0.2);
    atom_voxelizer.voxelize(gas_grid, {&atom, 1});
    const auto initial_summary = ExteriorClassifier{}.classify(gas_grid);
    const auto states_before = copy_states(gas_grid);

    const AdsorptionUpdater updater(0.2);
    const auto repeated_result = updater.apply_adsorption(gas_grid, {&atom, 1});
    REQUIRE(!repeated_result.geometry_changed());
    REQUIRE(!repeated_result.used_full_reclassification());
    REQUIRE(repeated_result.newly_solid_count == 0);
    REQUIRE(repeated_result.changed_voxel_ids.empty());
    REQUIRE(summaries_equal(repeated_result.classification, initial_summary));
    REQUIRE(copy_states(gas_grid) == states_before);

    const auto empty_result = updater.apply_adsorption(gas_grid, {nullptr, 0});
    REQUIRE(!empty_result.geometry_changed());
    REQUIRE(!empty_result.used_full_reclassification());
    REQUIRE(empty_result.changed_voxel_ids.empty());
    REQUIRE(summaries_equal(empty_result.classification, initial_summary));
    REQUIRE(copy_states(gas_grid) == states_before);
}

void test_locally_safe_adsorption_skips_full_reclassification()
{
    auto grid_spec = make_grid_spec(5, 5, 5);
    grid_spec.reservoir_faces.z_high = true;
    GasGrid actual_grid(grid_spec);
    GasGrid reference_grid(grid_spec);
    ExteriorClassifier{}.classify(actual_grid);
    ExteriorClassifier{}.classify(reference_grid);

    const Atom atom{{2.5, 2.5, 2.5}, 0.0};
    const auto actual_result = AdsorptionUpdater(0.0).apply_adsorption(
        actual_grid,
        {&atom, 1});
    const auto reference_result = apply_reference_update(
        reference_grid,
        AtomVoxelizer(0.0),
        {&atom, 1});

    REQUIRE(actual_result.geometry_changed());
    REQUIRE(!actual_result.used_full_reclassification());
    REQUIRE(actual_result.changed_voxel_ids == std::vector<VoxelId>({
        actual_grid.voxel_id({2, 2, 2})
    }));
    REQUIRE(summaries_equal(
        actual_result.classification,
        reference_result.classification));
    require_same_states(actual_grid, reference_grid);
}

void test_trench_pinch_off()
{
    auto grid_spec = make_grid_spec(5, 1, 5);
    grid_spec.reservoir_faces.z_high = true;
    GasGrid gas_grid(grid_spec);
    ExteriorClassifier{}.classify(gas_grid);

    const AdsorptionUpdater updater(0.0);
    for (std::int64_t z = 0; z < 4; ++z) {
        const std::vector<Atom> sidewall_atoms{
            {{1.5, 0.5, static_cast<double>(z) + 0.5}, 0.0},
            {{3.5, 0.5, static_cast<double>(z) + 0.5}, 0.0}
        };
        const auto sidewall_result = updater.apply_adsorption(
            gas_grid,
            {sidewall_atoms.data(), sidewall_atoms.size()});
        REQUIRE(sidewall_result.newly_solid_count == 2);
        REQUIRE(sidewall_result.used_affected_region_repair());
        REQUIRE(!sidewall_result.used_full_reclassification());
        REQUIRE(GasAccessibilityQuery(gas_grid).is_site_accessible({2.5, 0.5, 1.5}));
    }

    const GasAccessibilityQuery before_query(gas_grid);
    REQUIRE(before_query.is_site_accessible({2.5, 0.5, 1.5}));

    const Atom roof_atom{{2.5, 0.5, 3.5}, 0.0};
    const auto result = updater.apply_adsorption(gas_grid, {&roof_atom, 1});

    REQUIRE(result.geometry_changed());
    REQUIRE(result.used_affected_region_repair());
    REQUIRE(!result.used_full_reclassification());
    REQUIRE(result.repair_closed_voxel_count == 3);
    REQUIRE(result.repair_visited_voxel_count == 3);
    REQUIRE(result.newly_solid_count == 1);
    REQUIRE(result.classification.solid_count == 9);
    REQUIRE(result.classification.outside_accessible_count == 13);
    REQUIRE(result.classification.closed_void_count == 3);
    REQUIRE(result.changed_voxel_ids == std::vector<VoxelId>({
        gas_grid.voxel_id({2, 0, 0}),
        gas_grid.voxel_id({2, 0, 1}),
        gas_grid.voxel_id({2, 0, 2}),
        gas_grid.voxel_id({2, 0, 3})
    }));
    require_valid_changed_ids(result);

    const GasAccessibilityQuery after_query(gas_grid);
    REQUIRE(!after_query.is_site_accessible({2.5, 0.5, 1.5}));
}

void test_closure_across_periodic_seam()
{
    auto grid_spec = make_grid_spec(5, 1, 1);
    grid_spec.periodic.x = true;
    grid_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid gas_grid(grid_spec);
    set_solid(gas_grid, {{1, 0, 0}});
    ExteriorClassifier{}.classify(gas_grid);
    REQUIRE(gas_grid.gas_state({2, 0, 0}) == GasState::OutsideAccessible);

    const Atom seam_atom{{4.5, 0.5, 0.5}, 0.0};
    const auto result = AdsorptionUpdater(0.0).apply_adsorption(
        gas_grid,
        {&seam_atom, 1});

    REQUIRE(result.newly_solid_count == 1);
    REQUIRE(result.used_affected_region_repair());
    REQUIRE(!result.used_full_reclassification());
    REQUIRE(result.repair_closed_voxel_count == 2);
    REQUIRE(result.classification.solid_count == 2);
    REQUIRE(result.classification.outside_accessible_count == 1);
    REQUIRE(result.classification.closed_void_count == 2);
    REQUIRE(result.changed_voxel_ids == std::vector<VoxelId>({
        gas_grid.voxel_id({2, 0, 0}),
        gas_grid.voxel_id({3, 0, 0}),
        gas_grid.voxel_id({4, 0, 0})
    }));
}

void test_adsorption_blocks_only_explicit_source()
{
    auto grid_spec = make_grid_spec(3, 1, 1);
    grid_spec.periodic.x = true;
    grid_spec.explicit_source_voxels = {{1, 0, 0}};
    GasGrid gas_grid(grid_spec);
    ExteriorClassifier{}.classify(gas_grid);

    const Atom source_atom{{1.5, 0.5, 0.5}, 0.0};
    const auto result = AdsorptionUpdater(0.0).apply_adsorption(
        gas_grid,
        {&source_atom, 1});

    REQUIRE(result.newly_solid_count == 1);
    REQUIRE(result.used_affected_region_repair());
    REQUIRE(!result.used_full_reclassification());
    REQUIRE(result.repair_closed_voxel_count == 2);
    REQUIRE(result.changed_voxel_ids.size() == 3);
    REQUIRE(result.classification.solid_count == 1);
    REQUIRE(result.classification.outside_accessible_count == 0);
    REQUIRE(result.classification.closed_void_count == 2);
    REQUIRE(gas_grid.gas_state({1, 0, 0}) == GasState::Solid);
    REQUIRE(gas_grid.gas_state({0, 0, 0}) == GasState::ClosedVoid);
    REQUIRE(gas_grid.gas_state({2, 0, 0}) == GasState::ClosedVoid);
}

void test_random_sequences_match_explicit_reference()
{
    std::mt19937_64 random_engine(0x3c6ef372fe94f82bULL);
    std::uniform_int_distribution<int> x_distribution(0, 4);
    std::uniform_int_distribution<int> y_distribution(0, 3);
    std::uniform_int_distribution<int> z_distribution(0, 2);
    std::uniform_int_distribution<int> batch_distribution(1, 2);
    constexpr double precursor_radius = 0.15;
    std::size_t local_safe_count = 0;
    std::size_t affected_region_repair_count = 0;

    for (unsigned int periodic_mask = 0; periodic_mask < 8U; ++periodic_mask) {
        auto grid_spec = make_grid_spec(5, 4, 3);
        grid_spec.periodic = {
            (periodic_mask & 1U) != 0U,
            (periodic_mask & 2U) != 0U,
            (periodic_mask & 4U) != 0U
        };
        grid_spec.explicit_source_voxels = {{0, 0, 0}};

        GasGrid actual_grid(grid_spec);
        GasGrid reference_grid(grid_spec);
        ExteriorClassifier{}.classify(actual_grid);
        ExteriorClassifier{}.classify(reference_grid);
        const AdsorptionUpdater updater(precursor_radius);
        const AtomVoxelizer reference_voxelizer(precursor_radius);

        for (int event = 0; event < 20; ++event) {
            std::vector<Atom> atoms;
            const int batch_size = batch_distribution(random_engine);
            atoms.reserve(static_cast<std::size_t>(batch_size));
            for (int atom_index = 0; atom_index < batch_size; ++atom_index) {
                atoms.push_back({{
                    static_cast<double>(x_distribution(random_engine)) + 0.5,
                    static_cast<double>(y_distribution(random_engine)) + 0.5,
                    static_cast<double>(z_distribution(random_engine)) + 0.5
                }, 0.0});
            }

            const AtomView atom_view{atoms.data(), atoms.size()};
            const auto actual_result = updater.apply_adsorption(actual_grid, atom_view);
            const auto reference_result = apply_reference_update(
                reference_grid,
                reference_voxelizer,
                atom_view);

            REQUIRE(actual_result.newly_solid_count
                == reference_result.newly_solid_count);
            REQUIRE(actual_result.changed_voxel_ids
                == reference_result.changed_voxel_ids);
            REQUIRE(summaries_equal(
                actual_result.classification,
                reference_result.classification));
            require_valid_changed_ids(actual_result);
            require_same_states(actual_grid, reference_grid);
            if (actual_result.geometry_changed()) {
                if (actual_result.used_full_reclassification()) {
                    throw TestFailure("default updater unexpectedly used full classification");
                } else if (actual_result.used_affected_region_repair()) {
                    ++affected_region_repair_count;
                } else {
                    ++local_safe_count;
                }
            }
        }
    }
    REQUIRE(local_safe_count != 0);
    REQUIRE(affected_region_repair_count != 0);
}

void test_forced_reference_mode()
{
    auto grid_spec = make_grid_spec(5, 5, 5);
    grid_spec.reservoir_faces.z_high = true;
    GasGrid gas_grid(grid_spec);
    ExteriorClassifier{}.classify(gas_grid);
    const AdsorptionUpdater updater(
        0.0,
        ConnectivityRepairMode::FullReclassification);
    REQUIRE(updater.repair_mode() == ConnectivityRepairMode::FullReclassification);

    const Atom atom{{2.5, 2.5, 2.5}, 0.0};
    const auto result = updater.apply_adsorption(gas_grid, {&atom, 1});
    REQUIRE(result.used_full_reclassification());
    REQUIRE(!result.used_affected_region_repair());
    REQUIRE(result.repair_visited_voxel_count == 0);
    REQUIRE(result.repair_closed_voxel_count == 0);

    REQUIRE_THROWS_AS(
        AdsorptionUpdater(0.0, static_cast<ConnectivityRepairMode>(99)),
        std::invalid_argument);
}

}  // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"classified-grid and batch validation",
         test_requires_classified_grid_and_valid_batch},
        {"no-op and overlapping adsorption", test_no_op_and_overlapping_adsorption},
        {"locally safe adsorption skips full reclassification",
         test_locally_safe_adsorption_skips_full_reclassification},
        {"trench pinch-off", test_trench_pinch_off},
        {"periodic-seam closure", test_closure_across_periodic_seam},
        {"blocking the only explicit source", test_adsorption_blocks_only_explicit_source},
        {"random sequences match explicit reference",
         test_random_sequences_match_explicit_reference},
        {"forced full-reference mode", test_forced_reference_mode}
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
