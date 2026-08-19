#include "gasaccess/affected_region_repair.hpp"
#include "gasaccess/deposition_updater.hpp"
#include "gasaccess/exterior_classifier.hpp"

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

using gasaccess::AffectedRegionRepair;
using gasaccess::Atom;
using gasaccess::ConnectivityRepairMode;
using gasaccess::DepositionUpdateResult;
using gasaccess::DepositionUpdater;
using gasaccess::ExteriorClassifier;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridSpec;
using gasaccess::RemovedVoxel;
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
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {x, y, z};
    return grid_spec;
}

void set_solid(GasGrid& gas_grid, const VoxelCoord& voxel_coord)
{
    gas_grid.set_gas_state(gas_grid.voxel_id(voxel_coord), GasState::Solid);
}

void require_same_result(
    const GasGrid& actual_grid,
    const DepositionUpdateResult& actual_result,
    const GasGrid& reference_grid,
    const DepositionUpdateResult& reference_result)
{
    REQUIRE(actual_result.newly_solid_count == reference_result.newly_solid_count);
    REQUIRE(actual_result.changed_voxel_ids == reference_result.changed_voxel_ids);
    REQUIRE(actual_result.classification.solid_count
        == reference_result.classification.solid_count);
    REQUIRE(actual_result.classification.outside_accessible_count
        == reference_result.classification.outside_accessible_count);
    REQUIRE(actual_result.classification.closed_void_count
        == reference_result.classification.closed_void_count);
    REQUIRE(actual_grid.voxel_count() == reference_grid.voxel_count());
    for (VoxelId voxel_id = 0; voxel_id < actual_grid.voxel_count(); ++voxel_id) {
        REQUIRE(actual_grid.gas_state(voxel_id) == reference_grid.gas_state(voxel_id));
    }
}

void apply_and_compare(
    GasGrid& actual_grid,
    GasGrid& reference_grid,
    const std::vector<Atom>& atoms,
    DepositionUpdateResult& actual_result)
{
    const DepositionUpdater actual_updater(0.0);
    const DepositionUpdater reference_updater(
        0.0,
        ConnectivityRepairMode::FullReclassification);
    actual_result = actual_updater.apply_deposition(
        actual_grid,
        {atoms.data(), atoms.size()});
    const auto reference_result = reference_updater.apply_deposition(
        reference_grid,
        {atoms.data(), atoms.size()});
    REQUIRE(reference_result.used_full_reclassification());
    require_same_result(actual_grid, actual_result, reference_grid, reference_result);
}

void test_input_validation()
{
    auto grid_spec = make_grid_spec(3, 1, 1);
    grid_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid gas_grid(grid_spec);
    ExteriorClassifier{}.classify(gas_grid);
    AffectedRegionRepair repair;

    const auto empty_result = repair.repair(gas_grid, {nullptr, 0});
    REQUIRE(empty_result.visited_voxel_count == 0);
    REQUIRE(empty_result.newly_closed_voxel_ids.empty());
    REQUIRE_THROWS_AS(
        repair.repair(gas_grid, {nullptr, 1}),
        std::invalid_argument);

    RemovedVoxel removed_voxel{
        gas_grid.voxel_id({1, 0, 0}),
        GasState::OutsideAccessible
    };
    REQUIRE_THROWS_AS(
        repair.repair(gas_grid, {&removed_voxel, 1}),
        std::invalid_argument);
}

void test_trench_cavity_repair_is_local()
{
    auto grid_spec = make_grid_spec(9, 1, 9);
    grid_spec.reservoir_faces.z_high = true;
    GasGrid actual_grid(grid_spec);
    GasGrid reference_grid(grid_spec);
    for (std::int64_t z = 0; z < 8; ++z) {
        set_solid(actual_grid, {2, 0, z});
        set_solid(actual_grid, {6, 0, z});
        set_solid(reference_grid, {2, 0, z});
        set_solid(reference_grid, {6, 0, z});
    }
    ExteriorClassifier{}.classify(actual_grid);
    ExteriorClassifier{}.classify(reference_grid);

    const std::vector<Atom> roof_atoms{
        {{3.5, 0.5, 7.5}, 0.0},
        {{4.5, 0.5, 7.5}, 0.0},
        {{5.5, 0.5, 7.5}, 0.0}
    };
    DepositionUpdateResult result{};
    apply_and_compare(actual_grid, reference_grid, roof_atoms, result);

    REQUIRE(result.used_affected_region_repair());
    REQUIRE(!result.used_full_reclassification());
    REQUIRE(result.repair_closed_voxel_count == 21);
    REQUIRE(result.repair_visited_voxel_count == 21);
    REQUIRE(result.changed_voxel_ids.size() == 24);
    REQUIRE(result.repair_visited_voxel_count < actual_grid.voxel_count() / 3);
}

void test_removing_only_source_closes_all_remaining_gas()
{
    auto grid_spec = make_grid_spec(5, 1, 1);
    grid_spec.explicit_source_voxels = {{2, 0, 0}};
    GasGrid actual_grid(grid_spec);
    GasGrid reference_grid(grid_spec);
    ExteriorClassifier{}.classify(actual_grid);
    ExteriorClassifier{}.classify(reference_grid);

    const std::vector<Atom> atoms{{{2.5, 0.5, 0.5}, 0.0}};
    DepositionUpdateResult result{};
    apply_and_compare(actual_grid, reference_grid, atoms, result);
    REQUIRE(result.used_affected_region_repair());
    REQUIRE(result.repair_closed_voxel_count == 4);
    REQUIRE(result.repair_visited_voxel_count == 4);
}

void test_second_source_preserves_connectivity()
{
    auto grid_spec = make_grid_spec(7, 1, 1);
    grid_spec.explicit_source_voxels = {{0, 0, 0}, {6, 0, 0}};
    GasGrid actual_grid(grid_spec);
    GasGrid reference_grid(grid_spec);
    ExteriorClassifier{}.classify(actual_grid);
    ExteriorClassifier{}.classify(reference_grid);

    const std::vector<Atom> atoms{{{0.5, 0.5, 0.5}, 0.0}};
    DepositionUpdateResult result{};
    apply_and_compare(actual_grid, reference_grid, atoms, result);
    REQUIRE(result.used_affected_region_repair());
    REQUIRE(result.repair_closed_voxel_count == 0);
    REQUIRE(result.classification.outside_accessible_count == 6);
}

void test_multi_voxel_cut_closes_only_middle_component()
{
    auto grid_spec = make_grid_spec(11, 1, 1);
    grid_spec.explicit_source_voxels = {{0, 0, 0}, {10, 0, 0}};
    GasGrid actual_grid(grid_spec);
    GasGrid reference_grid(grid_spec);
    ExteriorClassifier{}.classify(actual_grid);
    ExteriorClassifier{}.classify(reference_grid);

    const std::vector<Atom> atoms{
        {{3.5, 0.5, 0.5}, 0.0},
        {{7.5, 0.5, 0.5}, 0.0}
    };
    DepositionUpdateResult result{};
    apply_and_compare(actual_grid, reference_grid, atoms, result);
    REQUIRE(result.used_affected_region_repair());
    REQUIRE(result.repair_closed_voxel_count == 3);
    REQUIRE(actual_grid.gas_state({4, 0, 0}) == GasState::ClosedVoid);
    REQUIRE(actual_grid.gas_state({5, 0, 0}) == GasState::ClosedVoid);
    REQUIRE(actual_grid.gas_state({6, 0, 0}) == GasState::ClosedVoid);
    REQUIRE(actual_grid.gas_state({2, 0, 0}) == GasState::OutsideAccessible);
    REQUIRE(actual_grid.gas_state({8, 0, 0}) == GasState::OutsideAccessible);
}

void test_periodic_seam_component_repair()
{
    auto grid_spec = make_grid_spec(7, 1, 1);
    grid_spec.periodic.x = true;
    grid_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid actual_grid(grid_spec);
    GasGrid reference_grid(grid_spec);
    set_solid(actual_grid, {2, 0, 0});
    set_solid(reference_grid, {2, 0, 0});
    ExteriorClassifier{}.classify(actual_grid);
    ExteriorClassifier{}.classify(reference_grid);

    const std::vector<Atom> atoms{{{6.5, 0.5, 0.5}, 0.0}};
    DepositionUpdateResult result{};
    apply_and_compare(actual_grid, reference_grid, atoms, result);
    REQUIRE(result.used_affected_region_repair());
    REQUIRE(result.repair_closed_voxel_count == 3);
    REQUIRE(actual_grid.gas_state({0, 0, 0}) == GasState::OutsideAccessible);
    for (std::int64_t x = 3; x < 6; ++x) {
        REQUIRE(actual_grid.gas_state({x, 0, 0}) == GasState::ClosedVoid);
    }
}

}  // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"input validation", test_input_validation},
        {"trench cavity repair is local", test_trench_cavity_repair_is_local},
        {"only-source removal", test_removing_only_source_closes_all_remaining_gas},
        {"second source preserves gas", test_second_source_preserves_connectivity},
        {"multi-voxel cut", test_multi_voxel_cut_closes_only_middle_component},
        {"periodic seam repair", test_periodic_seam_component_repair}
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
