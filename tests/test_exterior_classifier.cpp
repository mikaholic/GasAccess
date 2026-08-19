#include "gasaccess/exterior_classifier.hpp"

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

using gasaccess::ClassificationSummary;
using gasaccess::ExteriorClassifier;
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

#define REQUIRE(condition) \
    require_condition((condition), #condition, __FILE__, __LINE__)

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

void set_solid(GasGrid& gas_grid, const std::vector<VoxelCoord>& voxel_coords)
{
    for (const auto& voxel_coord : voxel_coords) {
        gas_grid.set_gas_state(gas_grid.voxel_id(voxel_coord), GasState::Solid);
    }
}

void require_summary_consistent(
    const GasGrid& gas_grid,
    const ClassificationSummary& summary)
{
    VoxelId solid_count = 0;
    VoxelId outside_count = 0;
    VoxelId closed_count = 0;

    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        switch (gas_grid.gas_state(voxel_id)) {
        case GasState::Solid:
            ++solid_count;
            break;
        case GasState::OutsideAccessible:
            ++outside_count;
            break;
        case GasState::ClosedVoid:
            ++closed_count;
            break;
        case GasState::Unclassified:
            REQUIRE(false);
            break;
        }
    }

    REQUIRE(summary.solid_count == solid_count);
    REQUIRE(summary.outside_accessible_count == outside_count);
    REQUIRE(summary.closed_void_count == closed_count);
    REQUIRE(solid_count + outside_count + closed_count == gas_grid.voxel_count());
}

std::vector<bool> reference_source_reachability(const GasGrid& gas_grid)
{
    std::vector<bool> reached(
        static_cast<std::size_t>(gas_grid.voxel_count()),
        false);
    std::vector<VoxelId> frontier;

    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        if (gas_grid.gas_state(voxel_id) != GasState::Solid
            && gas_grid.is_reservoir_source(voxel_id)) {
            reached[static_cast<std::size_t>(voxel_id)] = true;
            frontier.push_back(voxel_id);
        }
    }

    std::size_t frontier_index = 0;
    while (frontier_index < frontier.size()) {
        const auto voxel_id = frontier[frontier_index];
        ++frontier_index;
        const auto neighbors = gas_grid.neighbors(voxel_id);
        for (std::size_t index = 0; index < neighbors.count; ++index) {
            const auto neighbor_id = neighbors.ids[index];
            const auto local_index = static_cast<std::size_t>(neighbor_id);
            if (reached[local_index]
                || gas_grid.gas_state(neighbor_id) == GasState::Solid) {
                continue;
            }
            reached[local_index] = true;
            frontier.push_back(neighbor_id);
        }
    }
    return reached;
}

void require_connectivity_invariants(const GasGrid& gas_grid)
{
    const auto reached = reference_source_reachability(gas_grid);
    for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
        const auto state = gas_grid.gas_state(voxel_id);
        if (state == GasState::Solid) {
            REQUIRE(!reached[static_cast<std::size_t>(voxel_id)]);
            continue;
        }
        REQUIRE((state == GasState::OutsideAccessible)
            == reached[static_cast<std::size_t>(voxel_id)]);
    }
}

void test_empty_solid_and_sourceless_domains()
{
    ExteriorClassifier classifier;

    auto empty_spec = make_grid_spec(3, 2, 2);
    empty_spec.reservoir_faces.z_high = true;
    GasGrid empty_grid(empty_spec);
    const auto empty_summary = classifier.classify(empty_grid);
    REQUIRE(empty_summary.solid_count == 0);
    REQUIRE(empty_summary.outside_accessible_count == 12);
    REQUIRE(empty_summary.closed_void_count == 0);
    require_summary_consistent(empty_grid, empty_summary);
    require_connectivity_invariants(empty_grid);

    GasGrid solid_grid(empty_spec);
    solid_grid.fill_gas_state(GasState::Solid);
    const auto solid_summary = classifier.classify(solid_grid);
    REQUIRE(solid_summary.solid_count == 12);
    REQUIRE(solid_summary.outside_accessible_count == 0);
    REQUIRE(solid_summary.closed_void_count == 0);
    require_summary_consistent(solid_grid, solid_summary);

    GasGrid sourceless_grid(make_grid_spec(3, 2, 2));
    sourceless_grid.fill_gas_state(GasState::OutsideAccessible);
    const auto sourceless_summary = classifier.classify(sourceless_grid);
    REQUIRE(sourceless_summary.solid_count == 0);
    REQUIRE(sourceless_summary.outside_accessible_count == 0);
    REQUIRE(sourceless_summary.closed_void_count == 12);
    require_summary_consistent(sourceless_grid, sourceless_summary);
    require_connectivity_invariants(sourceless_grid);
}

void test_sources_and_reclassification()
{
    ExteriorClassifier classifier;

    auto grid_spec = make_grid_spec();
    grid_spec.explicit_source_voxels = {{1, 1, 1}};
    GasGrid blocked_source_grid(grid_spec);
    blocked_source_grid.set_gas_state(
        blocked_source_grid.voxel_id({1, 1, 1}),
        GasState::Solid);
    const auto blocked_summary = classifier.classify(blocked_source_grid);
    REQUIRE(blocked_summary.solid_count == 1);
    REQUIRE(blocked_summary.outside_accessible_count == 0);
    REQUIRE(blocked_summary.closed_void_count == 26);
    require_connectivity_invariants(blocked_source_grid);

    grid_spec = make_grid_spec(1, 1, 1);
    grid_spec.reservoir_faces = {true, true, true, true, true, true};
    GasGrid repeated_face_grid(grid_spec);
    const auto repeated_face_summary = classifier.classify(repeated_face_grid);
    REQUIRE(repeated_face_summary.outside_accessible_count == 1);
    REQUIRE(repeated_face_summary.closed_void_count == 0);

    grid_spec = make_grid_spec(4, 3, 2);
    grid_spec.periodic = {true, true, true};
    grid_spec.explicit_source_voxels = {{2, 1, 1}};
    GasGrid fully_periodic_grid(grid_spec);
    fully_periodic_grid.fill_gas_state(GasState::ClosedVoid);
    const auto first_summary = classifier.classify(fully_periodic_grid);
    REQUIRE(first_summary.outside_accessible_count == 24);
    fully_periodic_grid.set_gas_state(
        fully_periodic_grid.voxel_id({0, 0, 0}),
        GasState::Solid);
    const auto second_summary = classifier.classify(fully_periodic_grid);
    REQUIRE(second_summary.solid_count == 1);
    REQUIRE(second_summary.outside_accessible_count == 23);
    REQUIRE(second_summary.closed_void_count == 0);
    require_summary_consistent(fully_periodic_grid, second_summary);
    require_connectivity_invariants(fully_periodic_grid);
}

void test_open_and_sealed_trench()
{
    auto grid_spec = make_grid_spec(5, 1, 5);
    grid_spec.reservoir_faces.z_high = true;
    ExteriorClassifier classifier;

    const std::vector<VoxelCoord> sidewalls{
        {1, 0, 0}, {3, 0, 0},
        {1, 0, 1}, {3, 0, 1},
        {1, 0, 2}, {3, 0, 2},
        {1, 0, 3}, {3, 0, 3}
    };

    GasGrid open_trench(grid_spec);
    set_solid(open_trench, sidewalls);
    const auto open_summary = classifier.classify(open_trench);
    REQUIRE(open_summary.solid_count == 8);
    REQUIRE(open_summary.outside_accessible_count == 17);
    REQUIRE(open_summary.closed_void_count == 0);
    for (std::int64_t z = 0; z < 5; ++z) {
        REQUIRE(open_trench.gas_state(
            open_trench.voxel_id({2, 0, z})) == GasState::OutsideAccessible);
    }
    require_connectivity_invariants(open_trench);

    GasGrid sealed_trench(grid_spec);
    set_solid(sealed_trench, sidewalls);
    set_solid(sealed_trench, {{2, 0, 3}});
    const auto sealed_summary = classifier.classify(sealed_trench);
    REQUIRE(sealed_summary.solid_count == 9);
    REQUIRE(sealed_summary.outside_accessible_count == 13);
    REQUIRE(sealed_summary.closed_void_count == 3);
    for (std::int64_t z = 0; z < 3; ++z) {
        REQUIRE(sealed_trench.gas_state(
            sealed_trench.voxel_id({2, 0, z})) == GasState::ClosedVoid);
    }
    REQUIRE(sealed_trench.gas_state(
        sealed_trench.voxel_id({2, 0, 4})) == GasState::OutsideAccessible);
    require_summary_consistent(sealed_trench, sealed_summary);
    require_connectivity_invariants(sealed_trench);
}

void test_enclosed_cavity_and_face_channel()
{
    auto grid_spec = make_grid_spec();
    grid_spec.reservoir_faces.x_low = true;
    ExteriorClassifier classifier;

    const std::vector<VoxelCoord> shell{
        {0, 1, 1}, {2, 1, 1},
        {1, 0, 1}, {1, 2, 1},
        {1, 1, 0}, {1, 1, 2}
    };

    GasGrid enclosed_grid(grid_spec);
    set_solid(enclosed_grid, shell);
    const auto enclosed_summary = classifier.classify(enclosed_grid);
    REQUIRE(enclosed_summary.closed_void_count == 1);
    REQUIRE(enclosed_grid.gas_state(
        enclosed_grid.voxel_id({1, 1, 1})) == GasState::ClosedVoid);
    require_connectivity_invariants(enclosed_grid);

    GasGrid channel_grid(grid_spec);
    set_solid(channel_grid, {
        {2, 1, 1},
        {1, 0, 1}, {1, 2, 1},
        {1, 1, 0}, {1, 1, 2}
    });
    const auto channel_summary = classifier.classify(channel_grid);
    REQUIRE(channel_summary.closed_void_count == 0);
    REQUIRE(channel_grid.gas_state(
        channel_grid.voxel_id({1, 1, 1})) == GasState::OutsideAccessible);
    require_connectivity_invariants(channel_grid);
}

void test_diagonal_contacts_do_not_connect()
{
    ExteriorClassifier classifier;

    auto edge_spec = make_grid_spec(3, 3, 1);
    edge_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid edge_grid(edge_spec);
    edge_grid.fill_gas_state(GasState::Solid);
    edge_grid.set_gas_state(edge_grid.voxel_id({0, 0, 0}), GasState::Unclassified);
    edge_grid.set_gas_state(edge_grid.voxel_id({1, 1, 0}), GasState::Unclassified);
    const auto edge_summary = classifier.classify(edge_grid);
    REQUIRE(edge_summary.outside_accessible_count == 1);
    REQUIRE(edge_summary.closed_void_count == 1);
    REQUIRE(edge_grid.gas_state(
        edge_grid.voxel_id({1, 1, 0})) == GasState::ClosedVoid);

    auto corner_spec = make_grid_spec();
    corner_spec.explicit_source_voxels = {{0, 0, 0}};
    GasGrid corner_grid(corner_spec);
    corner_grid.fill_gas_state(GasState::Solid);
    corner_grid.set_gas_state(
        corner_grid.voxel_id({0, 0, 0}),
        GasState::Unclassified);
    corner_grid.set_gas_state(
        corner_grid.voxel_id({1, 1, 1}),
        GasState::Unclassified);
    const auto corner_summary = classifier.classify(corner_grid);
    REQUIRE(corner_summary.outside_accessible_count == 1);
    REQUIRE(corner_summary.closed_void_count == 1);
    REQUIRE(corner_grid.gas_state(
        corner_grid.voxel_id({1, 1, 1})) == GasState::ClosedVoid);
}

void require_periodic_seam_path(
    const GridSpec& grid_spec,
    const VoxelCoord& source,
    const VoxelCoord& blocked_1,
    const VoxelCoord& blocked_2,
    const VoxelCoord& target)
{
    GasGrid gas_grid(grid_spec);
    set_solid(gas_grid, {blocked_1, blocked_2});
    const auto summary = ExteriorClassifier{}.classify(gas_grid);
    REQUIRE(summary.solid_count == 2);
    REQUIRE(summary.outside_accessible_count == 2);
    REQUIRE(summary.closed_void_count == 0);
    REQUIRE(gas_grid.is_reservoir_source(source));
    REQUIRE(gas_grid.gas_state(
        gas_grid.voxel_id(target)) == GasState::OutsideAccessible);
    require_connectivity_invariants(gas_grid);
}

void test_periodic_paths_across_each_axis()
{
    auto x_spec = make_grid_spec(4, 1, 1);
    x_spec.periodic.x = true;
    x_spec.explicit_source_voxels = {{0, 0, 0}};
    require_periodic_seam_path(
        x_spec,
        {0, 0, 0},
        {1, 0, 0},
        {2, 0, 0},
        {3, 0, 0});

    auto y_spec = make_grid_spec(1, 4, 1);
    y_spec.periodic.y = true;
    y_spec.explicit_source_voxels = {{0, 0, 0}};
    require_periodic_seam_path(
        y_spec,
        {0, 0, 0},
        {0, 1, 0},
        {0, 2, 0},
        {0, 3, 0});

    auto z_spec = make_grid_spec(1, 1, 4);
    z_spec.periodic.z = true;
    z_spec.explicit_source_voxels = {{0, 0, 0}};
    require_periodic_seam_path(
        z_spec,
        {0, 0, 0},
        {0, 0, 1},
        {0, 0, 2},
        {0, 0, 3});
}

void test_randomized_classification_invariants()
{
    std::mt19937_64 random_engine(0xbb67ae8584caa73bULL);
    std::bernoulli_distribution solid_distribution(0.38);
    ExteriorClassifier classifier;

    for (unsigned int periodic_mask = 0; periodic_mask < 8U; ++periodic_mask) {
        for (int sample = 0; sample < 16; ++sample) {
            auto grid_spec = make_grid_spec(5, 4, 3);
            grid_spec.periodic = {
                (periodic_mask & 1U) != 0U,
                (periodic_mask & 2U) != 0U,
                (periodic_mask & 4U) != 0U
            };
            if (!grid_spec.periodic.x) {
                grid_spec.reservoir_faces.x_low = true;
            }
            if (!grid_spec.periodic.z) {
                grid_spec.reservoir_faces.z_high = true;
            }
            grid_spec.explicit_source_voxels = {{2, 2, 1}};

            GasGrid gas_grid(grid_spec);
            for (VoxelId voxel_id = 0; voxel_id < gas_grid.voxel_count(); ++voxel_id) {
                if (solid_distribution(random_engine)) {
                    gas_grid.set_gas_state(voxel_id, GasState::Solid);
                } else if ((voxel_id % 2U) == 0U) {
                    gas_grid.set_gas_state(voxel_id, GasState::OutsideAccessible);
                } else {
                    gas_grid.set_gas_state(voxel_id, GasState::ClosedVoid);
                }
            }

            const auto summary = classifier.classify(gas_grid);
            require_summary_consistent(gas_grid, summary);
            require_connectivity_invariants(gas_grid);
        }
    }
}

}  // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"empty, solid, and sourceless domains",
         test_empty_solid_and_sourceless_domains},
        {"sources and reclassification", test_sources_and_reclassification},
        {"open and sealed trench", test_open_and_sealed_trench},
        {"enclosed cavity and face channel", test_enclosed_cavity_and_face_channel},
        {"diagonal contacts remain disconnected", test_diagonal_contacts_do_not_connect},
        {"periodic paths across every axis", test_periodic_paths_across_each_axis},
        {"randomized classification invariants", test_randomized_classification_invariants}
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
