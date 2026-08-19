#include "gasaccess/atom_voxelizer.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using gasaccess::Atom;
using gasaccess::AtomView;
using gasaccess::AtomVoxelizer;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridSpec;
using gasaccess::Point3;
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
    std::uint64_t x = 5,
    std::uint64_t y = 5,
    std::uint64_t z = 5)
{
    GridSpec grid_spec{};
    grid_spec.spacing = {1.0, 1.0, 1.0};
    grid_spec.dimensions = {x, y, z};
    return grid_spec;
}

std::set<VoxelId> solid_ids(const GasGrid& gas_grid)
{
    std::set<VoxelId> result;
    for (VoxelId id = 0; id < gas_grid.voxel_count(); ++id) {
        if (gas_grid.gas_state(id) == GasState::Solid) {
            result.insert(id);
        }
    }
    return result;
}

std::set<VoxelId> ids_for(
    const GasGrid& gas_grid,
    const std::vector<VoxelCoord>& voxel_coords)
{
    std::set<VoxelId> result;
    for (const auto& voxel_coord : voxel_coords) {
        result.insert(gas_grid.voxel_id(voxel_coord));
    }
    return result;
}

double reference_axis_distance(
    double lhs,
    double rhs,
    double length,
    bool periodic)
{
    if (!periodic) {
        return std::abs(lhs - rhs);
    }
    return std::abs(std::remainder(lhs - rhs, length));
}

bool reference_is_blocked(
    const GasGrid& gas_grid,
    VoxelId voxel_id,
    const std::vector<Atom>& atoms,
    double precursor_radius)
{
    const auto& grid_spec = gas_grid.grid_spec();
    const Point3 lengths{
        grid_spec.spacing.x * static_cast<double>(grid_spec.dimensions.x),
        grid_spec.spacing.y * static_cast<double>(grid_spec.dimensions.y),
        grid_spec.spacing.z * static_cast<double>(grid_spec.dimensions.z)
    };
    const auto center = gas_grid.voxel_center(voxel_id);

    for (const auto& atom : atoms) {
        const double excluded_radius = atom.radius + precursor_radius;
        const double x_distance = reference_axis_distance(
            center.x,
            atom.position.x,
            lengths.x,
            grid_spec.periodic.x);
        const double y_distance = reference_axis_distance(
            center.y,
            atom.position.y,
            lengths.y,
            grid_spec.periodic.y);
        const double z_distance = reference_axis_distance(
            center.z,
            atom.position.z,
            lengths.z,
            grid_spec.periodic.z);
        const double distance_squared = x_distance * x_distance
            + y_distance * y_distance
            + z_distance * z_distance;
        if (distance_squared <= excluded_radius * excluded_radius) {
            return true;
        }
    }
    return false;
}

void test_single_atom_spherical_exclusion()
{
    GasGrid gas_grid(make_grid_spec());
    const Atom atom{{2.5, 2.5, 2.5}, 0.4};
    const AtomVoxelizer voxelizer(0.6);

    const auto newly_solid = voxelizer.voxelize(gas_grid, {&atom, 1});
    REQUIRE(newly_solid == 7);
    REQUIRE(solid_ids(gas_grid) == ids_for(gas_grid, {
        {2, 2, 2},
        {1, 2, 2}, {3, 2, 2},
        {2, 1, 2}, {2, 3, 2},
        {2, 2, 1}, {2, 2, 3}
    }));

    for (VoxelId id = 0; id < gas_grid.voxel_count(); ++id) {
        if (gas_grid.gas_state(id) != GasState::Solid) {
            REQUIRE(gas_grid.gas_state(id) == GasState::Unclassified);
        }
    }
}

void test_newly_solid_change_capture()
{
    GasGrid gas_grid(make_grid_spec());
    gas_grid.fill_gas_state(GasState::OutsideAccessible);
    const auto center_id = gas_grid.voxel_id({2, 2, 2});
    gas_grid.set_gas_state(center_id, GasState::ClosedVoid);

    const Atom atom{{2.5, 2.5, 2.5}, 0.4};
    const AtomVoxelizer voxelizer(0.6);
    std::vector<RemovedVoxel> removed_voxels{{
        gas_grid.voxel_count(),
        GasState::Unclassified
    }};
    const auto newly_solid = voxelizer.voxelize(
        gas_grid,
        {&atom, 1},
        removed_voxels);

    REQUIRE(newly_solid == 7);
    REQUIRE(removed_voxels.size() == 7);
    std::set<VoxelId> captured_ids;
    for (const auto& removed_voxel : removed_voxels) {
        REQUIRE(captured_ids.insert(removed_voxel.voxel_id).second);
        REQUIRE(gas_grid.gas_state(removed_voxel.voxel_id) == GasState::Solid);
        if (removed_voxel.voxel_id == center_id) {
            REQUIRE(removed_voxel.previous_state == GasState::ClosedVoid);
        } else {
            REQUIRE(removed_voxel.previous_state == GasState::OutsideAccessible);
        }
    }
    REQUIRE(captured_ids == ids_for(gas_grid, {
        {2, 2, 2},
        {1, 2, 2}, {3, 2, 2},
        {2, 1, 2}, {2, 3, 2},
        {2, 2, 1}, {2, 2, 3}
    }));

    REQUIRE(voxelizer.voxelize(
        gas_grid,
        {&atom, 1},
        removed_voxels) == 0);
    REQUIRE(removed_voxels.empty());
}

void test_additive_idempotent_and_order_independent()
{
    const std::vector<Atom> atoms{
        {{2.5, 2.5, 2.5}, 0.8},
        {{3.5, 2.5, 2.5}, 0.5},
        {{1.5, 1.5, 2.5}, 0.3}
    };
    auto reversed_atoms = atoms;
    std::reverse(reversed_atoms.begin(), reversed_atoms.end());
    const AtomVoxelizer voxelizer(0.35);

    GasGrid first_grid(make_grid_spec(6, 5, 5));
    const auto first_count = voxelizer.voxelize(
        first_grid,
        {atoms.data(), atoms.size()});
    REQUIRE(first_count == solid_ids(first_grid).size());
    REQUIRE(voxelizer.voxelize(first_grid, {atoms.data(), atoms.size()}) == 0);

    GasGrid second_grid(make_grid_spec(6, 5, 5));
    const auto second_count = voxelizer.voxelize(
        second_grid,
        {reversed_atoms.data(), reversed_atoms.size()});
    REQUIRE(second_count == first_count);
    REQUIRE(solid_ids(second_grid) == solid_ids(first_grid));

    GasGrid preclassified_grid(make_grid_spec());
    preclassified_grid.fill_gas_state(GasState::OutsideAccessible);
    const Atom center_atom{{2.5, 2.5, 2.5}, 0.1};
    REQUIRE(voxelizer.voxelize(preclassified_grid, {&center_atom, 1}) == 1);
    REQUIRE(preclassified_grid.gas_state(
        preclassified_grid.voxel_id({2, 2, 2})) == GasState::Solid);
    REQUIRE(preclassified_grid.gas_state(
        preclassified_grid.voxel_id({0, 0, 0})) == GasState::OutsideAccessible);
}

void test_periodic_seam_exclusion()
{
    auto periodic_spec = make_grid_spec(4, 3, 3);
    periodic_spec.periodic.x = true;
    GasGrid periodic_grid(periodic_spec);

    const Atom atom{{0.1, 1.5, 1.5}, 0.2};
    const AtomVoxelizer voxelizer(0.5);
    voxelizer.voxelize(periodic_grid, {&atom, 1});
    REQUIRE(periodic_grid.gas_state(
        periodic_grid.voxel_id({3, 1, 1})) == GasState::Solid);

    GasGrid nonperiodic_grid(make_grid_spec(4, 3, 3));
    voxelizer.voxelize(nonperiodic_grid, {&atom, 1});
    REQUIRE(nonperiodic_grid.gas_state(
        nonperiodic_grid.voxel_id({3, 1, 1})) == GasState::Unclassified);

    auto fully_periodic_spec = make_grid_spec(4, 4, 4);
    fully_periodic_spec.periodic = {true, true, true};
    GasGrid fully_periodic_grid(fully_periodic_spec);
    const Atom outside_primary_box{{-0.5, 4.5, 8.5}, 0.0};
    const AtomVoxelizer point_voxelizer(0.0);
    REQUIRE(point_voxelizer.voxelize(
        fully_periodic_grid,
        {&outside_primary_box, 1}) == 1);
    REQUIRE(fully_periodic_grid.gas_state(
        fully_periodic_grid.voxel_id({3, 0, 0})) == GasState::Solid);
}

void test_steric_pore_threshold()
{
    const std::vector<Atom> atoms{
        {{1.0, 0.5, 0.5}, 0.5},
        {{4.0, 0.5, 0.5}, 0.5}
    };
    const auto grid_spec = make_grid_spec(5, 1, 1);

    GasGrid below_threshold_grid(grid_spec);
    AtomVoxelizer below_threshold_voxelizer(0.99);
    below_threshold_voxelizer.voxelize(
        below_threshold_grid,
        {atoms.data(), atoms.size()});
    REQUIRE(below_threshold_grid.gas_state(
        below_threshold_grid.voxel_id({2, 0, 0})) == GasState::Unclassified);

    GasGrid at_threshold_grid(grid_spec);
    AtomVoxelizer at_threshold_voxelizer(1.0);
    at_threshold_voxelizer.voxelize(
        at_threshold_grid,
        {atoms.data(), atoms.size()});
    REQUIRE(at_threshold_grid.gas_state(
        at_threshold_grid.voxel_id({2, 0, 0})) == GasState::Solid);

    GasGrid above_threshold_grid(grid_spec);
    AtomVoxelizer above_threshold_voxelizer(1.01);
    above_threshold_voxelizer.voxelize(
        above_threshold_grid,
        {atoms.data(), atoms.size()});
    REQUIRE(above_threshold_grid.gas_state(
        above_threshold_grid.voxel_id({2, 0, 0})) == GasState::Solid);
}

void test_large_origin_and_small_spacing()
{
    auto grid_spec = make_grid_spec(3, 1, 1);
    grid_spec.origin = {1.0e8, -2.0e8, 5.0e7};
    grid_spec.spacing = {1.0e-3, 1.0e-3, 1.0e-3};
    GasGrid gas_grid(grid_spec);

    const Atom atom{gas_grid.voxel_center({1, 0, 0}), 0.0};
    const AtomVoxelizer voxelizer(4.0e-4);
    REQUIRE(voxelizer.voxelize(gas_grid, {&atom, 1}) == 1);
    REQUIRE(gas_grid.gas_state(gas_grid.voxel_id({1, 0, 0})) == GasState::Solid);
    REQUIRE(gas_grid.gas_state(
        gas_grid.voxel_id({0, 0, 0})) == GasState::Unclassified);
    REQUIRE(gas_grid.gas_state(
        gas_grid.voxel_id({2, 0, 0})) == GasState::Unclassified);
}

void test_invalid_voxelization_inputs()
{
    REQUIRE_THROWS_AS(AtomVoxelizer voxelizer(-0.1), std::invalid_argument);
    REQUIRE_THROWS_AS(
        AtomVoxelizer voxelizer(std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        AtomVoxelizer voxelizer(std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);

    GasGrid gas_grid(make_grid_spec());
    const AtomVoxelizer voxelizer(0.2);
    REQUIRE(voxelizer.voxelize(gas_grid, {nullptr, 0}) == 0);
    REQUIRE_THROWS_AS(voxelizer.voxelize(gas_grid, {nullptr, 1}), std::invalid_argument);

    Atom atom{{2.5, 2.5, 2.5}, -0.1};
    REQUIRE_THROWS_AS(voxelizer.voxelize(gas_grid, {&atom, 1}), std::invalid_argument);

    atom = {{2.5, 2.5, 2.5}, std::numeric_limits<double>::infinity()};
    REQUIRE_THROWS_AS(voxelizer.voxelize(gas_grid, {&atom, 1}), std::invalid_argument);

    atom = {{std::numeric_limits<double>::quiet_NaN(), 2.5, 2.5}, 0.1};
    REQUIRE_THROWS_AS(voxelizer.voxelize(gas_grid, {&atom, 1}), std::invalid_argument);

    atom = {{5.0, 2.5, 2.5}, 0.1};
    REQUIRE_THROWS_AS(voxelizer.voxelize(gas_grid, {&atom, 1}), std::invalid_argument);

    atom = {{2.5, 2.5, 2.5}, std::numeric_limits<double>::max()};
    REQUIRE_THROWS_AS(voxelizer.voxelize(gas_grid, {&atom, 1}), std::invalid_argument);

    GasGrid unchanged_grid(make_grid_spec());
    const std::vector<Atom> partially_invalid_atoms{
        {{2.5, 2.5, 2.5}, 0.1},
        {{2.5, 2.5, 2.5}, -0.1}
    };
    REQUIRE_THROWS_AS(
        voxelizer.voxelize(
            unchanged_grid,
            {partially_invalid_atoms.data(), partially_invalid_atoms.size()}),
        std::invalid_argument);
    REQUIRE(solid_ids(unchanged_grid).empty());
}

void test_randomized_against_brute_force()
{
    std::mt19937_64 random_engine(0x6a09e667f3bcc909ULL);
    std::uniform_real_distribution<double> unit_distribution(0.0, 1.0);
    std::uniform_real_distribution<double> radius_distribution(0.0, 0.55);
    constexpr double precursor_radius = 0.3;

    for (unsigned int mask = 0; mask < 8U; ++mask) {
        for (int sample = 0; sample < 12; ++sample) {
            GridSpec grid_spec{};
            grid_spec.origin = {-1.0, 0.25, 5.0};
            grid_spec.spacing = {
                0.31 + 0.01 * static_cast<double>(sample % 3),
                0.47 + 0.02 * static_cast<double>(sample % 4),
                0.73 + 0.03 * static_cast<double>(sample % 5)
            };
            grid_spec.dimensions = {5, 4, 3};
            grid_spec.periodic = {
                (mask & 1U) != 0U,
                (mask & 2U) != 0U,
                (mask & 4U) != 0U
            };

            const Point3 lengths{
                grid_spec.spacing.x * static_cast<double>(grid_spec.dimensions.x),
                grid_spec.spacing.y * static_cast<double>(grid_spec.dimensions.y),
                grid_spec.spacing.z * static_cast<double>(grid_spec.dimensions.z)
            };

            std::vector<Atom> atoms;
            atoms.reserve(6);
            for (int atom_index = 0; atom_index < 6; ++atom_index) {
                Point3 position{
                    grid_spec.origin.x + unit_distribution(random_engine) * lengths.x,
                    grid_spec.origin.y + unit_distribution(random_engine) * lengths.y,
                    grid_spec.origin.z + unit_distribution(random_engine) * lengths.z
                };

                if (grid_spec.periodic.x && atom_index == 0) {
                    position.x += lengths.x;
                }
                if (grid_spec.periodic.y && atom_index == 1) {
                    position.y -= lengths.y;
                }
                if (grid_spec.periodic.z && atom_index == 2) {
                    position.z += 2.0 * lengths.z;
                }
                atoms.push_back({position, radius_distribution(random_engine)});
            }

            GasGrid gas_grid(grid_spec);
            const AtomVoxelizer voxelizer(precursor_radius);
            const auto newly_solid = voxelizer.voxelize(
                gas_grid,
                {atoms.data(), atoms.size()});

            VoxelId expected_solid_count = 0;
            for (VoxelId id = 0; id < gas_grid.voxel_count(); ++id) {
                const bool expected_solid = reference_is_blocked(
                    gas_grid,
                    id,
                    atoms,
                    precursor_radius);
                REQUIRE((gas_grid.gas_state(id) == GasState::Solid) == expected_solid);
                if (expected_solid) {
                    ++expected_solid_count;
                }
            }
            REQUIRE(newly_solid == expected_solid_count);
        }
    }
}

}  // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"single-atom spherical exclusion", test_single_atom_spherical_exclusion},
        {"newly solid change capture", test_newly_solid_change_capture},
        {"additive and order-independent voxelization",
         test_additive_idempotent_and_order_independent},
        {"periodic seam exclusion", test_periodic_seam_exclusion},
        {"steric pore threshold", test_steric_pore_threshold},
        {"large origin and small spacing", test_large_origin_and_small_spacing},
        {"invalid voxelization inputs", test_invalid_voxelization_inputs},
        {"randomized brute-force differential", test_randomized_against_brute_force}
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
