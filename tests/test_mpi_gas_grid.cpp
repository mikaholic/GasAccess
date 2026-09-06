#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/gas_grid.hpp"
#include "gasaccess/mpi_gas_grid.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using gasaccess::AlignedGridGeometry;
using gasaccess::Atom;
using gasaccess::AtomChangeBatch;
using gasaccess::AtomVoxelizer;
using gasaccess::DistributedGasAccessibilityQuery;
using gasaccess::DistributedGasGrid;
using gasaccess::DistributedVoxelOccupancyChange;
using gasaccess::Face;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridDimensions;
using gasaccess::GridSpec;
using gasaccess::MpiDecompositionSpec;
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
#define REQUIRE_NEAR(actual, expected, tolerance) \
    require_near((actual), (expected), (tolerance), __FILE__, __LINE__)
#define REQUIRE_THROWS_AS(expression, exception_type) \
    require_throws<exception_type>( \
        [&]() { expression; }, #expression, __FILE__, __LINE__)

VoxelCoord process_location(int rank, const GridDimensions& process_grid)
{
    const auto unsigned_rank = static_cast<std::uint64_t>(rank);
    return {
        static_cast<std::int64_t>(unsigned_rank % process_grid.x),
        static_cast<std::int64_t>(
            (unsigned_rank / process_grid.x) % process_grid.y),
        static_cast<std::int64_t>(
            unsigned_rank / (process_grid.x * process_grid.y))
    };
}

GridSpec make_grid_spec(bool periodic_z = false)
{
    GridSpec grid_spec{};
    grid_spec.origin = {-1.0, 2.0, 4.0};
    grid_spec.spacing = {0.5, 0.75, 1.25};
    grid_spec.dimensions = {8, 8, 8};
    grid_spec.periodic = {true, true, periodic_z};
    if (!periodic_z) {
        grid_spec.reservoir_faces.z_high = true;
    }
    return grid_spec;
}

MpiDecompositionSpec make_decomposition_spec(
    const GridSpec& grid_spec,
    const GridDimensions& process_grid,
    int rank)
{
    const auto location = process_location(rank, process_grid);
    const Point3 global_upper{
        grid_spec.origin.x
            + grid_spec.spacing.x * static_cast<double>(grid_spec.dimensions.x),
        grid_spec.origin.y
            + grid_spec.spacing.y * static_cast<double>(grid_spec.dimensions.y),
        grid_spec.origin.z
            + grid_spec.spacing.z * static_cast<double>(grid_spec.dimensions.z)
    };
    const Point3 process_lengths{
        (global_upper.x - grid_spec.origin.x)
            / static_cast<double>(process_grid.x),
        (global_upper.y - grid_spec.origin.y)
            / static_cast<double>(process_grid.y),
        (global_upper.z - grid_spec.origin.z)
            / static_cast<double>(process_grid.z)
    };

    MpiDecompositionSpec decomposition_spec{};
    decomposition_spec.communicator = MPI_COMM_WORLD;
    decomposition_spec.global_lower = grid_spec.origin;
    decomposition_spec.global_upper = global_upper;
    decomposition_spec.local_lower = {
        grid_spec.origin.x + static_cast<double>(location.x) * process_lengths.x,
        grid_spec.origin.y + static_cast<double>(location.y) * process_lengths.y,
        grid_spec.origin.z + static_cast<double>(location.z) * process_lengths.z
    };
    decomposition_spec.local_upper = {
        decomposition_spec.local_lower.x + process_lengths.x,
        decomposition_spec.local_lower.y + process_lengths.y,
        decomposition_spec.local_lower.z + process_lengths.z
    };
    decomposition_spec.process_grid = process_grid;
    decomposition_spec.process_location = location;
    decomposition_spec.atom_ghost_distance = 2.0;
    decomposition_spec.maximum_excluded_radius = 1.0;
    return decomposition_spec;
}

GasState expected_state(const VoxelCoord& coordinate)
{
    const auto value = (coordinate.x
        + 2 * coordinate.y
        + 3 * coordinate.z) % 3;
    if (value == 0) {
        return GasState::Solid;
    }
    if (value == 1) {
        return GasState::OutsideAccessible;
    }
    return GasState::ClosedVoid;
}

std::optional<VoxelCoord> normalize_coordinate(
    VoxelCoord coordinate,
    const GridSpec& grid_spec)
{
    const auto normalize_axis = [](std::int64_t& value,
                                   std::uint64_t dimension,
                                   bool periodic) {
        if (value < 0) {
            if (!periodic) {
                return false;
            }
            value = static_cast<std::int64_t>(dimension - 1);
        } else if (static_cast<std::uint64_t>(value) >= dimension) {
            if (!periodic) {
                return false;
            }
            value = 0;
        }
        return true;
    };
    if (!normalize_axis(coordinate.x, grid_spec.dimensions.x, grid_spec.periodic.x)
        || !normalize_axis(
            coordinate.y,
            grid_spec.dimensions.y,
            grid_spec.periodic.y)
        || !normalize_axis(
            coordinate.z,
            grid_spec.dimensions.z,
            grid_spec.periodic.z)) {
        return std::nullopt;
    }
    return coordinate;
}

bool expected_accessibility(const VoxelCoord& coordinate, const GridSpec& grid_spec)
{
    if (expected_state(coordinate) == GasState::OutsideAccessible) {
        return true;
    }
    const std::array<VoxelCoord, 6> candidates{{
        {coordinate.x - 1, coordinate.y, coordinate.z},
        {coordinate.x + 1, coordinate.y, coordinate.z},
        {coordinate.x, coordinate.y - 1, coordinate.z},
        {coordinate.x, coordinate.y + 1, coordinate.z},
        {coordinate.x, coordinate.y, coordinate.z - 1},
        {coordinate.x, coordinate.y, coordinate.z + 1}
    }};
    for (const auto& candidate : candidates) {
        const auto normalized = normalize_coordinate(candidate, grid_spec);
        if (normalized && *normalized != coordinate
            && expected_state(*normalized) == GasState::OutsideAccessible) {
            return true;
        }
    }
    return false;
}

void fill_expected_states(DistributedGasGrid& grid)
{
    const auto& range = grid.owned_range();
    for (std::int64_t z = range.begin.z; z < range.end.z; ++z) {
        for (std::int64_t y = range.begin.y; y < range.end.y; ++y) {
            for (std::int64_t x = range.begin.x; x < range.end.x; ++x) {
                const VoxelCoord coordinate{x, y, z};
                grid.set_owned_gas_state(coordinate, expected_state(coordinate));
            }
        }
    }
}

void test_alignment_and_validation(int rank, int size)
{
    const GridDimensions process_grid{
        static_cast<std::uint64_t>(size),
        1,
        1
    };
    const AlignedGridGeometry geometry = gasaccess::make_aligned_grid_geometry(
        {-1.0, 2.0, 4.0},
        {9.0, 15.0, 21.0},
        process_grid,
        {0.9, 1.1, 1.4});
    REQUIRE(geometry.dimensions.x % process_grid.x == 0);
    REQUIRE(geometry.dimensions.y % process_grid.y == 0);
    REQUIRE(geometry.dimensions.z % process_grid.z == 0);
    REQUIRE_NEAR(
        geometry.spacing.x * static_cast<double>(geometry.dimensions.x),
        10.0,
        1.0e-12);
    REQUIRE_NEAR(
        geometry.spacing.y * static_cast<double>(geometry.dimensions.y),
        13.0,
        1.0e-12);
    REQUIRE_NEAR(
        geometry.spacing.z * static_cast<double>(geometry.dimensions.z),
        17.0,
        1.0e-12);
    REQUIRE(geometry.spacing.x != geometry.spacing.y);

    gasaccess::validate_atom_ghost_coverage(1.0, 1.0);
    REQUIRE_THROWS_AS(
        gasaccess::validate_atom_ghost_coverage(0.99, 1.0),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        gasaccess::validate_atom_ghost_coverage(
            std::numeric_limits<double>::quiet_NaN(),
            1.0),
        std::invalid_argument);

    auto grid_spec = make_grid_spec();
    auto decomposition_spec = make_decomposition_spec(grid_spec, process_grid, rank);
    decomposition_spec.local_lower.x += 0.1;
    REQUIRE_THROWS_AS(
        DistributedGasGrid grid(grid_spec, decomposition_spec),
        std::invalid_argument);

    decomposition_spec = make_decomposition_spec(grid_spec, process_grid, rank);
    decomposition_spec.atom_ghost_distance = 0.9;
    REQUIRE_THROWS_AS(
        DistributedGasGrid grid(grid_spec, decomposition_spec),
        std::invalid_argument);
}

void verify_halo_layout(
    const GridDimensions& process_grid,
    bool periodic_z,
    int rank)
{
    auto grid_spec = make_grid_spec(periodic_z);
    DistributedGasGrid grid(
        grid_spec,
        make_decomposition_spec(grid_spec, process_grid, rank));
    const DistributedGasAccessibilityQuery query(grid);
    fill_expected_states(grid);
    grid.exchange_ghost_states();

    std::uint64_t global_owned_count = 0;
    const auto local_owned_count = grid.owned_voxel_count();
    MPI_Allreduce(
        &local_owned_count,
        &global_owned_count,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    REQUIRE(global_owned_count == 8U * 8U * 8U);

    const auto& range = grid.owned_range();
    for (std::int64_t z = range.begin.z; z < range.end.z; ++z) {
        for (std::int64_t y = range.begin.y; y < range.end.y; ++y) {
            for (std::int64_t x = range.begin.x; x < range.end.x; ++x) {
                const VoxelCoord coordinate{x, y, z};
                const std::array<VoxelCoord, 6> candidates{{
                    {x - 1, y, z}, {x + 1, y, z},
                    {x, y - 1, z}, {x, y + 1, z},
                    {x, y, z - 1}, {x, y, z + 1}
                }};
                for (const auto& candidate : candidates) {
                    const auto normalized = normalize_coordinate(candidate, grid_spec);
                    if (normalized) {
                        REQUIRE(grid.gas_state(*normalized)
                            == expected_state(*normalized));
                    }
                }
                REQUIRE(query.is_site_accessible(grid.voxel_center(coordinate))
                    == expected_accessibility(coordinate, grid_spec));
            }
        }
    }

    if (!periodic_z && grid.decomposition().spec().process_location.z == 0) {
        REQUIRE(grid.decomposition().neighbor_rank(Face::ZLow) == MPI_PROC_NULL);
    }
    if (!periodic_z
        && static_cast<std::uint64_t>(
            grid.decomposition().spec().process_location.z)
            == process_grid.z - 1) {
        REQUIRE(grid.decomposition().neighbor_rank(Face::ZHigh) == MPI_PROC_NULL);
    }
}

void test_halo_exchange_and_queries(int rank, int size)
{
    const auto process_count = static_cast<std::uint64_t>(size);
    verify_halo_layout({process_count, 1, 1}, false, rank);
    verify_halo_layout({1, process_count, 1}, false, rank);
    verify_halo_layout({1, 1, process_count}, false, rank);
    verify_halo_layout({1, 1, process_count}, true, rank);
    if (size == 4) {
        verify_halo_layout({2, 2, 1}, false, rank);
    }
}

void test_distributed_non_cubic_voxelization(int rank, int size)
{
    const GridDimensions process_grid{
        static_cast<std::uint64_t>(size),
        1,
        1
    };
    auto grid_spec = make_grid_spec();
    std::vector<Atom> valid_atoms{
        {{0.95, 4.9, 7.1}, 0.45},
        {{2.05, 5.0, 8.0}, 0.35},
        {{grid_spec.origin.x - 0.1, 3.0, 6.0}, 0.25}
    };
    const Atom ignored_nonperiodic_image{
        {1.0, 3.0, 14.1},
        0.4
    };
    std::vector<Atom> distributed_atoms = valid_atoms;
    distributed_atoms.push_back(ignored_nonperiodic_image);
    constexpr double precursor_radius = 0.3;

    DistributedGasGrid distributed_grid(
        grid_spec,
        make_decomposition_spec(grid_spec, process_grid, rank));
    const auto local_newly_solid = distributed_grid.voxelize_owned_atoms(
        {distributed_atoms.data(), distributed_atoms.size()},
        precursor_radius);

    GasGrid serial_grid(grid_spec);
    const auto serial_newly_solid = AtomVoxelizer(precursor_radius).voxelize(
        serial_grid,
        {valid_atoms.data(), valid_atoms.size()});
    const auto& range = distributed_grid.owned_range();
    for (std::int64_t z = range.begin.z; z < range.end.z; ++z) {
        for (std::int64_t y = range.begin.y; y < range.end.y; ++y) {
            for (std::int64_t x = range.begin.x; x < range.end.x; ++x) {
                const VoxelCoord coordinate{x, y, z};
                REQUIRE(distributed_grid.gas_state(coordinate)
                    == serial_grid.gas_state(coordinate));
            }
        }
    }

    std::uint64_t global_newly_solid = 0;
    MPI_Allreduce(
        &local_newly_solid,
        &global_newly_solid,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    REQUIRE(global_newly_solid == serial_newly_solid);

    const Atom oversized_atom{{1.0, 3.0, 6.0}, 0.8};
    REQUIRE_THROWS_AS(
        distributed_grid.voxelize_owned_atoms(
            {&oversized_atom, 1},
            precursor_radius),
        std::invalid_argument);
}

void test_owned_blocker_count_invariant(int rank, int size)
{
    const GridDimensions process_grid{
        static_cast<std::uint64_t>(size),
        1,
        1
    };
    const auto grid_spec = make_grid_spec();
    DistributedGasGrid grid(
        grid_spec,
        make_decomposition_spec(grid_spec, process_grid, rank));
    const auto coordinate = grid.owned_range().begin;
    REQUIRE(grid.owned_blocker_count(coordinate) == 0);

    grid.set_owned_blocker_count(coordinate, 2);
    REQUIRE(grid.owned_blocker_count(coordinate) == 2);
    REQUIRE(grid.gas_state(coordinate) == GasState::Solid);
    REQUIRE(grid.owned_gas_state_count(GasState::Solid) == 1);
    grid.set_owned_blocker_count(coordinate, 1);
    REQUIRE(grid.owned_gas_state_count(GasState::Solid) == 1);
    grid.set_owned_blocker_count(coordinate, 0);
    REQUIRE(grid.gas_state(coordinate) == GasState::Unclassified);

    grid.set_owned_gas_state(coordinate, GasState::Solid);
    REQUIRE(grid.owned_blocker_count(coordinate) == 1);
    grid.set_owned_gas_state(coordinate, GasState::ClosedVoid);
    REQUIRE(grid.owned_blocker_count(coordinate) == 0);

    grid.fill_owned_gas_state(GasState::Solid);
    const auto& range = grid.owned_range();
    for (std::int64_t z = range.begin.z; z < range.end.z; ++z) {
        for (std::int64_t y = range.begin.y; y < range.end.y; ++y) {
            for (std::int64_t x = range.begin.x; x < range.end.x; ++x) {
                REQUIRE(grid.owned_blocker_count({x, y, z}) == 1);
            }
        }
    }
    grid.fill_owned_gas_state(GasState::OutsideAccessible);
    REQUIRE(grid.owned_blocker_count(coordinate) == 0);

    if (size > 1) {
        VoxelCoord non_owned = coordinate;
        non_owned.x = rank == 0
            ? range.end.x
            : range.begin.x - 1;
        REQUIRE(!grid.owns(non_owned));
        REQUIRE_THROWS_AS(
            grid.owned_blocker_count(non_owned),
            std::out_of_range);
    }
}

void test_distributed_mixed_change_parity_and_rollback(int rank, int size)
{
    const GridDimensions process_grid{
        static_cast<std::uint64_t>(size),
        1,
        1
    };
    auto grid_spec = make_grid_spec();
    constexpr double precursor_radius = 0.0;
    const Atom first{{0.75, 3.875, 8.375}, 0.0};
    const Atom second{{2.25, 3.875, 8.375}, 0.0};
    const Atom periodic_third{{-1.25, 3.875, 8.375}, 0.0};
    const std::vector<Atom> initial_atoms{first, first, second};
    const std::vector<Atom> additions{first, periodic_third};
    const std::vector<Atom> removals{first, second};

    DistributedGasGrid distributed_grid(
        grid_spec,
        make_decomposition_spec(grid_spec, process_grid, rank));
    distributed_grid.voxelize_owned_atoms(
        {initial_atoms.data(), initial_atoms.size()},
        precursor_radius);

    GasGrid serial_grid(grid_spec);
    const AtomVoxelizer voxelizer(precursor_radius);
    voxelizer.voxelize(
        serial_grid,
        {initial_atoms.data(), initial_atoms.size()});
    const AtomChangeBatch batch{
        {additions.data(), additions.size()},
        {removals.data(), removals.size()}
    };
    const auto serial_result = voxelizer.apply_atom_changes(serial_grid, batch);

    std::vector<DistributedVoxelOccupancyChange> changes;
    const auto local_result = distributed_grid.apply_owned_atom_changes(
        batch,
        precursor_radius,
        changes);
    std::array<std::uint64_t, 3> local_counts{{
        local_result.blocker_count_changed_voxel_count,
        local_result.newly_solid_count,
        local_result.newly_gas_count
    }};
    std::array<std::uint64_t, 3> global_counts{};
    MPI_Allreduce(
        local_counts.data(),
        global_counts.data(),
        static_cast<int>(global_counts.size()),
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    REQUIRE(global_counts[0]
        == serial_result.blocker_count_changed_voxel_count);
    REQUIRE(global_counts[1] == serial_result.newly_solid_count);
    REQUIRE(global_counts[2] == serial_result.newly_gas_count);

    const auto& range = distributed_grid.owned_range();
    for (std::int64_t z = range.begin.z; z < range.end.z; ++z) {
        for (std::int64_t y = range.begin.y; y < range.end.y; ++y) {
            for (std::int64_t x = range.begin.x; x < range.end.x; ++x) {
                const VoxelCoord coordinate{x, y, z};
                REQUIRE(distributed_grid.owned_blocker_count(coordinate)
                    == serial_grid.blocker_count(coordinate));
                REQUIRE(distributed_grid.gas_state(coordinate)
                    == serial_grid.gas_state(coordinate));
            }
        }
    }

    const Atom invalid{{0.75, 3.875, 8.375}, -0.1};
    changes.assign(1, {{-1, -1, -1}, 7, 8, GasState::ClosedVoid});
    REQUIRE_THROWS_AS(
        distributed_grid.apply_owned_atom_changes(
            {{&invalid, 1}, {&periodic_third, 1}},
            precursor_radius,
            changes),
        std::invalid_argument);
    REQUIRE(changes.size() == 1);
    REQUIRE((changes[0].voxel_coord == VoxelCoord{-1, -1, -1}));

    const VoxelCoord absent_coordinate{2, 2, 3};
    if (distributed_grid.owns(absent_coordinate)) {
        const Atom absent{distributed_grid.voxel_center(absent_coordinate), 0.0};
        REQUIRE(distributed_grid.owned_blocker_count(absent_coordinate) == 0);
        REQUIRE_THROWS_AS(
            distributed_grid.apply_owned_atom_changes(
                {{nullptr, 0}, {&absent, 1}},
                precursor_radius,
                changes),
            std::underflow_error);
        REQUIRE(distributed_grid.owned_blocker_count(absent_coordinate) == 0);
        REQUIRE(changes.size() == 1);
    }
}

}  // namespace

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_failure_count = 0;
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"alignment and validation", [&]() {
             test_alignment_and_validation(rank, size);
         }},
        {"halo exchange and queries", [&]() {
             test_halo_exchange_and_queries(rank, size);
         }},
        {"distributed non-cubic voxelization", [&]() {
             test_distributed_non_cubic_voxelization(rank, size);
         }},
        {"owned blocker-count invariant", [&]() {
             test_owned_blocker_count_invariant(rank, size);
         }},
        {"distributed mixed-change parity and rollback", [&]() {
             test_distributed_mixed_change_parity_and_rollback(rank, size);
         }}
    };

    for (const auto& test : tests) {
        try {
            test.second();
            if (rank == 0) {
                std::cout << "[PASS] " << test.first << '\n';
            }
        } catch (const std::exception& error) {
            ++local_failure_count;
            std::cerr << "[rank " << rank << "] [FAIL] " << test.first
                      << ": " << error.what() << '\n';
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int global_failure_count = 0;
    MPI_Allreduce(
        &local_failure_count,
        &global_failure_count,
        1,
        MPI_INT,
        MPI_SUM,
        MPI_COMM_WORLD);
    if (rank == 0 && global_failure_count == 0) {
        std::cout << tests.size() << " MPI test groups passed on "
                  << size << " rank(s)\n";
    }
    MPI_Finalize();
    return global_failure_count == 0 ? 0 : 1;
}
