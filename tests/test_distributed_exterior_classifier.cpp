#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/distributed_exterior_classifier.hpp"
#include "gasaccess/exterior_classifier.hpp"
#include "gasaccess/gas_grid.hpp"
#include "gasaccess/mpi_gas_grid.hpp"

#include <mpi.h>

#include <array>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using gasaccess::DistributedClassificationSummary;
using gasaccess::DistributedExteriorClassifier;
using gasaccess::DistributedGasAccessibilityQuery;
using gasaccess::DistributedGasGrid;
using gasaccess::ExteriorClassifier;
using gasaccess::GasAccessibilityQuery;
using gasaccess::GasGrid;
using gasaccess::GasState;
using gasaccess::GridDimensions;
using gasaccess::GridSpec;
using gasaccess::MpiDecompositionSpec;
using gasaccess::Point3;
using gasaccess::VoxelCoord;

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

GridSpec make_grid_spec()
{
    GridSpec grid_spec{};
    grid_spec.origin = {-1.0, 2.0, 4.0};
    grid_spec.spacing = {0.5, 0.75, 1.25};
    grid_spec.dimensions = {8, 8, 8};
    grid_spec.periodic = {true, true, false};
    grid_spec.reservoir_faces.z_high = true;
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

std::vector<GridDimensions> test_layouts(int size)
{
    const auto process_count = static_cast<std::uint64_t>(size);
    std::vector<GridDimensions> layouts{
        {process_count, 1, 1},
        {1, process_count, 1},
        {1, 1, process_count}};
    if (size == 4) {
        layouts.push_back({2, 2, 1});
    }
    return layouts;
}

std::vector<VoxelCoord> make_barrier(bool leave_opening)
{
    std::vector<VoxelCoord> solid_voxels;
    for (std::int64_t y = 0; y < 8; ++y) {
        for (std::int64_t x = 0; x < 8; ++x) {
            if (leave_opening && x == 3 && y == 3) {
                continue;
            }
            solid_voxels.push_back({x, y, 4});
        }
    }
    return solid_voxels;
}

std::vector<VoxelCoord> make_deterministic_solids()
{
    std::vector<VoxelCoord> solid_voxels;
    for (std::int64_t z = 0; z < 8; ++z) {
        for (std::int64_t y = 0; y < 8; ++y) {
            for (std::int64_t x = 0; x < 8; ++x) {
                const auto pattern = (17 * x + 11 * y + 5 * z + x * y) % 13;
                if (pattern == 0 || pattern == 3 || pattern == 8) {
                    solid_voxels.push_back({x, y, z});
                }
            }
        }
    }
    return solid_voxels;
}

DistributedClassificationSummary verify_against_serial(
    const GridSpec& grid_spec,
    const GridDimensions& process_grid,
    const std::vector<VoxelCoord>& solid_voxels,
    int rank)
{
    GasGrid serial_grid(grid_spec);
    DistributedGasGrid distributed_grid(
        grid_spec,
        make_decomposition_spec(grid_spec, process_grid, rank));

    for (const auto& solid_voxel : solid_voxels) {
        serial_grid.set_gas_state(
            serial_grid.voxel_id(solid_voxel),
            GasState::Solid);
        if (distributed_grid.owns(solid_voxel)) {
            distributed_grid.set_owned_gas_state(
                solid_voxel,
                GasState::Solid);
        }
    }

    const auto serial_summary = ExteriorClassifier().classify(serial_grid);
    DistributedExteriorClassifier classifier;
    const auto distributed_summary = classifier.classify(distributed_grid);

    std::uint64_t counted_solid = 0;
    std::uint64_t counted_outside = 0;
    std::uint64_t counted_closed = 0;
    const auto& range = distributed_grid.owned_range();
    const DistributedGasAccessibilityQuery distributed_query(distributed_grid);
    const GasAccessibilityQuery serial_query(serial_grid);
    for (auto z = range.begin.z; z < range.end.z; ++z) {
        for (auto y = range.begin.y; y < range.end.y; ++y) {
            for (auto x = range.begin.x; x < range.end.x; ++x) {
                const VoxelCoord voxel_coord{x, y, z};
                const auto state = distributed_grid.gas_state(voxel_coord);
                REQUIRE(state == serial_grid.gas_state(voxel_coord));
                if (state == GasState::Solid) {
                    ++counted_solid;
                } else if (state == GasState::OutsideAccessible) {
                    ++counted_outside;
                } else if (state == GasState::ClosedVoid) {
                    ++counted_closed;
                } else {
                    REQUIRE(false);
                }

                const std::array<VoxelCoord, 6> neighbors{{
                    {x - 1, y, z},
                    {x + 1, y, z},
                    {x, y - 1, z},
                    {x, y + 1, z},
                    {x, y, z - 1},
                    {x, y, z + 1}}};
                for (const auto& neighbor : neighbors) {
                    const auto normalized = normalize_coordinate(neighbor, grid_spec);
                    if (normalized.has_value()) {
                        REQUIRE(distributed_grid.gas_state(*normalized)
                            == serial_grid.gas_state(*normalized));
                    }
                }

                const auto center = distributed_grid.voxel_center(voxel_coord);
                REQUIRE(distributed_query.is_site_accessible(center)
                    == serial_query.is_site_accessible(center));
            }
        }
    }

    REQUIRE(distributed_summary.local_solid_count == counted_solid);
    REQUIRE(distributed_summary.local_outside_accessible_count == counted_outside);
    REQUIRE(distributed_summary.local_closed_void_count == counted_closed);
    REQUIRE(distributed_summary.local_visited_voxel_count == counted_outside);
    REQUIRE(counted_solid + counted_outside + counted_closed
        == distributed_grid.owned_voxel_count());

    std::uint64_t global_solid = 0;
    std::uint64_t global_outside = 0;
    std::uint64_t global_closed = 0;
    std::uint64_t global_sent = 0;
    std::uint64_t global_received = 0;
    MPI_Allreduce(
        &counted_solid,
        &global_solid,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        &counted_outside,
        &global_outside,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        &counted_closed,
        &global_closed,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        &distributed_summary.sent_frontier_entry_count,
        &global_sent,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        &distributed_summary.received_frontier_entry_count,
        &global_received,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    REQUIRE(global_solid == serial_summary.solid_count);
    REQUIRE(global_outside == serial_summary.outside_accessible_count);
    REQUIRE(global_closed == serial_summary.closed_void_count);
    REQUIRE(global_sent == global_received);

    return distributed_summary;
}

void test_boundary_sources_and_inactive_ranks(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    for (const auto& layout : test_layouts(size)) {
        const auto summary = verify_against_serial(grid_spec, layout, {}, rank);
        REQUIRE(summary.local_closed_void_count == 0);
        REQUIRE(summary.communication_round_count > 0);
    }

    auto no_source_spec = make_grid_spec();
    no_source_spec.reservoir_faces = {};
    for (const auto& layout : test_layouts(size)) {
        const auto summary = verify_against_serial(
            no_source_spec,
            layout,
            {},
            rank);
        REQUIRE(summary.local_outside_accessible_count == 0);
        REQUIRE(summary.communication_round_count == 0);
    }
}

void test_sealed_and_open_regions(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    for (const auto& layout : test_layouts(size)) {
        const auto sealed_summary = verify_against_serial(
            grid_spec,
            layout,
            make_barrier(false),
            rank);
        std::uint64_t global_closed = 0;
        MPI_Allreduce(
            &sealed_summary.local_closed_void_count,
            &global_closed,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            MPI_COMM_WORLD);
        REQUIRE(global_closed == 8U * 8U * 4U);

        const auto open_summary = verify_against_serial(
            grid_spec,
            layout,
            make_barrier(true),
            rank);
        REQUIRE(open_summary.local_closed_void_count == 0);
    }
}

void test_periodic_explicit_source(int rank, int size)
{
    auto grid_spec = make_grid_spec();
    grid_spec.periodic = {true, true, true};
    grid_spec.reservoir_faces = {};
    grid_spec.explicit_source_voxels = {{1, 0, 0}};
    for (const auto& layout : test_layouts(size)) {
        const auto summary = verify_against_serial(
            grid_spec,
            layout,
            make_deterministic_solids(),
            rank);
        REQUIRE(summary.communication_round_count > 0);
    }
}

void test_deterministic_obstacles(int rank, int size)
{
    const auto grid_spec = make_grid_spec();
    for (const auto& layout : test_layouts(size)) {
        verify_against_serial(
            grid_spec,
            layout,
            make_deterministic_solids(),
            rank);
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
        {"boundary sources and inactive ranks", [&]() {
             test_boundary_sources_and_inactive_ranks(rank, size);
         }},
        {"sealed and open regions", [&]() {
             test_sealed_and_open_regions(rank, size);
         }},
        {"periodic explicit source", [&]() {
             test_periodic_explicit_source(rank, size);
         }},
        {"deterministic obstacles", [&]() {
             test_deterministic_obstacles(rank, size);
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
        std::cout << tests.size() << " distributed classification test groups "
                  << "passed on " << size << " rank(s)\n";
    }
    MPI_Finalize();
    return global_failure_count == 0 ? 0 : 1;
}
