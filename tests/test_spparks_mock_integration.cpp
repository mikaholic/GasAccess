#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/distributed_exterior_classifier.hpp"
#include "gasaccess/exterior_classifier.hpp"
#include "gasaccess/spparks_adapter.hpp"
#include "mock_spparks.hpp"
#include "mock_spparks_fixture.hpp"

#include <mpi.h>

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
using gasaccess::AtomView;
using gasaccess::AtomVoxelizer;
using gasaccess::DistributedExteriorClassifier;
using gasaccess::DistributedGasAccessibilityQuery;
using gasaccess::DistributedGasGrid;
using gasaccess::ExteriorClassifier;
using gasaccess::GasAccessibilityQuery;
using gasaccess::GasGrid;
using gasaccess::GridDimensions;
using gasaccess::SpparksAdapterConfig;
using gasaccess::SpparksAtomBuffer;
using gasaccess::VoxelCoord;
using gasaccess::testing::MockScenario;
using gasaccess::testing::MockSpparksApp;
using gasaccess::testing::MockSpparksDomain;

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

std::vector<GridDimensions> test_layouts(int process_count)
{
    const auto count = static_cast<std::uint64_t>(process_count);
    std::vector<GridDimensions> layouts{
        {count, 1, 1},
        {1, count, 1},
        {1, 1, count}};
    if (process_count == 4) {
        layouts.push_back({2, 2, 1});
    }
    return layouts;
}

std::vector<Atom> make_serial_atoms(
    MockScenario scenario,
    const gasaccess::testing::MockFixtureSpec& fixture_spec)
{
    const auto mock_atoms = gasaccess::testing::make_mock_atoms(
        scenario, fixture_spec);
    std::vector<Atom> atoms;
    atoms.reserve(mock_atoms.size());
    for (const auto& atom : mock_atoms) {
        atoms.push_back({atom.position, atom.radius});
    }
    return atoms;
}

void compare_owned_states(
    const DistributedGasGrid& distributed_grid,
    const GasGrid& serial_grid)
{
    const auto& range = distributed_grid.owned_range();
    for (auto z = range.begin.z; z < range.end.z; ++z) {
        for (auto y = range.begin.y; y < range.end.y; ++y) {
            for (auto x = range.begin.x; x < range.end.x; ++x) {
                const VoxelCoord coordinate{x, y, z};
                REQUIRE(distributed_grid.gas_state(coordinate)
                        == serial_grid.gas_state(coordinate));
            }
        }
    }
}

void compare_face_ghost_states(
    const DistributedGasGrid& distributed_grid,
    const GasGrid& serial_grid)
{
    const auto& range = distributed_grid.owned_range();
    const auto& spec = distributed_grid.global_grid_spec();
    const auto compare = [&](VoxelCoord ghost_coordinate) {
        VoxelCoord serial_coordinate = ghost_coordinate;
        const auto normalize = [](std::int64_t& coordinate,
                                  std::uint64_t dimension,
                                  bool periodic) {
            if (coordinate < 0) {
                if (!periodic) {
                    return false;
                }
                coordinate = static_cast<std::int64_t>(dimension - 1);
            } else if (static_cast<std::uint64_t>(coordinate) >= dimension) {
                if (!periodic) {
                    return false;
                }
                coordinate = 0;
            }
            return true;
        };
        if (!normalize(serial_coordinate.x, spec.dimensions.x, spec.periodic.x)
            || !normalize(
                serial_coordinate.y, spec.dimensions.y, spec.periodic.y)
            || !normalize(
                serial_coordinate.z, spec.dimensions.z, spec.periodic.z)) {
            return;
        }
        REQUIRE(distributed_grid.gas_state(serial_coordinate)
                == serial_grid.gas_state(serial_coordinate));
    };

    for (auto z = range.begin.z; z < range.end.z; ++z) {
        for (auto y = range.begin.y; y < range.end.y; ++y) {
            compare({range.begin.x - 1, y, z});
            compare({range.end.x, y, z});
        }
    }
    for (auto z = range.begin.z; z < range.end.z; ++z) {
        for (auto x = range.begin.x; x < range.end.x; ++x) {
            compare({x, range.begin.y - 1, z});
            compare({x, range.end.y, z});
        }
    }
    for (auto y = range.begin.y; y < range.end.y; ++y) {
        for (auto x = range.begin.x; x < range.end.x; ++x) {
            compare({x, y, range.begin.z - 1});
            compare({x, y, range.end.z});
        }
    }
}

void verify_scenario(
    MockScenario scenario,
    const GridDimensions& process_grid,
    int process_count)
{
    const auto fixture_spec =
        gasaccess::testing::make_mock_fixture_spec(scenario);
    MockSpparksDomain domain(
        MPI_COMM_WORLD,
        fixture_spec.global_lower,
        fixture_spec.global_upper,
        process_grid,
        fixture_spec.spparks_periodic);
    MockSpparksApp app(domain);
    app.set_owned_atoms(gasaccess::testing::make_mock_atoms(
        scenario, fixture_spec, &domain));
    app.synchronize_ghost_atoms(fixture_spec.ghost_distance);

    SpparksAdapterConfig adapter_config{};
    adapter_config.requested_spacing = fixture_spec.requested_spacing;
    adapter_config.gas_periodic = fixture_spec.gas_periodic;
    adapter_config.reservoir_faces = fixture_spec.reservoir_faces;
    adapter_config.atom_ghost_distance = fixture_spec.ghost_distance;
    adapter_config.maximum_excluded_radius =
        fixture_spec.atom_radius + fixture_spec.precursor_radius;

    const auto grid_spec = gasaccess::make_spparks_grid_spec(
        domain, adapter_config);
    REQUIRE(grid_spec.dimensions.x == fixture_spec.dimensions.x);
    REQUIRE(grid_spec.dimensions.y == fixture_spec.dimensions.y);
    REQUIRE(grid_spec.dimensions.z == fixture_spec.dimensions.z);
    REQUIRE(grid_spec.periodic.x && grid_spec.periodic.y);
    REQUIRE(!grid_spec.periodic.z);
    REQUIRE(domain.zperiodic == 1);

    const auto decomposition_spec = gasaccess::make_spparks_decomposition_spec(
        domain, domain.world, adapter_config);
    REQUIRE(decomposition_spec.process_location.x == domain.myloc[0]);
    REQUIRE(decomposition_spec.process_location.y == domain.myloc[1]);
    REQUIRE(decomposition_spec.process_location.z == domain.myloc[2]);
    REQUIRE(decomposition_spec.local_lower.x == domain.subxlo);
    REQUIRE(decomposition_spec.local_upper.z == domain.subzhi);

    DistributedGasGrid distributed_grid(grid_spec, decomposition_spec);
    SpparksAtomBuffer atom_buffer;
    atom_buffer.assign(
        app,
        [&app](std::size_t atom_index) {
            return app.radius[atom_index];
        });
    distributed_grid.voxelize_owned_atoms(
        atom_buffer.atom_view(), fixture_spec.precursor_radius);
    DistributedExteriorClassifier().classify(distributed_grid);

    GasGrid serial_grid(grid_spec);
    const auto serial_atoms = make_serial_atoms(scenario, fixture_spec);
    AtomVoxelizer(fixture_spec.precursor_radius).voxelize(
        serial_grid, AtomView{serial_atoms.data(), serial_atoms.size()});
    ExteriorClassifier().classify(serial_grid);

    compare_owned_states(distributed_grid, serial_grid);
    compare_face_ghost_states(distributed_grid, serial_grid);

    const DistributedGasAccessibilityQuery distributed_query(distributed_grid);
    const GasAccessibilityQuery serial_query(serial_grid);
    std::uint64_t local_inner_probe_count = 0;
    std::uint64_t local_accessible_inner_probe_count = 0;
    for (int atom_index = 0; atom_index < app.nlocal; ++atom_index) {
        const gasaccess::Point3 atom_position{
            app.xyz[atom_index][0],
            app.xyz[atom_index][1],
            app.xyz[atom_index][2]};
        const bool distributed_accessible =
            distributed_query.is_site_accessible(atom_position);
        REQUIRE(distributed_accessible
                == serial_query.is_site_accessible(atom_position));
        const auto coordinate = serial_grid.locate_voxel(atom_position);
        REQUIRE(coordinate.has_value());
        if (gasaccess::testing::is_inner_trench_probe(
                scenario, coordinate.value())) {
            ++local_inner_probe_count;
            if (distributed_accessible) {
                ++local_accessible_inner_probe_count;
            }
        }
    }

    std::uint64_t global_inner_probe_count = 0;
    std::uint64_t global_accessible_inner_probe_count = 0;
    const auto local_ghost_count = static_cast<std::uint64_t>(app.nghost);
    std::uint64_t global_ghost_count = 0;
    MPI_Allreduce(
        &local_inner_probe_count,
        &global_inner_probe_count,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        &local_accessible_inner_probe_count,
        &global_accessible_inner_probe_count,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    MPI_Allreduce(
        &local_ghost_count,
        &global_ghost_count,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD);
    REQUIRE(global_inner_probe_count == 32U);
    REQUIRE(global_accessible_inner_probe_count
            == (scenario == MockScenario::OpenTrench ? 32U : 0U));
    if (process_count > 1) {
        REQUIRE(global_ghost_count > 0);
    }
}

void test_static_spparks_workflow(int process_count)
{
    for (const auto& process_grid : test_layouts(process_count)) {
        verify_scenario(
            MockScenario::OpenTrench, process_grid, process_count);
        verify_scenario(
            MockScenario::SealedTrench, process_grid, process_count);
    }
}

void test_ghost_cutoff_guard(int rank, int process_count)
{
    const auto fixture_spec = gasaccess::testing::make_mock_fixture_spec(
        MockScenario::OpenTrench);
    const GridDimensions process_grid{
        static_cast<std::uint64_t>(process_count), 1, 1};
    MockSpparksDomain domain(
        MPI_COMM_WORLD,
        fixture_spec.global_lower,
        fixture_spec.global_upper,
        process_grid,
        fixture_spec.spparks_periodic);
    SpparksAdapterConfig adapter_config{};
    adapter_config.requested_spacing = fixture_spec.requested_spacing;
    adapter_config.gas_periodic = fixture_spec.gas_periodic;
    adapter_config.reservoir_faces = fixture_spec.reservoir_faces;
    adapter_config.atom_ghost_distance = 0.49;
    adapter_config.maximum_excluded_radius = 0.50;
    bool rejected = false;
    try {
        DistributedGasGrid grid(
            gasaccess::make_spparks_grid_spec(domain, adapter_config),
            gasaccess::make_spparks_decomposition_spec(
                domain, domain.world, adapter_config));
        (void)grid;
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    REQUIRE(rejected);
    (void)rank;
}

}  // namespace

int main(int argument_count, char** arguments)
{
    MPI_Init(&argument_count, &arguments);
    int rank = 0;
    int process_count = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);

    int local_failure_count = 0;
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"static SPPARKS workflow", [&]() {
             test_static_spparks_workflow(process_count);
         }},
        {"ghost cutoff guard", [&]() {
             test_ghost_cutoff_guard(rank, process_count);
         }}};
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
        std::cout << tests.size() << " SPPARKS mock integration test groups "
                  << "passed on " << process_count << " rank(s)\n";
    }
    MPI_Finalize();
    return global_failure_count == 0 ? 0 : 1;
}
