#include "gasaccess/accessibility_query.hpp"
#include "gasaccess/atom_change_event_buffer.hpp"
#include "gasaccess/atom_change_updater.hpp"
#include "gasaccess/atom_voxelizer.hpp"
#include "gasaccess/distributed_atom_change_updater.hpp"
#include "gasaccess/distributed_exterior_classifier.hpp"
#include "gasaccess/exterior_classifier.hpp"
#include "gasaccess/spparks_adapter.hpp"
#include "mock_spparks.hpp"
#include "mock_spparks_fixture.hpp"

#include <mpi.h>

#include <algorithm>
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
using gasaccess::AccessibilityRepairKind;
using gasaccess::AtomChangeEventBuffer;
using gasaccess::AtomChangeUpdateResult;
using gasaccess::AtomChangeUpdater;
using gasaccess::AtomView;
using gasaccess::AtomVoxelizer;
using gasaccess::DistributedAtomChangeUpdateResult;
using gasaccess::DistributedAtomChangeUpdater;
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
using gasaccess::testing::MockSpparksAtom;
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
                REQUIRE(distributed_grid.owned_blocker_count(coordinate)
                        == serial_grid.blocker_count(coordinate));
                REQUIRE(distributed_grid.gas_state(coordinate)
                        == serial_grid.gas_state(coordinate));
            }
        }
    }
}

void compare_face_ghost_states(
    const DistributedGasGrid& distributed_grid,
    const GasGrid& serial_grid);

struct MockKmcStep {
    std::vector<MockSpparksAtom> added_atoms{};
    std::vector<MockSpparksAtom> removed_atoms{};
    AccessibilityRepairKind expected_repair_kind =
        AccessibilityRepairKind::None;
};

std::vector<MockSpparksAtom> locally_owned_events(
    const MockSpparksDomain& domain,
    const std::vector<MockSpparksAtom>& events)
{
    std::vector<MockSpparksAtom> owned_events;
    for (const auto& event : events) {
        if (domain.owns(event.position)) {
            owned_events.push_back(event);
        }
    }
    return owned_events;
}

void append_mock_app_events(
    const MockSpparksApp& app,
    bool additions,
    AtomChangeEventBuffer& event_buffer)
{
    const auto atom_count = static_cast<std::size_t>(app.nlocal + app.nghost);
    for (std::size_t index = 0; index < atom_count; ++index) {
        const Atom atom{
            {app.xyz[index][0], app.xyz[index][1], app.xyz[index][2]},
            app.radius[index]};
        if (additions) {
            event_buffer.record_adsorption(atom);
        } else {
            event_buffer.record_desorption(atom);
        }
    }
}

AtomChangeEventBuffer synchronize_event_step(
    const MockSpparksDomain& domain,
    const MockKmcStep& step,
    double ghost_distance)
{
    MockSpparksApp added_events(domain);
    added_events.set_owned_atoms(locally_owned_events(
        domain,
        step.added_atoms));
    added_events.synchronize_ghost_atoms(ghost_distance);

    // This separate application-style buffer is intentional: removed atoms
    // carry their old records even after the KMC atom list has discarded them.
    MockSpparksApp removed_events(domain);
    removed_events.set_owned_atoms(locally_owned_events(
        domain,
        step.removed_atoms));
    removed_events.synchronize_ghost_atoms(ghost_distance);

    AtomChangeEventBuffer event_buffer;
    event_buffer.reserve(
        static_cast<std::size_t>(
            added_events.nlocal + added_events.nghost),
        static_cast<std::size_t>(
            removed_events.nlocal + removed_events.nghost));
    append_mock_app_events(added_events, true, event_buffer);
    append_mock_app_events(removed_events, false, event_buffer);
    return event_buffer;
}

AtomChangeEventBuffer global_event_step(const MockKmcStep& step)
{
    AtomChangeEventBuffer event_buffer;
    for (const auto& atom : step.added_atoms) {
        event_buffer.record_adsorption({atom.position, atom.radius});
    }
    for (const auto& atom : step.removed_atoms) {
        event_buffer.record_desorption({atom.position, atom.radius});
    }
    return event_buffer;
}

void apply_to_atom_source(
    std::vector<MockSpparksAtom>& atoms,
    const MockKmcStep& step)
{
    for (const auto& removed_atom : step.removed_atoms) {
        const auto iterator = std::find_if(
            atoms.begin(),
            atoms.end(),
            [&removed_atom](const MockSpparksAtom& atom) {
                return atom.id == removed_atom.id;
            });
        REQUIRE(iterator != atoms.end());
        REQUIRE(iterator->position.x == removed_atom.position.x);
        REQUIRE(iterator->position.y == removed_atom.position.y);
        REQUIRE(iterator->position.z == removed_atom.position.z);
        REQUIRE(iterator->radius == removed_atom.radius);
        atoms.erase(iterator);
    }
    for (const auto& added_atom : step.added_atoms) {
        const auto duplicate = std::find_if(
            atoms.begin(),
            atoms.end(),
            [&added_atom](const MockSpparksAtom& atom) {
                return atom.id == added_atom.id;
            });
        REQUIRE(duplicate == atoms.end());
        atoms.push_back(added_atom);
    }
}

std::vector<Atom> convert_atoms(const std::vector<MockSpparksAtom>& mock_atoms)
{
    std::vector<Atom> atoms;
    atoms.reserve(mock_atoms.size());
    for (const auto& atom : mock_atoms) {
        atoms.push_back({atom.position, atom.radius});
    }
    return atoms;
}

void compare_changed_coordinates(
    const DistributedGasGrid& distributed_grid,
    const DistributedAtomChangeUpdateResult& distributed_result,
    const GasGrid& serial_grid,
    const AtomChangeUpdateResult& serial_result)
{
    std::vector<VoxelCoord> expected_owned_changes;
    for (const auto voxel_id : serial_result.changed_voxel_ids) {
        const auto coordinate = serial_grid.voxel_coord(voxel_id);
        if (distributed_grid.owns(coordinate)) {
            expected_owned_changes.push_back(coordinate);
        }
    }
    REQUIRE(distributed_result.changed_owned_voxel_coords
            == expected_owned_changes);
}

void test_reversible_tkmc_event_sequence(int process_count)
{
    constexpr double atom_radius = 1.1;
    constexpr double precursor_radius = 0.0;
    constexpr double ghost_distance = atom_radius + precursor_radius;
    const GridDimensions process_grid{
        static_cast<std::uint64_t>(process_count), 1, 1};
    const gasaccess::Point3 global_lower{0.0, 0.0, 0.0};
    const gasaccess::Point3 global_upper{32.0, 1.0, 1.0};
    MockSpparksDomain domain(
        MPI_COMM_WORLD,
        global_lower,
        global_upper,
        process_grid,
        {false, false, false});

    gasaccess::SpparksAdapterConfig adapter_config{};
    adapter_config.requested_spacing = {1.0, 1.0, 1.0};
    adapter_config.reservoir_faces.x_low = true;
    adapter_config.atom_ghost_distance = ghost_distance;
    adapter_config.maximum_excluded_radius = ghost_distance;
    const auto grid_spec = gasaccess::make_spparks_grid_spec(
        domain,
        adapter_config);
    const auto decomposition_spec =
        gasaccess::make_spparks_decomposition_spec(
            domain,
            domain.world,
            adapter_config);

    std::vector<MockSpparksAtom> active_atoms{
        {1, {7.5, 0.5, 0.5}, atom_radius},
        {2, {23.5, 0.5, 0.5}, atom_radius}};
    MockSpparksApp app(domain);
    app.set_owned_atoms(locally_owned_events(domain, active_atoms));
    app.synchronize_ghost_atoms(ghost_distance);
    SpparksAtomBuffer initial_atom_buffer;
    initial_atom_buffer.assign(
        app,
        [&app](std::size_t atom_index) {
            return app.radius[atom_index];
        });

    DistributedGasGrid distributed_grid(grid_spec, decomposition_spec);
    distributed_grid.voxelize_owned_atoms(
        initial_atom_buffer.atom_view(),
        precursor_radius);
    DistributedExteriorClassifier{}.classify(distributed_grid);

    GasGrid serial_grid(grid_spec);
    auto initial_serial_atoms = convert_atoms(active_atoms);
    AtomVoxelizer(precursor_radius).voxelize(
        serial_grid,
        {initial_serial_atoms.data(), initial_serial_atoms.size()});
    ExteriorClassifier{}.classify(serial_grid);
    compare_owned_states(distributed_grid, serial_grid);
    compare_face_ghost_states(distributed_grid, serial_grid);

    const MockSpparksAtom atom_2_old = active_atoms[1];
    const MockSpparksAtom atom_2_new{
        2, {24.5, 0.5, 0.5}, atom_radius};
    const MockSpparksAtom atom_3{
        3, {15.5, 0.5, 0.5}, atom_radius};
    const MockSpparksAtom atom_4{
        4, {11.5, 0.5, 0.5}, atom_radius};
    const std::vector<MockKmcStep> steps{
        {{atom_3}, {}, AccessibilityRepairKind::Closing},
        {{}, {active_atoms[0]}, AccessibilityRepairKind::Opening},
        {{atom_2_new}, {atom_2_old}, AccessibilityRepairKind::Mixed},
        {{atom_4}, {atom_3}, AccessibilityRepairKind::Mixed},
        {{}, {}, AccessibilityRepairKind::None},
        {{atom_4}, {atom_4}, AccessibilityRepairKind::None}};

    DistributedAtomChangeUpdater distributed_updater(precursor_radius);
    AtomChangeUpdater serial_updater(precursor_radius);
    for (std::size_t step_index = 0; step_index < steps.size(); ++step_index) {
        const auto& step = steps[step_index];
        auto distributed_events = synchronize_event_step(
            domain,
            step,
            ghost_distance);
        auto serial_events = global_event_step(step);

        DistributedAtomChangeUpdateResult distributed_result{};
        AtomChangeUpdateResult serial_result{};
        if (step_index == 0) {
            distributed_result = distributed_updater.apply_adsorption(
                distributed_grid,
                distributed_events.added_atoms());
            serial_result = serial_updater.apply_adsorption(
                serial_grid,
                serial_events.added_atoms());
        } else if (step_index == 1) {
            distributed_result = distributed_updater.apply_desorption(
                distributed_grid,
                distributed_events.removed_atoms());
            serial_result = serial_updater.apply_desorption(
                serial_grid,
                serial_events.removed_atoms());
        } else {
            distributed_result = distributed_updater.apply_atom_changes(
                distributed_grid,
                distributed_events.atom_changes());
            serial_result = serial_updater.apply_atom_changes(
                serial_grid,
                serial_events.atom_changes());
        }

        REQUIRE(distributed_result.repair_kind
                == step.expected_repair_kind);
        REQUIRE(serial_result.repair_kind == step.expected_repair_kind);
        REQUIRE(distributed_result.global_newly_solid_count
                == serial_result.newly_solid_count);
        REQUIRE(distributed_result.global_newly_gas_count
                == serial_result.newly_gas_count);
        compare_changed_coordinates(
            distributed_grid,
            distributed_result,
            serial_grid,
            serial_result);

        apply_to_atom_source(active_atoms, step);
        GasGrid reconstructed_grid(grid_spec);
        const auto reconstructed_atoms = convert_atoms(active_atoms);
        AtomVoxelizer(precursor_radius).voxelize(
            reconstructed_grid,
            {reconstructed_atoms.data(), reconstructed_atoms.size()});
        ExteriorClassifier{}.classify(reconstructed_grid);
        compare_owned_states(distributed_grid, reconstructed_grid);
        compare_face_ghost_states(distributed_grid, reconstructed_grid);
        compare_owned_states(distributed_grid, serial_grid);

        const DistributedGasAccessibilityQuery distributed_query(
            distributed_grid);
        const GasAccessibilityQuery serial_query(reconstructed_grid);
        const auto& owned_range = distributed_grid.owned_range();
        for (auto x = owned_range.begin.x; x < owned_range.end.x; ++x) {
            const auto center = distributed_grid.voxel_center({x, 0, 0});
            REQUIRE(distributed_query.is_site_accessible(center)
                    == serial_query.is_site_accessible(center));
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
         }},
        {"reversible tKMC event sequence", [&]() {
             test_reversible_tkmc_event_sequence(process_count);
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
