# GasAccess

GasAccess is a C++17 library for geometric gas-accessibility analysis on a
Cartesian voxel grid. Development follows the recorded phased plan in
[`docs/plans/GAS_ACCESS_DEVELOPMENT_PLAN.md`](docs/plans/GAS_ACCESS_DEVELOPMENT_PLAN.md).

The current implementation includes the original serial and distributed
geometry/connectivity/query baseline, distributed incremental repair for
monotonic adsorption, and Phase R3 serial and distributed incremental
desorption repair. Phase R4 adds atomic mixed adsorption/desorption repair in
both the serial and distributed C++ layers. Phase R5 exposes the reversible
path through the C API and a KMC-oriented owning event buffer.

- dense one-byte gas-state storage;
- checked 64-bit voxel identifiers;
- independent finite positive x/y/z voxel spacings;
- voxel coordinate/index and world-coordinate conversion;
- allocation-free six-face neighbor lookup;
- independently periodic x, y, and z axes;
- non-periodic reservoir faces and explicit reservoir source voxels;
- non-owning in-memory atom input with double-precision positions and radii;
- additive spherical steric exclusion using `R_atom + R_precursor`;
- bounded per-atom candidate traversal with periodic-image distance handling;
- full serial six-neighbor flood-fill from configured reservoir sources;
- cached `Solid`, `OutsideAccessible`, and `ClosedVoid` voxel states;
- incrementally maintained solid/outside/closed state counts;
- checked per-voxel atom-blocker counts, stored only for owned voxels under MPI;
- allocation-free cached queries using either a seven-voxel default stencil or
  a caller-supplied stencil of at most 27 voxel IDs;
- a Phase 9 KMC contract requiring only
  `query.is_site_accessible(atom_position)` inside the site loop;
- fixed-capacity `3x3x3` local connectivity checks for single-voxel removal;
- conservative escalation for sources, suspected pinch-offs, and multi-voxel
  changes;
- source-aware component repair from surviving neighbors of newly solid voxels;
- exact newly-solid voxel capture without a production full-grid discovery scan;
- reusable epoch/frontier traversal storage with periodic-boundary support;
- selectable incremental or forced full-reference update modes;
- sorted changed-voxel reporting and post-update classification counts;
- serial and distributed incremental opening repair after atom desorption,
  with overlap-safe occupancy and forced full-reclassification reference modes;
- atomic mixed atom-change batches with incremental closing-then-opening
  repair and one final MPI ghost-state synchronization;
- an owning KMC event buffer that retains added and old removed-atom records;
- an exception-safe C99 API with separate adsorption-specific and reversible
  update-result handles;
- deterministic open-trench, sealed-trench, and bulk workload generation;
- machine-readable correctness checksums, timings, and memory reporting;
- MPI-aligned global dimensions and exact integer owned ranges;
- one layer of face-connected gas ghost states with reusable communication
  buffers;
- local-only distributed site queries after gas-state halo exchange;
- defensive atom-ghost-distance validation;
- distributed six-face flood-fill using local breadth-first traversal and
  compact face-frontier messages;
- global termination detection without gathering or replicating the gas grid;
- final gas-state halo synchronization before cached site queries;
- distributed change-capturing atom voxelization with owned state counters;
- conservative constant-size topology filtering with rank-boundary escalation;
- affected-component MPI repair using compact face-frontier messages;
- forced distributed full reclassification for debugging and differential
  validation;
- templated mapping from the public SPPARKS `Domain` and `App` field contracts;
- independent SPPARKS physical and GasAccess reservoir-boundary settings;
- a neighbor-directed mock SPPARKS atom halo, with no global atom gather;
- static open/sealed trench acceptance and a million-atom MPI benchmark;
- optional direct compile checking against an available SPPARKS source tree.

The Phase 3 classifier is the correctness-reference implementation that later
incremental update algorithms are tested against. The focused KMC query
contract is documented in
[`docs/integration/KMC_SITE_QUERY.md`](docs/integration/KMC_SITE_QUERY.md).

## Requirements

- CMake 3.16 or newer
- A C++17 compiler
- A C99 compiler when building the C API test client
- An MPI C++ implementation when `GASACCESS_ENABLE_MPI=ON`

The serial core has no third-party library dependencies. MPI is the only
dependency of the optional `gasaccess_mpi` target.

## Compile

Configure and compile an optimized build:

```sh
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON \
    -DGASACCESS_ENABLE_MPI=ON \
    -DGASACCESS_SPPARKS_SOURCE_DIR=/home/mikaholic/project/spparks/src
cmake --build build --parallel
```

Leave `GASACCESS_ENABLE_MPI` off, its default, for a serial-only build.
`GASACCESS_SPPARKS_SOURCE_DIR` is optional; when provided, an adapter target is
compiled against the real SPPARKS base-class headers without linking SPPARKS.

With MPI enabled this produces `build/libgasaccess.a`,
`build/libgasaccess_mpi.a`, the standalone reference driver, the eleven serial
test executables, `build/gasaccess_mpi_grid_tests`,
`build/gasaccess_distributed_classifier_tests`,
`build/gasaccess_distributed_updater_tests`,
`build/gasaccess_distributed_desorption_tests`, and
`build/gasaccess_distributed_atom_change_tests`. Phase 14 also produces
`build/gasaccess_spparks_mock_driver` and
`build/gasaccess_spparks_mock_integration_tests`. The lifecycle efficiency
matrix is provided by `build/gasaccess_mpi_efficiency_driver`.

- `build/gasaccess_grid_tests`
- `build/gasaccess_atom_voxelizer_tests`
- `build/gasaccess_exterior_classifier_tests`
- `build/gasaccess_accessibility_query_tests`
- `build/gasaccess_kmc_query_integration_tests`
- `build/gasaccess_c_api_tests`
- `build/gasaccess_adsorption_updater_tests`
- `build/gasaccess_desorption_updater_tests`
- `build/gasaccess_atom_change_updater_tests`
- `build/gasaccess_local_topology_filter_tests`
- `build/gasaccess_affected_region_repair_tests`

## Test

```sh
ctest --test-dir build --output-on-failure
```

To display every individual test-group result directly:

```sh
./build/gasaccess_grid_tests
./build/gasaccess_atom_voxelizer_tests
./build/gasaccess_exterior_classifier_tests
./build/gasaccess_accessibility_query_tests
./build/gasaccess_kmc_query_integration_tests
./build/gasaccess_c_api_tests
./build/gasaccess_adsorption_updater_tests
./build/gasaccess_desorption_updater_tests
./build/gasaccess_atom_change_updater_tests
./build/gasaccess_local_topology_filter_tests
./build/gasaccess_affected_region_repair_tests
```

The registered MPI tests exercise one, two, four, and eight ranks by default:

```sh
ctest --test-dir build --output-on-failure -L mpi
```

Override the rank matrix at configure time when a machine has a different MPI
slot budget:

```sh
cmake -S . -B build \
    -DGASACCESS_ENABLE_MPI=ON \
    -DGASACCESS_MPI_TEST_RANK_COUNTS="1;2;4;8"
```

The MPI data contract, SPPARKS field mapping, boundary override, aligned-grid
helper, and update/query ordering are documented in
[`docs/integration/MPI_GRID.md`](docs/integration/MPI_GRID.md).
The concrete adapter, mock application, and remaining production handoff are
documented in
[`docs/integration/SPPARKS_STATIC_INTEGRATION.md`](docs/integration/SPPARKS_STATIC_INTEGRATION.md).

## Incremental desorption

The serial C++ layer can remove one or more atoms using their old positions
and radii:

```cpp
gasaccess::DesorptionUpdater updater(precursor_radius);
const auto result = updater.apply_desorption(
    gas_grid,
    {removed_atoms, removed_atom_count});
```

Blocker counts prevent a voxel from becoming gas while another atom still
overlaps it. When the final blocker is removed, GasAccess seeds only newly gas
voxels connected to a reservoir or existing accessible gas and flood-fills the
reachable `ClosedVoid` components. Empty batches and overlap-only removals skip
the flood fill. Details and correctness results are recorded in
[`docs/benchmarks/PHASE_R2_SERIAL_DESORPTION.md`](docs/benchmarks/PHASE_R2_SERIAL_DESORPTION.md).

The MPI layer provides the same operation collectively:

```cpp
gasaccess::DistributedDesorptionUpdater updater(precursor_radius);
const auto result = updater.apply_desorption(
    distributed_grid,
    {removed_atoms, removed_atom_count});
```

Every rank calls it in the same order, including ranks with an empty local
atom view. Each rank updates only owned blocker counts; opening propagation
crosses subdomains through compact face-frontier offsets, terminates with a
global activity check, and synchronizes final face ghosts once. Phase R3's MPI
contract and 1/2/4/8-rank correctness results are recorded in
[`docs/benchmarks/PHASE_R3_DISTRIBUTED_DESORPTION.md`](docs/benchmarks/PHASE_R3_DISTRIBUTED_DESORPTION.md).

## Atomic mixed atom changes

Phase R4 accepts adsorbed and desorbed atoms in one transaction. Removed
records use their old positions and radii; an atom move is represented by one
old record and one new record in the same batch:

```cpp
gasaccess::AtomChangeUpdater updater(precursor_radius);
const auto result = updater.apply_atom_changes(
    gas_grid,
    {{added_atoms, added_atom_count},
     {removed_atoms, removed_atom_count}});
```

The MPI equivalent is collective:

```cpp
gasaccess::DistributedAtomChangeUpdater updater(precursor_radius);
const auto result = updater.apply_atom_changes(
    distributed_grid,
    {{added_atoms, added_atom_count},
     {removed_atoms, removed_atom_count}});
```

GasAccess aggregates both directions before committing blocker counts. It then
runs closing repair for newly solid voxels followed by opening repair for newly
gas voxels, so the final connectivity—not either temporary event ordering—is
reported. Under MPI, compact boundary queries bridge the two passes without
an intermediate full ghost exchange; final face ghosts are synchronized once.
`AtomChangeRepairMode::FullReclassification` retains a debugging and
differential-testing reference path. Details are recorded in
[`docs/benchmarks/PHASE_R4_MIXED_ATOM_CHANGES.md`](docs/benchmarks/PHASE_R4_MIXED_ATOM_CHANGES.md).

Callers may retain one updater for all event kinds. The owning event buffer
keeps old desorption records valid through the update:

```cpp
gasaccess::AtomChangeEventBuffer events;
events.record_adsorption(added_atom);
events.record_desorption(removed_atom_at_old_position);
events.record_move(old_atom, new_atom);

const auto result = updater.apply_atom_changes(
    distributed_grid,
    events.atom_changes());
events.clear();
```

Under MPI, the host application must synchronize both event directions far
enough to cover every affected owned voxel before recording the rank-local
views. Every rank then calls the updater collectively, including ranks whose
event buffer is empty.

## Reversible C API

The C99 interface exposes `ga_apply_adsorption()` with an adsorption-specific
result contract. Reversible clients use `ga_atom_change_batch` and the separate
`ga_atom_change_result` lifecycle:

```c
ga_atom_change_batch changes = {
    {added_atoms, added_atom_count},
    {removed_atoms_at_old_positions, removed_atom_count}
};
ga_atom_change_result* result = NULL;

ga_status status = ga_apply_atom_changes(
    grid, &changes, precursor_radius, &result);
ga_atom_change_update_summary summary;
ga_atom_change_result_get_summary(result, &summary);
ga_atom_change_result_destroy(result);
```

`ga_apply_desorption()` is the pure-removal convenience entry point. Both
operations also have `_with_mode` variants for forced full reclassification.
The Phase R5 interface and mock KMC acceptance results are documented in
[`docs/benchmarks/PHASE_R5_PUBLIC_INTEGRATION.md`](docs/benchmarks/PHASE_R5_PUBLIC_INTEGRATION.md).

## SPPARKS-style static acceptance

Run the open trench, sealed trench, and million-atom slab through the same
owned/ghost atom, distributed voxelization, flood-fill, and site-query sequence
that the real application will use:

```sh
mpiexec -n 1 ./build/gasaccess_spparks_mock_driver --scenario open-trench
mpiexec -n 4 ./build/gasaccess_spparks_mock_driver --scenario sealed-trench
mpiexec -n 4 ./build/gasaccess_spparks_mock_driver --scenario million-slab
```

The driver does not invoke adsorption or a KMC event. It adds a static atom
structure, synchronizes mock SPPARKS-style ghosts, constructs and classifies
the gas grid, then calls
`query.is_site_accessible(atom_position)` for each owned atom. Phase 14 timing
results are recorded in
[`docs/benchmarks/PHASE14_SPPARKS_MOCK.md`](docs/benchmarks/PHASE14_SPPARKS_MOCK.md).

## MPI lifecycle efficiency driver

The repeated efficiency driver covers KMC initialization, cached coordinate
queries, and collective adsorption, desorption, or mixed atom-change repair:

```sh
mpiexec -n 8 ./build/gasaccess_mpi_efficiency_driver \
    --operation initialization --scenario million-slab

mpiexec -n 8 ./build/gasaccess_mpi_efficiency_driver \
    --operation query --scenario million-slab --query-count 1000000

mpiexec -n 8 ./build/gasaccess_mpi_efficiency_driver \
    --operation repair --change-kind desorption \
    --case worst --nx 128 --ny 128 --nz 128

mpiexec -n 8 ./build/gasaccess_mpi_efficiency_driver \
    --operation repair --change-kind mixed \
    --case medium --nx 128 --ny 128 --nz 128
```

Fast operations repeat until the requested cumulative measured duration is
reached. The primary result is average slowest-rank operation time. Each repair
case is checked against a paired forced-full reclassification, and closing and
opening traversal metrics are reported separately. Adsorption is the default
`--change-kind`. The original lifecycle driver is
documented in
[`docs/benchmarks/PHASE15_LIFECYCLE_EFFICIENCY.md`](docs/benchmarks/PHASE15_LIFECYCLE_EFFICIENCY.md);
the current best/medium/worst MPI scaling results are in
[`docs/benchmarks/ATOM_CHANGE_MPI_SCALING.md`](docs/benchmarks/ATOM_CHANGE_MPI_SCALING.md).
The original Phase R6 matrix is retained in
[`docs/benchmarks/PHASE_R6_REVERSIBLE_EFFICIENCY.md`](docs/benchmarks/PHASE_R6_REVERSIBLE_EFFICIENCY.md).

## Reference driver

Display all driver options:

```sh
./build/gasaccess_reference_driver --help
```

Run the default open-, sealed-, and adsorption-driven pinch-off fixtures:

```sh
./build/gasaccess_reference_driver --scenario open-trench
./build/gasaccess_reference_driver --scenario sealed-trench
./build/gasaccess_reference_driver \
    --scenario pinch-off --nx 64 --ny 32 --nz 64 --update-count 640
```

Run a deterministic million-atom synthetic structure:

```sh
./build/gasaccess_reference_driver \
    --scenario bulk --nx 128 --ny 64 --nz 128 \
    --atom-count 1000000 --query-count 1000000 --update-count 3 \
    --seed 12345
```

The output uses one `key=value` field per line so it can be archived or parsed
by benchmark automation. Timings cover grid construction, atom generation,
voxelization, initial classification, cached queries, and adsorption updates.
It also reports median, p95, and maximum update latency separated into locally
safe, affected-region repair, and full-reference paths. Memory output separates
the exact one-byte-per-voxel persistent state from a traversal-frontier
logical-payload bound and whole-process peak RSS.

The recorded Phase 6 environment, commands, results, and interpretation are in
[`docs/benchmarks/PHASE6_BASELINE.md`](docs/benchmarks/PHASE6_BASELINE.md).
The first local-filter comparison is in
[`docs/benchmarks/PHASE7_FILTER.md`](docs/benchmarks/PHASE7_FILTER.md).
Phase 8 repair locality and correctness results are in
[`docs/benchmarks/PHASE8_REPAIR.md`](docs/benchmarks/PHASE8_REPAIR.md).
Phase 10 serial optimization and scale results are in
[`docs/benchmarks/PHASE10_SERIAL_OPTIMIZATION.md`](docs/benchmarks/PHASE10_SERIAL_OPTIMIZATION.md).
