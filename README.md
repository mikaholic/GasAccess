# GasAccess

GasAccess is a C++17 library for geometric gas-accessibility analysis on a
Cartesian voxel grid. Development follows the recorded phased plan in
[`docs/plans/GAS_ACCESS_DEVELOPMENT_PLAN.md`](docs/plans/GAS_ACCESS_DEVELOPMENT_PLAN.md).

The current standalone implementation is complete through Phase 14: the serial
geometry/connectivity/update/query baseline, scale optimization, non-cubic
voxel support, distributed initial classification, and distributed incremental
repair for monotonic deposition, plus a SPPARKS-compatible adapter and static
MPI acceptance application.

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
- an exception-safe C99 API with opaque grid and update-result handles;
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
`build/libgasaccess_mpi.a`, the standalone reference driver, the nine serial
test executables, `build/gasaccess_mpi_grid_tests`,
`build/gasaccess_distributed_classifier_tests`, and
`build/gasaccess_distributed_updater_tests`. Phase 14 also produces
`build/gasaccess_spparks_mock_driver` and
`build/gasaccess_spparks_mock_integration_tests`.

- `build/gasaccess_grid_tests`
- `build/gasaccess_atom_voxelizer_tests`
- `build/gasaccess_exterior_classifier_tests`
- `build/gasaccess_accessibility_query_tests`
- `build/gasaccess_kmc_query_integration_tests`
- `build/gasaccess_c_api_tests`
- `build/gasaccess_deposition_updater_tests`
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
./build/gasaccess_deposition_updater_tests
./build/gasaccess_local_topology_filter_tests
./build/gasaccess_affected_region_repair_tests
```

The registered MPI tests exercise one, two, and four ranks:

```sh
ctest --test-dir build --output-on-failure -L mpi
```

The MPI data contract, SPPARKS field mapping, boundary override, aligned-grid
helper, and update/query ordering are documented in
[`docs/integration/MPI_GRID.md`](docs/integration/MPI_GRID.md).
The concrete adapter, mock application, and remaining production handoff are
documented in
[`docs/integration/SPPARKS_STATIC_INTEGRATION.md`](docs/integration/SPPARKS_STATIC_INTEGRATION.md).

## SPPARKS-style static acceptance

Run the open trench, sealed trench, and million-atom slab through the same
owned/ghost atom, distributed voxelization, flood-fill, and site-query sequence
that the real application will use:

```sh
mpiexec -n 1 ./build/gasaccess_spparks_mock_driver --scenario open-trench
mpiexec -n 4 ./build/gasaccess_spparks_mock_driver --scenario sealed-trench
mpiexec -n 4 ./build/gasaccess_spparks_mock_driver --scenario million-slab
```

The driver does not invoke deposition or a KMC event. It adds a static atom
structure, synchronizes mock SPPARKS-style ghosts, constructs and classifies
the gas grid, then calls
`query.is_site_accessible(atom_position)` for each owned atom. Phase 14 timing
results are recorded in
[`docs/benchmarks/PHASE14_SPPARKS_MOCK.md`](docs/benchmarks/PHASE14_SPPARKS_MOCK.md).

## Reference driver

Display all driver options:

```sh
./build/gasaccess_reference_driver --help
```

Run the default open-, sealed-, and deposition-driven pinch-off fixtures:

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
voxelization, initial classification, cached queries, and deposition updates.
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
