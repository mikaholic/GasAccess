# GasAccess

GasAccess is a C++17 library for geometric gas-accessibility analysis on a
Cartesian voxel grid. Development follows the recorded phased plan in
[`docs/plans/GAS_ACCESS_DEVELOPMENT_PLAN.md`](docs/plans/GAS_ACCESS_DEVELOPMENT_PLAN.md).

The current implementation provides Phase 1 grid geometry/topology, Phase 2
static atom voxelization, Phase 3 exterior classification, Phase 4 cached site
queries with C interoperability, the Phase 5 deposition-update baseline, and
the Phase 6 standalone reference driver:

- dense one-byte gas-state storage;
- checked 64-bit voxel identifiers;
- voxel coordinate/index and world-coordinate conversion;
- allocation-free six-face neighbor lookup;
- independently periodic x, y, and z axes;
- non-periodic reservoir faces and explicit reservoir source voxels;
- non-owning in-memory atom input with double-precision positions and radii;
- additive spherical steric exclusion using `R_atom + R_precursor`;
- bounded per-atom candidate traversal with periodic-image distance handling;
- full serial six-neighbor flood-fill from configured reservoir sources;
- cached `Solid`, `OutsideAccessible`, and `ClosedVoid` voxel states;
- solid/outside/closed classification summary counts;
- allocation-free cached queries using either a seven-voxel default stencil or
  a caller-supplied stencil of at most 27 voxel IDs;
- full connectivity recomputation after a caller-supplied deposition batch;
- sorted changed-voxel reporting and post-update classification counts;
- an exception-safe C99 API with opaque grid and update-result handles;
- deterministic open-trench, sealed-trench, and bulk workload generation;
- machine-readable correctness checksums, timings, and memory reporting.

The Phase 3 classifier is the correctness-reference implementation that later
incremental update algorithms will be tested against.

## Requirements

- CMake 3.16 or newer
- A C++17 compiler
- A C99 compiler when building the C API test client

The core and its current tests have no third-party library dependencies.

## Compile

Configure and compile an optimized build:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
```

This produces the static library `build/libgasaccess.a`, the standalone
`build/gasaccess_reference_driver`, and six dedicated test executables:

- `build/gasaccess_grid_tests`
- `build/gasaccess_atom_voxelizer_tests`
- `build/gasaccess_exterior_classifier_tests`
- `build/gasaccess_accessibility_query_tests`
- `build/gasaccess_c_api_tests`
- `build/gasaccess_deposition_updater_tests`

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
./build/gasaccess_c_api_tests
./build/gasaccess_deposition_updater_tests
```

## Reference driver

Display all driver options:

```sh
./build/gasaccess_reference_driver --help
```

Run the default open- and sealed-trench fixtures:

```sh
./build/gasaccess_reference_driver --scenario open-trench
./build/gasaccess_reference_driver --scenario sealed-trench
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
voxelization, initial classification, cached queries, and Phase 5 full-update
recomputation. Memory output separates the exact one-byte-per-voxel persistent
state from a traversal-frontier logical-payload bound and whole-process peak
RSS.

The recorded Phase 6 environment, commands, results, and interpretation are in
[`docs/benchmarks/PHASE6_BASELINE.md`](docs/benchmarks/PHASE6_BASELINE.md).
