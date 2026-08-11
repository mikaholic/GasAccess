# GasAccess

GasAccess is a C++17 library for geometric gas-accessibility analysis on a
Cartesian voxel grid. Development follows the recorded phased plan in
[`docs/plans/GAS_ACCESS_DEVELOPMENT_PLAN.md`](docs/plans/GAS_ACCESS_DEVELOPMENT_PLAN.md).

The current implementation provides Phase 1 grid geometry/topology and Phase 2
static atom voxelization:

- dense one-byte gas-state storage;
- checked 64-bit voxel identifiers;
- voxel coordinate/index and world-coordinate conversion;
- allocation-free six-face neighbor lookup;
- independently periodic x, y, and z axes;
- non-periodic reservoir faces and explicit reservoir source voxels;
- non-owning in-memory atom input with double-precision positions and radii;
- additive spherical steric exclusion using `R_atom + R_precursor`;
- bounded per-atom candidate traversal with periodic-image distance handling.

Gas-connectivity flood-fill is Phase 3 and is not implemented yet.

## Requirements

- CMake 3.16 or newer
- A C++17 compiler

The core and its current tests have no third-party dependencies.

## Compile

Configure and compile an optimized build:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
```

This produces the static library `build/libgasaccess.a` and two test executables:

- `build/gasaccess_grid_tests`
- `build/gasaccess_atom_voxelizer_tests`

## Test

```sh
ctest --test-dir build --output-on-failure
```

To display every individual test-group result directly:

```sh
./build/gasaccess_grid_tests
./build/gasaccess_atom_voxelizer_tests
```
