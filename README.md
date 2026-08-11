# GasAccess

GasAccess is a C++17 library for geometric gas-accessibility analysis on a
Cartesian voxel grid. Development follows the recorded phased plan in
[`docs/plans/GAS_ACCESS_DEVELOPMENT_PLAN.md`](docs/plans/GAS_ACCESS_DEVELOPMENT_PLAN.md).

The current implementation provides Phase 1 grid geometry and topology:

- dense one-byte gas-state storage;
- checked 64-bit voxel identifiers;
- voxel coordinate/index and world-coordinate conversion;
- allocation-free six-face neighbor lookup;
- independently periodic x, y, and z axes;
- non-periodic reservoir faces and explicit reservoir source voxels.

Atom voxelization and gas-connectivity flood-fill are later phases and are not
implemented yet.

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

This produces the static library `build/libgasaccess.a` and the test executable
`build/gasaccess_grid_tests`.

## Test

```sh
ctest --test-dir build --output-on-failure
```

To display every individual test-group result directly:

```sh
./build/gasaccess_grid_tests
```
