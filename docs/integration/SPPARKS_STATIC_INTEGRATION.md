# SPPARKS-style static integration

Phase 14 supplies a thin adapter and a self-contained MPI application that
exercise the GasAccess workflow without requiring the unavailable production
KMC application. The mock deliberately mirrors the public fields used from
SPPARKS `Domain` and `App`; it does not reimplement the SPPARKS simulator.

## What the adapter consumes

`include/gasaccess/spparks_adapter.hpp` maps these existing public fields:

- `Domain::boxxlo` through `boxzhi`;
- `Domain::subxlo` through `subzhi`;
- `Domain::procgrid` and `Domain::myloc`;
- `Domain::box_exist` and `Domain::dimension`;
- `App::nlocal`, `App::nghost`, and `App::xyz`.

Atom radius is intentionally supplied by an application callback because the
SPPARKS base `App` class has no universal atom-radius field. The caller also
supplies nominal voxel spacing, GasAccess periodic axes, reservoir faces,
atom ghost distance, and the maximum excluded radius. GasAccess does not copy
the SPPARKS boundary flags automatically.

The optional CMake setting below compiles the adapter directly against the real
SPPARKS headers as a compatibility check:

```sh
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON \
    -DGASACCESS_ENABLE_MPI=ON \
    -DGASACCESS_SPPARKS_SOURCE_DIR=/path/to/spparks/src
cmake --build build --parallel
```

No SPPARKS object file or runtime library is linked by this compile-only test.

## Production call sequence

After the host application has performed its normal owned/ghost atom
synchronization, the static initialization is:

```cpp
gasaccess::SpparksAdapterConfig config{};
config.requested_spacing = {dx, dy, dz};
config.gas_periodic = {true, true, false};
config.reservoir_faces.z_high = true;
config.atom_ghost_distance = ghost_distance;
config.maximum_excluded_radius = maximum_atom_radius + precursor_radius;

const auto grid_spec = gasaccess::make_spparks_grid_spec(*domain, config);
const auto decomposition = gasaccess::make_spparks_decomposition_spec(
    *domain, world, config);
gasaccess::DistributedGasGrid gas_grid(grid_spec, decomposition);

gasaccess::SpparksAtomBuffer atoms;
atoms.assign(*app, [&](std::size_t atom_index) {
    return application_atom_radius(atom_index);
});
gas_grid.voxelize_owned_atoms(atoms.atom_view(), precursor_radius);

gasaccess::DistributedExteriorClassifier classifier;
classifier.classify(gas_grid);
gasaccess::DistributedGasAccessibilityQuery query(gas_grid);
```

The KMC site loop needs only:

```cpp
const bool accessible = query.is_site_accessible(atom_position);
```

That call uses cached owned/face-ghost gas states and performs no MPI
communication. The KMC remains responsible for its ordinary rate scan.

## Reversible event sequence

At a caller-selected synchronization point, the host synchronizes added event
records at their new positions and removed event records at their old
positions. `AtomChangeEventBuffer` owns those records while GasAccess performs
the collective update:

```cpp
gasaccess::AtomChangeEventBuffer events;
for (const auto& atom : synchronized_added_events) {
    events.record_adsorption(atom);
}
for (const auto& atom : synchronized_removed_events) {
    events.record_desorption(atom);  // old position and radius
}

gasaccess::DistributedAtomChangeUpdater updater(precursor_radius);
const auto result = updater.apply_atom_changes(
    gas_grid,
    events.atom_changes());
events.clear();
```

The buffer deliberately does not choose synchronization timing or perform atom
communication. tKMC retains that flexibility and may accumulate several
sector-safe events before synchronizing. Every rank must enter the collective
in the same order, including ranks with empty buffers. Accessibility queries
remain local and communication-free after the updater returns.

The existing `DistributedAdsorptionUpdater` and
`DistributedDesorptionUpdater` interfaces remain available. The unified
updater also provides `apply_adsorption()` and `apply_desorption()` convenience
methods for applications that prefer one persistent update object.

## Mock application

`MockSpparksDomain` uses SPPARKS-compatible x-fastest Cartesian rank numbering,
`[sublo, subhi)` ownership, `procneigh`, global/local bounds, and independent
physical periodic flags. `MockSpparksApp` exposes owned atoms first and ghost
atoms second through `nlocal`, `nghost`, `id`, and `xyz` arrays. Its test-only
ghost synchronization sends atoms only to relevant spatial neighbors; it does
not gather the global structure.

The fixtures intentionally keep the mock SPPARKS domain periodic in x/y/z but
call GasAccess with periodic x/y, non-periodic z, and a top-z gas reservoir.
They cover:

- an open trench whose internal wall probes are accessible;
- the same trench with a cap, whose internal wall probes are inaccessible;
- a 128 x 128 x 64 solid slab containing 1,048,576 atoms beneath an equally
  sized gas region.
- a reversible six-step tKMC sequence containing adsorption, desorption, an
  atom move, a mixed batch, an empty batch, and a no-net-change batch.

The reversible sequence synchronizes added and old removed-event buffers
independently, including ghost copies for atom footprints that cross ownership
boundaries. After every collective update it compares owned blocker counts,
states, face ghosts, changed coordinates, and cached queries with serial and
freshly reconstructed references at 1, 2, 4, and 8 ranks.

Run them with:

```sh
mpiexec -n 1 ./build/gasaccess_spparks_mock_driver --scenario open-trench
mpiexec -n 4 ./build/gasaccess_spparks_mock_driver --scenario sealed-trench
mpiexec -n 4 ./build/gasaccess_spparks_mock_driver --scenario million-slab
```

The driver prints machine-readable counts, communication volume, peak RSS,
per-stage maximum rank time, total time, and query throughput. It exits with an
error if the expected accessibility result is not reproduced.

## Real-application handoff

The actual tKMC application is not present in this repository, so GasAccess
cannot name its atom-radius accessor, select its synchronization point, or
invoke its existing atom/event exchange. Those application-specific calls are
the handoff. Grid construction, classification, the owning reversible event
contract, collective repair, and site-loop queries are independent of the mock
types. The SPPARKS field adapter also compiles directly against an available
real SPPARKS source tree.
