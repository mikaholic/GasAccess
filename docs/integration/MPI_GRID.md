# MPI grid and halo contract

Phases 11 and 12 provide an optional C++ MPI layer in
`GasAccess::gasaccess_mpi`. It
reuses the host application's Cartesian decomposition; it does not create or
rebalance a second domain decomposition. Distributed initial flood-fill is
implemented; distributed incremental connectivity repair remains Phase 13
work.

## Build

```sh
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON \
    -DGASACCESS_ENABLE_MPI=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The serial target remains `GasAccess::gasaccess`. MPI clients link
`GasAccess::gasaccess_mpi`, which also links the serial core and `MPI::MPI_CXX`.

## Non-cubic aligned geometry

`GridSpec::spacing` contains independent x/y/z components. The helper
`make_aligned_grid_geometry()` takes global bounds, the MPI process grid, and a
requested per-axis resolution. It returns dimensions divisible by the process
count on each axis and the exact spacing obtained from
`box_length / dimension`.

The returned spacing is explicit; a directly supplied `GridSpec` is validated
and never silently modified. Different x/y/z spacings change geometric
resolution but not six-face connectivity. Atom exclusion remains spherical in
physical coordinates and uses an ordinary Euclidean distance test.

## SPPARKS adapter fields

The application constructs `MpiDecompositionSpec` from existing SPPARKS data:

| GasAccess field | SPPARKS source |
|---|---|
| `communicator` | existing `world` communicator |
| `global_lower`, `global_upper` | global domain box bounds |
| `local_lower`, `local_upper` | `[sublo, subhi)` bounds for this rank |
| `process_grid` | `Domain::procgrid` |
| `process_location` | `Domain::myloc` |
| `atom_ghost_distance` | minimum relevant off-lattice bin/ghost distance |
| `maximum_excluded_radius` | maximum `R_atom + R_precursor` for the grid |

The adapter validates the SPPARKS x-fastest rank ordering, communicator size,
global extent, local subdomain bounds, process-grid divisibility, and exact
voxel/subdomain-face alignment. The communicator and MPI runtime must remain
valid for the lifetime of the distributed grid.

## Boundary configuration

GasAccess uses `GridSpec::periodic` and `GridSpec::reservoir_faces`; it does not
copy the KMC periodic flags automatically. A typical all-periodic KMC can be
viewed by GasAccess as:

```cpp
grid_spec.periodic = {true, true, false};
grid_spec.reservoir_faces.z_high = true;
grid_spec.reservoir_faces.z_low = false;
```

This disables gas-state wrapping at both global z faces, makes the top face an
external source, and leaves the bottom face closed/non-source. Internal z rank
interfaces still exchange normally. Synchronized periodic image atoms outside
a GasAccess non-periodic global boundary are ignored during owned-voxel
voxelization.

## Owned and ghost state

Each rank stores its owned rectangular voxel brick plus one padded layer.
Only the six face regions are valid ghosts; edge and corner padding is not
communicated. `exchange_ghost_states()` packs reusable one-byte state buffers,
posts nonblocking point-to-point messages to spatial neighbors, handles
periodic self-neighbors by local copy, and leaves non-periodic global faces
without a neighbor.

One face layer is sufficient because exact alignment guarantees that a site
owned by the KMC rank lies in an owned gas voxel. Its containing voxel and six
face neighbors therefore require only owned or face-ghost state.

The required order is:

```text
SPPARKS synchronizes owned/ghost atoms
        |
        v
GasAccess updates owned occupancy/connectivity
        |
        v
GasAccess exchanges face ghost states
        |
        v
KMC calls query.is_site_accessible(atom_position)
```

`DistributedGasAccessibilityQuery::is_site_accessible()` performs no MPI
operation and allocates no memory. It requires the containing voxel to belong
to the calling rank.

## Atom ghost coverage

Before occupancy construction, GasAccess checks:

```text
atom_ghost_distance >= maximum(R_atom + R_precursor)
```

`voxelize_owned_atoms()` checks every supplied atom against both the declared
maximum and actual ghost distance before mutating state. This is deliberately
defensive even when the KMC bin setup already guarantees sufficient coverage.

The atom input may contain already synchronized owned and ghost atoms;
duplicate coverage is harmless because voxel solidification is idempotent.

## Distributed initial classification

`DistributedExteriorClassifier::classify()` treats owned `Solid` states as
occupancy and recomputes all other owned states. It seeds free voxels on the
configured non-periodic reservoir faces and any configured explicit source
voxels. Fully periodic GasAccess domains therefore need at least one explicit
source voxel when an exterior-connected component is desired.

Each rank exhausts a local six-face breadth-first frontier. When traversal
crosses an owned face, it sends only the tangential face offset to the owning
neighbor; no global voxel IDs, gas graph, or full-grid state are gathered.
Ranks exchange frontier counts and payloads, then use `MPI_Allreduce` only to
decide whether any rank has newly received work. Ranks with no local frontier
continue participating until global termination.

The classifier returns local counts and traversal/communication diagnostics in
`DistributedClassificationSummary`. Counts are local deliberately; an
application may reduce them when it needs global reporting. After termination,
the classifier calls `exchange_ghost_states()` so this sequence is sufficient:

```cpp
gasaccess::DistributedExteriorClassifier classifier;
const auto summary = classifier.classify(distributed_grid);

gasaccess::DistributedGasAccessibilityQuery query(distributed_grid);
const bool accessible = query.is_site_accessible(atom_position);
```

`classify()` is collective over the decomposition communicator. Every rank in
that communicator must call it in the same order. The subsequent position-based
query remains local, allocation-free, and communication-free.
