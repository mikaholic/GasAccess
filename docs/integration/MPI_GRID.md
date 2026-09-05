# MPI grid and halo contract

Phases 11 through 13 provide an optional C++ MPI layer in
`GasAccess::gasaccess_mpi`. It
reuses the host application's Cartesian decomposition; it does not create or
rebalance a second domain decomposition. Both distributed initial flood-fill
and incremental deposition/desorption connectivity repair are implemented.
Phase 14 adds the SPPARKS field adapter and static acceptance application
documented in
[`SPPARKS_STATIC_INTEGRATION.md`](SPPARKS_STATIC_INTEGRATION.md).

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

`make_spparks_decomposition_spec()` constructs `MpiDecompositionSpec` from
existing SPPARKS data:

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

The atom input may contain already synchronized owned and ghost atoms, but each
physical atom must appear exactly once in a rank's view. Blocker counts retain
overlap multiplicity so a duplicated atom record would require a matching
duplicate removal before its voxels become gas.

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

## Distributed incremental deposition update

`DistributedDepositionUpdater::apply_deposition()` is the production update
path for monotonic gas-to-solid changes. The caller supplies deposited atoms
after the normal KMC owned/ghost synchronization:

```cpp
gasaccess::DistributedDepositionUpdater updater(precursor_radius);

const auto result = updater.apply_deposition(
    distributed_grid,
    {deposited_atoms, deposited_atom_count});
```

Every rank in the grid communicator must call the updater in the same order,
including ranks whose local atom view produces no newly solid owned voxel. The
updater performs these steps:

1. Solidify owned voxels and retain their previous classified states.
2. Reduce the global changed and previously accessible voxel counts.
3. Skip connectivity work when only existing closed-void gas was removed.
4. For a single accessible removal, attempt a fixed `3x3x3` local proof.
5. Escalate multi-voxel changes and unprovable rank-boundary cases to
   distributed affected-component repair.
6. Synchronize final face ghosts before returning.

The local proof uses only owned data for its `3x3x3` traversal. Zero- and
one-neighbor cases are still provably safe using face ghosts. When a larger
proof would need an unavailable edge or corner ghost, the updater conservatively
repairs instead of declaring the event safe.

Repair synchronizes only the small removed-voxel metadata list, then explores
stale `OutsideAccessible` components from surviving neighbors. Traversal uses
compact tangential face offsets and termination/source reductions. Components
that reach a configured source stay accessible; exhausted components are
relabelled `ClosedVoid` on their owning ranks. No full gas grid is gathered.

`DistributedDepositionUpdateResult::changed_owned_voxel_coords` contains local
owned coordinates whose state changed. It is diagnostic; the KMC integration
does not need atom registration or selective propensity invalidation. After the
update, KMC continues to call only:

```cpp
const bool accessible = query.is_site_accessible(atom_position);
```

For debugging, construct the updater with
`ConnectivityRepairMode::FullReclassification`. Every geometry-changing call
then uses the Phase 12 classifier and reports all changed owned coordinates.

This deposition entry point remains restricted to gas-to-solid changes. Use
the desorption entry point below for atom removals, or the unified atom-change
entry point for additions and removals in one atomic batch.

## Distributed incremental desorption update

`DistributedDesorptionUpdater::apply_desorption()` handles synchronized
solid-to-gas changes using removed atoms at their old positions and radii:

```cpp
gasaccess::DistributedDesorptionUpdater updater(precursor_radius);

const auto result = updater.apply_desorption(
    distributed_grid,
    {removed_atoms, removed_atom_count});
```

Every rank calls this operation in the same order, including ranks with an
empty local atom view. The caller must keep each removed atom record available
until the collective returns and supply it exactly once in each synchronized
rank view where it is needed for owned-voxel coverage.

The updater performs these steps:

1. Subtract atom footprints from owned blocker counts transactionally.
2. Reduce global blocker-change and newly-gas counts.
3. Initialize newly gas owned voxels as `ClosedVoid`.
4. Seed newly gas voxels that are reservoir sources or touch existing
   `OutsideAccessible` gas.
5. Drain local six-neighbor frontiers through `ClosedVoid` voxels.
6. Exchange only compact tangential offsets at crossed MPI faces and repeat
   until every rank is inactive.
7. Synchronize final face ghosts once before returning.

Removing one of several overlapping blockers updates its count without
running accessibility repair. A newly gas component without an accessible
seed remains `ClosedVoid`. A seeded component is promoted incrementally, and
ranks without local atom changes still join frontier exchange and termination
collectives when propagation reaches them.

`DistributedDesorptionUpdateResult` reports local/global occupancy changes,
sorted changed owned coordinates, opened/visited voxel counts, participating
ranks, communication rounds, and sent/received frontier entries. Construct the
updater with `DesorptionRepairMode::FullReclassification` for the full
distributed correctness-reference path.

The incremental and reference implementations are differentially tested on
one, two, four, and eight ranks. The worst correctness fixture opens more than
half the grid and propagates through every rank. Details are in
[`../benchmarks/PHASE_R3_DISTRIBUTED_DESORPTION.md`](../benchmarks/PHASE_R3_DISTRIBUTED_DESORPTION.md).

## Distributed atomic mixed update

`DistributedAtomChangeUpdater::apply_atom_changes()` accepts additions and
removals together:

```cpp
gasaccess::DistributedAtomChangeUpdater updater(precursor_radius);

const auto result = updater.apply_atom_changes(
    distributed_grid,
    {{added_atoms, added_atom_count},
     {removed_atoms, removed_atom_count}});
```

`AtomChangeEventBuffer` can own synchronized records across the collective
call. Record removed atoms before the KMC atom list discards or moves them:

```cpp
gasaccess::AtomChangeEventBuffer events;
events.record_deposition(added_atom);
events.record_desorption(removed_atom_at_old_position);

const auto result = updater.apply_atom_changes(
    distributed_grid,
    events.atom_changes());
events.clear();
```

The event buffer owns records but does not perform MPI synchronization. The
host must exchange both added and old removed-atom records using the same
excluded-radius coverage required for initialization.

Both atom views are aggregated before any blocker count is committed. A move
is therefore represented as removal at the old position and addition at the
new position without exposing an intermediate hole or blocker. If any rank
rejects the batch, successful local commits are rolled back before all ranks
receive an exception.

After final occupancy is known, the updater performs distributed closing
repair around newly solid voxels, then distributed opening repair from newly
gas voxels. The ordering is intentional: the closing pass may conservatively
close gas that is reachable only through a newly opened route, and the opening
pass restores that route using the final solid mask.

A preceding closing pass can make face-ghost labels stale. Rather than
exchanging complete ghost faces between passes, newly gas boundary voxels send
compact seed queries to neighboring owners. Neighbors acknowledge only those
adjacencies whose post-closing owned state remains `OutsideAccessible`. Normal
opening frontiers then propagate through `ClosedVoid` voxels. One full face
ghost exchange occurs after both passes.

`DistributedAtomChangeUpdateResult` reports the two occupancy directions and
separate closing/opening traversal, rank-participation, communication-round,
and frontier-entry metrics. `AtomChangeRepairMode::FullReclassification`
selects the full distributed reference. Phase R4 differential tests cover one,
two, four, and eight ranks; see
[`../benchmarks/PHASE_R4_MIXED_ATOM_CHANGES.md`](../benchmarks/PHASE_R4_MIXED_ATOM_CHANGES.md).
