# Phase R5 public interfaces and KMC integration

Date: 2026-09-05

Phase R5 promotes the reversible atom-change path from an internal C++
algorithm to an application-facing integration contract. It retains the
adsorption-specific API, adds reversible C interfaces, and verifies a realistic
caller-controlled MPI event sequence.

## Stable C++ entry points

`AtomChangeUpdater` and `DistributedAtomChangeUpdater` retain
`apply_atom_changes()` as the primary operation. Both now also expose
`apply_adsorption()` and `apply_desorption()` convenience methods so one
persistent updater can serve every event kind and reuse its workspaces.

The established `AdsorptionUpdater`, `DesorptionUpdater`,
`DistributedAdsorptionUpdater`, and `DistributedDesorptionUpdater` classes and
method signatures remain unchanged. Existing adsorption clients therefore
compile and retain their original result contracts and topology-filter path.

## Owning event buffer

`AtomChangeEventBuffer` owns separate vectors for added and removed atoms. It
supports adsorption, desorption, movement, and batch append operations. A move
records the old atom in the removal direction and the new atom in the addition
direction.

The buffer's `atom_changes()` view stays valid until the buffer is modified or
destroyed. This prevents removed records from disappearing when KMC erases or
moves atoms in its live structure before the accessibility collective begins.
The buffer intentionally performs no MPI communication and does not dictate
when KMC synchronizes.

## C99 interface

The public C header adds:

- `ga_atom_view` and `ga_atom_change_batch` descriptors;
- incremental/full atom-change repair modes and repair-kind constants;
- `ga_apply_desorption()` and `ga_apply_desorption_with_mode()`;
- `ga_apply_atom_changes()` and `ga_apply_atom_changes_with_mode()`; and
- creation-by-update, summary access, changed-voxel access, and destruction
  functions for the opaque `ga_atom_change_result` handle.

`ga_atom_change_update_summary` reports blocker-count changes, newly solid and
newly gas counts, repair kind, closing/opening traversal and relabel counts,
method flags, final changed-voxel count, and final classification.

The reversible result is separate from `ga_update_result`. No existing C type,
constant, function signature, adsorption behavior, or result accessor was
changed. The per-grid C++ updater is cached by precursor radius and repair mode
so repeated C calls reuse repair workspaces.

## KMC synchronization contract

For each caller-selected synchronization point:

1. KMC completes its conflict-free sector events.
2. It retains added atoms at new positions and removed atoms at old positions.
3. It synchronizes both event directions far enough to cover every owned voxel
   within `R_atom + R_precursor`.
4. Every rank calls `apply_atom_changes()` once, including ranks with an empty
   local buffer.
5. The event buffer can be cleared after the collective returns.
6. KMC resumes communication-free cached accessibility queries.

## Mock integration acceptance

The SPPARKS-style MPI mock now runs a six-step reversible sequence:

1. adsorption-only;
2. desorption-only;
3. movement represented as remove-old plus add-new;
4. a mixed adsorption/desorption batch;
5. an empty collective batch; and
6. an exact no-net-change batch.

Added and removed events are synchronized independently through mock atom
halos. Event atoms have a nonzero excluded radius and are positioned so their
footprints cross MPI ownership boundaries. The removed-event application
retains old records separately, and the owning GasAccess buffer outlives those
temporary synchronization objects.

After every step, the test compares the distributed result with serial
incremental repair and then reconstructs occupancy and accessibility from the
current final atom set. It checks global transition counts, final changed
coordinates, every owned blocker count and state, face ghosts, and cached
coordinate queries.

## Correctness coverage

- The pure-C client now contains six test groups, including mixed movement,
  pure desorption, exact cancellation, full-reference mode, result accessors,
  invalid modes, null batches, and invalid removed-atom views.
- The serial atom-change executable contains eight groups, including ownership
  of copied old/new records and all three unified updater entry points.
- The SPPARKS-style mock contains three groups and runs the reversible event
  sequence at 1, 2, 4, and 8 MPI ranks.
- Existing adsorption-specific test coverage remains compiled and executed.

## Verification

The Release MPI-enabled build is warning-clean. All 70 registered CTest cases
pass, including the C99 client and the mock KMC integration matrix at 1, 2,
4, and 8 ranks.

## Scope boundary

The production KMC source is outside this repository, so choosing the exact
sector synchronization call site and mapping its atom storage are application
handoff tasks. The subsequent Phase R6 added the repeated adsorption,
desorption, and mixed efficiency/scaling benchmarks; its results are in
[`PHASE_R6_REVERSIBLE_EFFICIENCY.md`](PHASE_R6_REVERSIBLE_EFFICIENCY.md).

## Conclusion

Phase R5 passes its completion gate: existing adsorption interfaces remain
compatible, reversible serial C and C++ interfaces are available, and a mock
KMC path uses one collective update contract for every required event type at
all requested MPI rank counts.
