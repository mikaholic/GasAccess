# Reversible Atom-Change Development Plan

Status: Complete; Phases R1-R6 implemented and verified
Last updated: 2026-09-05

## 1. Objective

Extend GasAccess from monotonic deposition-only updates to reversible atomic
changes while preserving communication-free accessibility queries. A single
collective update must accept batched deposited and desorbed atoms, commit the
final voxel occupancy once, repair only accessibility that may have become
stale, and return with consistent owned and ghost states on every MPI rank.

The target runtime sequence is:

```text
batched added and removed atoms
               |
               v
aggregate net blocker-count changes
               |
               v
commit final voxel occupancy once
               |
               v
detect newly solid and newly gas voxels
               |
               +--> closing repair: OutsideAccessible -> ClosedVoid
               |
               +--> opening repair: ClosedVoid -> OutsideAccessible
               |
               v
synchronize final face ghosts once
               |
               v
resume allocation-free cached queries
```

The extension has four acceptance goals:

1. Pure deposition remains correct and does not suffer a material performance
   regression.
2. Pure desorption incrementally opens reservoir-connected voids without a
   full-grid scan.
3. A mixed deposition/desorption batch is interpreted atomically from its final
   occupancy, not as two observable intermediate structures.
4. Deposition, desorption, and mixed repair have comparable repeated benchmark
   coverage for one, two, four, and eight MPI ranks.

## 2. Starting constraints

Before Phase R1, the implementation was intentionally monotonic:

- `GasGrid` and `DistributedGasGrid` store `GasState`, but do not retain how
  many atoms block each voxel.
- `AtomVoxelizer` only changes gas voxels to `Solid`.
- `DepositionUpdater::apply_deposition()` and its distributed counterpart only
  report newly solid voxels.
- Affected-region repair only handles connectivity loss by changing stale
  `OutsideAccessible` voxels to `ClosedVoid`.
- The documented fallback for a solid-to-gas change is complete occupancy
  reconstruction followed by full distributed classification.

Desorption cannot be implemented correctly by setting every voxel covered by a
removed atom to gas. Another atom may still cover the same voxel. Reversible
occupancy accounting is therefore a prerequisite for incremental opening
repair.

## 3. API and behavior contracts

### 3.1 Primary C++ interface

The primary update interface will represent one atomic batch:

```cpp
struct AtomChangeBatch {
    AtomView added_atoms{};
    AtomView removed_atoms{};  // old positions and radii
};

AtomChangeUpdateResult apply_atom_changes(
    GasGrid& grid,
    const AtomChangeBatch& changes);
```

The MPI equivalent is collective over the grid communicator:

```cpp
DistributedAtomChangeUpdateResult apply_atom_changes(
    DistributedGasGrid& grid,
    const AtomChangeBatch& changes);
```

Existing entry points remain as compatibility wrappers:

```cpp
apply_deposition(grid, deposited_atoms);
apply_desorption(grid, removed_atoms);
```

### 3.2 Atomic batch semantics

- Validate all inputs and all prospective blocker-count changes before
  modifying the grid.
- Aggregate additions and removals per voxel before committing occupancy.
- Derive transitions from the old and final blocker counts; do not expose or
  repair an addition-only or removal-only intermediate structure.
- Treat an atom move as removal at its old position plus addition at its new
  position in the same batch.
- Require every added or removed atom record exactly once. The initial design
  does not add an internal atom registry; tKMC remains the atomic source of
  truth.
- Removed atoms must be passed with their old positions and radii. Under MPI,
  their synchronized event records must remain available until the collective
  update returns.
- All ranks call the collective in the same order, including ranks with an
  empty local batch.

### 3.3 Result contract

The unified result will report at least:

```cpp
enum class AccessibilityRepairKind : std::uint8_t {
    None,
    Closing,
    Opening,
    Mixed,
    FullReclassification
};

std::uint64_t newly_solid_count;
std::uint64_t newly_gas_count;
std::uint64_t repair_closed_voxel_count;
std::uint64_t repair_opened_voxel_count;
std::uint64_t closing_visited_voxel_count;
std::uint64_t opening_visited_voxel_count;
AccessibilityRepairKind repair_kind;
```

Changed voxel coordinates remain sorted and unique. Distributed results also
retain rank participation, communication-round, and frontier-entry metrics for
the closing and opening passes.

## 4. Phase R1 — Reversible voxel occupancy

Status: complete (2026-09-05). Implementation and measurements are recorded in
[`../benchmarks/PHASE_R1_REVERSIBLE_OCCUPANCY.md`](../benchmarks/PHASE_R1_REVERSIBLE_OCCUPANCY.md).

### 4.1 Implementation

- Add a checked per-voxel blocker count, initially `std::uint32_t`.
- Keep blocker counts in a separate structure-of-arrays allocation from
  `GasState`. Accessibility queries and flood fills must not read blocker
  counts.
- Store blocker counts only for MPI-owned voxels. Face ghosts continue to carry
  final `GasState`; blocker-count ghosts are not required.
- During initial voxelization, increment the blocker count for every atom whose
  excluded sphere contains the voxel center.
- Define occupancy from the count:

```text
blocker_count > 0  -> Solid
blocker_count == 0 -> gas
```

- Aggregate signed deltas for a batch and update each touched owned voxel once
  where practical.
- Reject decrement underflow and increment overflow before committing the
  batch.
- Report only net occupancy transitions:

```text
0 -> positive  = newly solid
positive -> 0  = newly gas
all other changes leave occupancy unchanged
```

- Preserve direct state-setting fixture support by keeping blocker counts and
  `GasState::Solid` consistent in test/setup APIs.
- Reuse scratch buffers and avoid per-event full-grid clearing.

Likely implementation areas include `gas_grid.*`, `mpi_gas_grid.*`,
`atom_voxelizer.*`, a new atom-change description header, and replacement or
extension of the current deposition-only updater types.

### 4.2 Correctness tests

- Two atoms overlap one voxel; removing either one leaves the voxel solid.
- Removing the final blocker changes the voxel to gas.
- Added and removed footprints cancel within one batch.
- Added and removed footprints partially overlap.
- Periodic atom footprints update counts once per physical atom and owned
  voxel.
- Empty batches and batches with no net transition are no-ops.
- Invalid removal, count underflow, count overflow, and null non-empty views
  fail before mutation.
- Serial and distributed occupancy match full voxelization of the final atom
  structure.
- All existing deposition tests remain passing.

### 4.3 Memory and performance expectations

A `std::uint32_t` blocker array adds four bytes per owned voxel. A `128^3` grid
contains 2,097,152 voxels, so the total owned blocker storage is 8 MiB before
MPI distribution. The expected runtime effects are:

| Operation | Expected impact |
|---|---|
| Cached accessibility query | None beyond unrelated cache pressure |
| Existing flood-fill traversal | No direct blocker-count access |
| Initialization | Additional count update for each covered voxel |
| Deposition | Additional count update even when a voxel was already solid |
| Desorption | Enables local work instead of mandatory full reconstruction |

Capture Release-mode deposition baselines before R1, then repeat them after R1
on the same machine and build configuration. Use the existing timing method:
untimed warmups, enough repetitions to exceed one measured second, and the
average of the slowest-rank elapsed time from each repetition.

The following are investigation thresholds rather than correctness failures:

| Measurement | Investigate when regression exceeds |
|---|---:|
| Cached query latency | 3% |
| Existing repair traversal | 5% |
| End-to-end deposition update | 10% |
| Initialization | 15% |

Do not introduce compressed or sparse counters before measurement. If the
straightforward array fails the memory or performance review, evaluate checked
16-bit counts or sparse storage for overlap counts as separate optimizations.

### 4.4 Completion gate

- Occupancy after every tested batch exactly matches full voxelization from the
  final atoms.
- Existing deposition correctness remains unchanged.
- Before/after runtime and peak-memory results are recorded for MPI ranks
  1, 2, 4, and 8.
- Any regression beyond an investigation threshold is explained or corrected
  before R2.

## 5. Phase R2 — Serial desorption and opening repair

Status: complete (2026-09-05). Implementation and verification are recorded in
[`../benchmarks/PHASE_R2_SERIAL_DESORPTION.md`](../benchmarks/PHASE_R2_SERIAL_DESORPTION.md).

### 5.1 Algorithm

After final blocker counts are committed:

1. Initialize newly gas voxels as `ClosedVoid`.
2. Identify newly gas voxels that are reservoir sources or touch an existing
   `OutsideAccessible` voxel.
3. Seed a six-neighbor flood fill from those voxels.
4. Traverse connected `ClosedVoid` and newly gas voxels.
5. Relabel every reached voxel `OutsideAccessible`.
6. Leave newly opened components without a source connection as `ClosedVoid`.

The following fast paths avoid unnecessary traversal:

- No newly gas voxel: return.
- Newly gas component has no source or accessible neighbor: keep it closed.
- Newly gas voxels touch accessible gas but no old closed component: promote
  only the newly gas cluster.

Add an opening-repair implementation with reusable BFS storage and epoch
markers analogous to the current closing repair.

### 5.2 Tests

- Removal leaves occupancy unchanged because another atom overlaps.
- A newly created void remains closed.
- A one-voxel opening exposes a small cavity.
- A multi-voxel opening exposes multiple connected closed regions.
- Desorption exposes a configured source voxel directly.
- Openings cross each periodic seam.
- Empty, repeated, and invalid event batches.
- Deterministic randomized deposition/desorption sequences compared after every
  event with full serial voxelization and `ExteriorClassifier`.

### 5.3 Completion gate

Every serial incremental result, state count, changed-coordinate list, and
query result matches full reconstruction and classification.

## 6. Phase R3 — Distributed MPI desorption repair

Status: complete (2026-09-05). Implementation and verification are recorded in
[`../benchmarks/PHASE_R3_DISTRIBUTED_DESORPTION.md`](../benchmarks/PHASE_R3_DISTRIBUTED_DESORPTION.md).

### 6.1 Algorithm

1. Apply synchronized atom-count deltas only to owned voxels.
2. Reduce global blocker-count changes and newly gas counts.
3. Prepare opening seeds from owned source/accessible boundaries.
4. Drain each local frontier through owned `ClosedVoid` voxels.
5. Exchange compact tangential offsets only when the frontier crosses an MPI
   face.
6. Use a global activity reduction for termination.
7. Relabel reached owned voxels `OutsideAccessible`.
8. Exchange final face-ghost states once before returning.

Reuse the current decomposition, periodic normalization, face-message format,
epoch workspaces, and slowest-rank metric conventions where practical.

### 6.2 MPI correctness matrix

Run at 1, 2, 4, and 8 ranks:

- Local opening contained within one rank.
- Opening located on an MPI boundary.
- Periodic-seam opening.
- Cavity spanning several ranks.
- Large cavity whose opening propagation visits every rank.
- Ranks with no local atom change but required collective participation.
- Empty and no-net-change collective batches.
- Differential comparison among incremental MPI, forced-full MPI, and serial
  reference results.
- Final owned counts, face ghosts, changed coordinates, and accessibility
  queries.
- Balanced sent and received frontier-entry totals.

### 6.3 Completion gate

All rank counts match the serial and forced-full references. The worst fixture
must visit and open voxels on every rank and finish with balanced frontier
traffic.

## 7. Phase R4 — Mixed deposition and desorption

Status: complete (2026-09-05). Implementation and verification are recorded in
[`../benchmarks/PHASE_R4_MIXED_ATOM_CHANGES.md`](../benchmarks/PHASE_R4_MIXED_ATOM_CHANGES.md).

### 7.1 Safe first implementation

Expose mixed batches as soon as occupancy accounting is available, but route a
batch containing both newly solid and newly gas voxels through full distributed
classification. This establishes correct API semantics before optimizing the
mixed path.

### 7.2 Incremental implementation

After committing the final solid mask:

1. Run closing repair around newly solid voxels.
2. Run opening repair from newly gas/source-connected voxels.
3. Let the opening pass restore any region conservatively closed by the first
   pass but connected through a newly opened route in the final geometry.
4. Deduplicate all changed coordinates.
5. Synchronize face ghosts once after both passes.

The full-classification fallback remains available as a debug/reference mode.
Do not remove the default mixed fallback until the incremental two-pass result
passes differential testing.

### 7.3 Tests

- Addition and removal cancel with no net occupancy transition.
- Atom move represented as remove-old plus add-new.
- Close one cavity while opening a different cavity.
- Partially overlapping added and removed atom footprints.
- Closing and opening changes owned by different ranks.
- One local pass combined with one all-rank pass.
- Both passes cross periodic seams.
- Deterministically shuffled mixed-event sequences.
- Random mixed sequences compared after every batch against full serial and
  distributed classification.

### 7.4 Completion gate

The incremental mixed path becomes the default only when all tested states,
counts, changed coordinates, ghosts, and queries exactly match the full
reference for MPI ranks 1, 2, 4, and 8.

## 8. Phase R5 — Public interfaces and tKMC integration

Status: complete (2026-09-05). Interfaces and integration verification are
recorded in
[`../benchmarks/PHASE_R5_PUBLIC_INTEGRATION.md`](../benchmarks/PHASE_R5_PUBLIC_INTEGRATION.md).

### 8.1 Library interfaces

- Stabilize the C++ `apply_atom_changes()` and `apply_desorption()` entry points.
- Keep `apply_deposition()` source-compatible as a wrapper.
- Add equivalent C API batch descriptors, result fields, and lifecycle
  functions.
- Preserve existing C deposition symbols and behavior.
- Update documentation from the monotonic-only contract to the reversible
  contract after R4 acceptance.
- Extend the SPPARKS/tKMC-style atom buffer or add an event buffer that retains
  old removed-atom records through the collective call.

### 8.2 tKMC call sequence

At each caller-selected synchronization point:

1. Complete the KMC sector's atomic events.
2. Synchronize owned and ghost added/removed event records far enough to cover
   every affected owned voxel.
3. Preserve removed atoms' old positions and radii.
4. Call `apply_atom_changes()` collectively once.
5. Resume communication-free accessibility queries after the call returns.

Add a mock KMC integration sequence containing deposition-only,
desorption-only, atom movement, mixed, empty, and no-net-change batches.

### 8.3 Completion gate

Existing deposition clients compile unchanged, while the mock tKMC path uses
one collective batch API for all supported atomic changes.

## 9. Phase R6 — Efficiency and scaling benchmarks

Completion: implemented and verified on 2026-09-05. All 48 repeated Release
configurations passed, and results are archived in
[`../benchmarks/PHASE_R6_REVERSIBLE_EFFICIENCY.md`](../benchmarks/PHASE_R6_REVERSIBLE_EFFICIENCY.md)
and
[`../benchmarks/PHASE_R6_RESULTS.csv`](../benchmarks/PHASE_R6_RESULTS.csv).

Extend `gasaccess_mpi_efficiency_driver` and the existing repair fixture rather
than introducing a separate timing program. Preserve the current deposition
CLI as the default and add:

```text
--operation repair
--change-kind deposition|desorption|mixed
--case baseline|best|medium|worst
```

### 9.1 Fixture definitions

#### Deposition

Retain the existing baseline, best, medium, and worst fixtures unchanged as
regression measurements.

#### Desorption

| Case | Geometry and event | Expected work |
|---|---|---|
| `baseline` | Remove one of two overlapping blockers | No solid-to-gas transition |
| `best` | Open an approximately `8x8x8` pocket on one rank | Small local promotion |
| `medium` | Open a cavity containing approximately 25% of the grid | Multi-rank opening flood |
| `worst` | Open a cavity containing more than half of the grid | Every rank visits and opens voxels |

The worst fixture must use a large closed cavity connected to the reservoir by
a removable plug. Removing the plug should propagate through all MPI ranks and
open a substantial, but not necessarily total, portion of the grid.

#### Mixed change

Use two disjoint cavities so the same batch exercises both state directions.

| Case | Geometry and event | Expected work |
|---|---|---|
| `baseline` | Added and removed footprints cancel | Detection only |
| `best` | Close one small pocket and open another on one rank | Small local closing and opening |
| `medium` | Close and open components totaling approximately 25% of the grid | Moderate distributed work |
| `worst` | Close one large component and open another; combined affected volume exceeds half the grid | Both passes involve every rank |

### 9.2 Benchmark matrix

Run the complete strong-scaling matrix:

```text
3 change kinds
x 4 cases
x 4 MPI rank counts (1, 2, 4, 8)
= 48 benchmark configurations
```

Use identical global grid dimensions for each rank count. The x dimension must
remain divisible by the process count, and the worst desorption and mixed
fixtures must cross every x-decomposed rank.

### 9.3 Timing method

- Run at least three untimed warmups by default.
- Run at least ten measured repetitions.
- Continue until cumulative measured time reaches at least one second.
- Time only the complete collective `apply_atom_changes()` call.
- Exclude fixture construction, initial classification, state restoration, and
  correctness validation from the timed interval.
- Use the maximum elapsed rank time for each repetition.
- Report the arithmetic mean of those maximum-rank samples as the primary
  value, together with minimum, maximum, standard deviation, repetition count,
  and cumulative measured time.
- Register short one-repetition CTest smoke cases, while release benchmark runs
  retain the repeated timing defaults.

### 9.4 Reported metrics

- End-to-end update time.
- Occupancy-update time when phase timing is enabled.
- Newly solid and newly gas voxel counts.
- Closing/opening visited voxel counts.
- Newly closed and newly opened voxel counts.
- Closing/opening participating-rank counts.
- Minimum and maximum per-rank traversal work.
- Communication rounds.
- Sent and received frontier entries.
- Time per visited voxel.
- Peak resident memory.
- Paired incremental versus forced-full reclassification time and speedup.

### 9.5 Benchmark correctness conditions

Timings are measurements, not fixed pass/fail thresholds. Each benchmark must
still assert:

- Incremental states and counts match forced full classification.
- Expected reservoir, cavity, and control probes have correct accessibility.
- Changed coordinates are sorted and unique.
- Sent and received frontier-entry totals balance.
- The best case remains local to one rank when the fixture permits it.
- Every rank participates in the relevant worst-case traversal.
- The worst desorption case opens more than half the grid.
- The worst mixed case performs both closing and opening repair.

### 9.6 Completion gate

- All 48 configurations pass correctness checks.
- Repeated Release-mode results are archived in a new benchmark report under
  `docs/benchmarks/`.
- Deposition before/after R1 comparisons are included.
- Desorption and mixed incremental/full speedups are summarized without
  claiming an advantage when the affected size approaches the full grid.

## 10. Planned implementation order

```text
R1 reversible occupancy and performance gate
 |
 v
R2 serial opening repair
 |
 v
R3 distributed opening repair
 |
 v
R4 mixed fallback, then incremental mixed repair
 |
 v
R5 stable public and tKMC integration interfaces
 |
 v
R6 complete benchmark matrix and report
```

Do not optimize mixed repair before serial and distributed desorption are
differentially correct. Do not interpret a nearly exhausted benchmark budget
or a fast smoke run as performance acceptance; release measurements must use
the repeated timing protocol above.
