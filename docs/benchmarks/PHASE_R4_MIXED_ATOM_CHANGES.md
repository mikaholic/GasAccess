# Phase R4 atomic mixed atom changes

Date: 2026-09-05

Phase R4 adds one serial and one collective MPI operation for applying added
and removed atoms as a single atomic occupancy transaction. It supports pure
adsorption, pure desorption, atom movement, and genuinely mixed batches while
preserving the final-state semantics required by tKMC.

## Public C++ path

The serial entry point is
`AtomChangeUpdater::apply_atom_changes()`. The distributed entry point is
`DistributedAtomChangeUpdater::apply_atom_changes()` and is collective over
the distributed grid communicator.

Each accepts an `AtomChangeBatch` containing added atoms at their new
positions and removed atoms at their old positions. A move is represented by
one record in each view. Both footprints are aggregated before blocker counts
are committed, so connectivity repair never observes an addition-only or
removal-only intermediate geometry.

The existing adsorption and desorption updaters remain source-compatible.
`AtomChangeRepairMode::FullReclassification` provides a serial or distributed
correctness-reference path for differential testing and debugging.

## Incremental repair sequence

After the final blocker counts have been committed, the incremental updater:

1. identifies net `gas -> solid` and `solid -> gas` transitions;
2. runs closing repair from the newly solid voxels;
3. runs opening repair from the newly gas voxels using the final solid mask;
4. reports sorted, unique voxels whose final state differs from its state
   before the batch; and
5. synchronizes MPI face-ghost states once after both repair passes.

Closing first is conservative. It may temporarily close a component whose old
route was blocked by the batch. The following opening pass restores any part
that remains connected through a newly opened route in the final geometry.
Voxels closed and reopened during these two passes are excluded from the final
changed-coordinate list when their final state did not change.

## Distributed communication

The distributed closing repair visits owned accessible voxels and exchanges
only compact face-frontier offsets when propagation crosses an ownership
boundary. The opening repair uses the same ownership and frontier rules.

Closing can make neighboring ghost labels stale before opening begins. A
mixed batch therefore performs a compact boundary seed query and
acknowledgment for newly gas voxels instead of exchanging complete ghost
faces between the passes. The neighboring owner answers from its post-closing
owned state. A single full face-state exchange establishes final ghost
consistency after both passes.

If any rank rejects an occupancy update, ranks that committed locally restore
their old blocker counts and states before the collective operation throws.
This prevents a partially committed distributed batch.

## Result metrics

The unified results report:

- local and global blocker-count, newly solid, and newly gas counts;
- the repair kind (`None`, `Closing`, `Opening`, `Mixed`, or
  `FullReclassification`);
- sorted, unique final changed voxel IDs or owned coordinates;
- separate closing and opening visits and changed-voxel counts; and
- separate distributed rank-participation, communication-round, and
  sent/received frontier-entry counts.

## Correctness coverage

The serial executable contains seven test groups covering:

- validation, no-op batches, exact cancellation, and failed-removal rollback;
- an atom move applied without exposing an intermediate structure;
- closing one cavity while opening another;
- partially overlapping added and removed footprints;
- a channel swap that closes and reopens transient voxels but reports only
  final changes;
- pure adsorption and pure desorption through the unified API; and
- 120 deterministic mixed batches checked after every event against both
  forced full reclassification and reconstruction from the final atom set.

The MPI executable runs eight test groups at 1, 2, 4, and 8 ranks:

- empty batches, cancellation, and a cross-rank move;
- collective rejection with rollback of successful local commits;
- closing and opening on different ranks;
- post-closing boundary seed validation, including the stale-ghost case;
- a channel swap whose transiently closed component is reopened;
- local closing combined with opening work that reaches every rank;
- closing and opening propagation across a periodic x seam; and
- a deterministic sequence of 24 mixed events with updater workspace reuse.

Every valid MPI event is compared with serial incremental and distributed
forced-full references. The checks include blocker counts, all owned states,
face ghosts, state totals, coordinate queries, final changed-coordinate lists,
repair counts, participating ranks, and globally balanced frontier traffic.

## Verification

The Release MPI-enabled build is warning-clean. All 70 registered CTest cases
pass, including the mixed atom-change executable at 1, 2, 4, and 8 ranks. The
complete suite also confirms that existing initialization, query, adsorption,
desorption, C API, static-integration, and efficiency tests remain passing.

## Scope boundary

Phase R4 establishes the C++ algorithm and correctness contract. Stable C API
descriptors and the tKMC-oriented event buffer were deferred to Phase R5 and
are now complete. Repeated adsorption, desorption, and mixed
best/medium/worst efficiency measurements remain Phase R6.

## Conclusion

Phase R4 passes its completion gate: incremental mixed updates match final
atom-set reconstruction and full serial/distributed classification, including
ownership boundaries and periodic seams, for MPI ranks 1, 2, 4, and 8.
