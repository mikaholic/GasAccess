# Phase R2 serial desorption and opening repair

Date: 2026-09-05

Phase R2 adds incremental serial accessibility repair after atom removal. It
uses the blocker counts introduced in Phase R1, so removing one atom does not
free a voxel that is still covered by another atom.

## Implemented path

`DesorptionUpdater::apply_desorption()` performs the following sequence:

1. Require a fully classified input grid.
2. Transactionally subtract all removed-atom footprints from blocker counts.
3. Detect only `positive -> 0` blocker transitions as newly gas voxels.
4. Initialize newly gas voxels as `ClosedVoid`.
5. Seed newly gas voxels that are reservoir sources or touch existing
   `OutsideAccessible` gas.
6. Flood through connected `ClosedVoid` voxels and relabel them
   `OutsideAccessible`.
7. Leave unseeded components closed and return sorted, unique state changes.

`OpeningRegionRepair` retains its frontier, seed, and epoch-marker storage
between calls. It does not clear or scan a full-grid visited array per event.
The query path is unchanged and still reads only cached gas states.

`DesorptionRepairMode::FullReclassification` is available as a serial
correctness reference. It is not the default.

## Fast paths

- Empty removal batch: no occupancy or accessibility work.
- Overlap-only removal (`2 -> 1` blockers): update one count and skip repair.
- Newly gas component with no source or accessible neighbor: keep it
  `ClosedVoid` without flood filling.
- Newly gas voxel surrounded by accessible gas: visit and open only that voxel.

In the locality fixture, opening an 11-voxel cavity in a 10,000-voxel grid
visits exactly 11 voxels. This verifies that the incremental path follows the
affected component rather than scanning the complete grid.

## Correctness coverage

The Phase R2 test executable covers:

- fully classified-grid and atom-view validation;
- empty, invalid, underflow, and repeated-removal behavior;
- removal of one versus the final overlapping blocker;
- creation of a new closed void with no reservoir path;
- a one-voxel opening that exposes a cavity;
- one batch opening two disconnected cavities;
- direct exposure of an explicit reservoir source;
- opening propagation across each periodic x, y, and z seam;
- the local 11-of-10,000-voxel traversal fixture; and
- 160 deterministic deposition/desorption events compared after every event
  with fresh atom voxelization and full `ExteriorClassifier` reconstruction.

The randomized differential check compares blocker counts, every gas state,
all four state counts, sorted changed-voxel IDs, classification summaries, and
cached accessibility-query results.

The final Release MPI-enabled build passes all 61 registered CTest cases. This
includes the unchanged deposition tests and the existing 1/2/4/8-rank MPI and
efficiency matrix.

## Scope boundary

Phase R2 is deliberately serial. It does not yet provide distributed opening
frontiers, MPI desorption synchronization, incremental mixed
deposition/desorption repair, or C API exposure. Distributed desorption repair
is Phase R3; mixed-event repair follows in Phase R4.

## Conclusion

Phase R2 passes its completion gate: every tested incremental serial result
matches complete reconstruction and classification, while local and
overlap-only events avoid full-grid connectivity work.
