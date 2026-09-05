# Phase R3 distributed desorption and opening repair

Date: 2026-09-05

Phase R3 extends the serial desorption path from Phase R2 across MPI ranks.
It preserves owned-only occupancy mutation and communication-free cached
queries while allowing a newly opened gas component to propagate through any
number of subdomains.

## Implemented path

`DistributedDesorptionUpdater::apply_desorption()` is collective over the gas
grid communicator and performs the following sequence:

1. Require a fully classified distributed grid on every rank.
2. Transactionally subtract synchronized removed-atom footprints from owned
   blocker counts.
3. Reduce blocker-count changes and `positive -> 0` occupancy transitions.
4. Return immediately for empty and overlap-only batches.
5. Initialize newly gas owned voxels as `ClosedVoid`.
6. Seed newly gas voxels that are reservoir sources or touch the existing
   `OutsideAccessible` region.
7. Flood locally through owned `ClosedVoid` voxels and exchange compact
   tangential offsets when the frontier crosses an MPI face.
8. Repeat until a global activity reduction says every frontier is empty.
9. Exchange final face-ghost states once before returning.

`DistributedOpeningRegionRepair` retains epoch markers, frontier storage, and
face buffers between calls. It does not clear a full owned-voxel marker array
for each event and does not gather or replicate the gas grid.

If any rank rejects an occupancy batch, successful ranks restore their prior
blocker counts and gas states before every rank receives an exception. This
prevents partial distributed occupancy commits and collective deadlock for
checked input failures.

`DesorptionRepairMode::FullReclassification` remains available as the MPI
correctness reference. It applies the same occupancy change and then invokes
the full `DistributedExteriorClassifier`.

## Result metrics

The distributed result reports:

- local and global blocker-count changes;
- local and global newly gas counts;
- sorted, unique changed owned coordinates;
- local opening visits and opened voxels;
- the global number of ranks participating in opening traversal;
- communication rounds and local sent/received frontier entries; and
- final local classification counts.

## Correctness coverage

The Phase R3 executable runs the same seven groups at 1, 2, 4, and 8 ranks:

- collective input rejection and rollback of rank-local successful changes;
- empty removal, overlap `2 -> 1`, and final-blocker `1 -> 0` behavior;
- an enclosed newly gas voxel that correctly remains `ClosedVoid`;
- a one-voxel local opening confined to one rank;
- an opening located on an MPI ownership boundary;
- opening through the periodic x seam;
- a large opening whose flood fill visits every rank; and
- a deterministic sequence of 16 removals with updater-workspace reuse.

After every valid event, the test compares incremental MPI, forced-full MPI,
and serial incremental results. It checks blocker counts, every owned state,
face ghosts, all state counts, cached coordinate queries, changed-coordinate
lists, opening counts, participating ranks, and globally balanced sent and
received frontier traffic.

The worst fixture uses a `32 x 8 x 8` grid. Removing one plug opens 1,537 of
2,048 voxels, so more than half the grid is updated. The opening traversal has
nonzero work on every tested rank and crosses every x-decomposed subdomain.

## Verification

The Release MPI-enabled build is warning-clean. All 65 registered CTest cases
pass, including the new distributed desorption executable at 1, 2, 4, and 8
ranks. The Open MPI runtime emits ignorable component-loader warnings for
unavailable `libfabric` transports on this WSL host and falls back to the
available local transport.

## Scope boundary

Phase R3 implements pure desorption correctness and traversal metrics. It does
not yet combine deposition and desorption in one atomic public update, expose
desorption through the C API, or add repeated desorption timing fixtures.
Those remain Phases R4, R5, and R6 respectively.

## Conclusion

Phase R3 passes its completion gate: incremental distributed states match the
serial and forced-full references at every tested rank count, the dominant
cavity opens on every rank, final ghosts are consistent, and global frontier
traffic balances exactly.
