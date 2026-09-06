# Phase R1 reversible voxel occupancy

Date: 2026-09-05

Phase R1 adds reversible occupancy accounting without yet adding desorption
accessibility repair. `GasGrid` now stores a checked `std::uint32_t` blocker
count for every voxel. `DistributedGasGrid` stores the same count only for
owned voxels; face ghosts continue to store only `GasState`.

## Implemented behavior

- Existing serial and distributed adsorption entry points remain available.
- Initial and adsorption voxelization increments every covered voxel, including
  overlap with voxels that were already solid.
- `AtomChangeBatch` accepts added and removed atom views in one transaction.
- Additions and removals are aggregated by voxel before mutation.
- Only net `0 -> positive` and `positive -> 0` occupancy transitions are
  reported as newly solid and newly gas.
- Invalid atom data, null non-empty views, blocker underflow, and blocker
  overflow leave the grid and caller-owned change output unchanged.
- A final removal changes the affected voxel to `Unclassified`. Opening and
  reclassification of accessibility are intentionally deferred to Phase R2.
- Direct state setters keep synthetic fixtures consistent: setting `Solid`
  establishes one blocker, while setting a non-solid state clears blockers.

The low-level mixed-change methods are
`AtomVoxelizer::apply_atom_changes()` and
`DistributedGasGrid::apply_owned_atom_changes()`. The distributed method is an
owned-voxel transaction; the collective high-level update and opening repair
remain Phase R3/R4 work.

## Correctness coverage

The added tests cover:

- one and multiple blockers per voxel;
- removing one overlap versus removing the final blocker;
- exact cancellation within a mixed batch;
- partially overlapping old/new atom footprints;
- periodic images updating one physical owned voxel;
- empty batches and no-net-change batches;
- null views, invalid geometry, underflow, overflow, and rollback;
- randomized blocker counts against brute-force voxel coverage; and
- serial/distributed state and blocker-count parity at 1, 2, 4, and 8 ranks.

After implementation, all 60 registered CTest cases passed in the Release MPI
build, including every 1/2/4/8-rank test.

## Performance method

The comparison used two Release builds from the same source tree:

- before: commit `0bf9e77`;
- after: the Phase R1 working tree;
- compiler: GCC 8.5.0 with `-O3 -DNDEBUG`;
- MPI: Open MPI/OpenRTE 4.1.1;
- host: WSL2 on an Intel Core i7-13700, 24 logical CPUs.

Both builds used the `128x128x128` fixtures and the same commands. Initialization
and query ran for at least one measured second. Detection-only repair also ran
for at least one measured second. Best repair used 50 measured repetitions;
medium and worst repair used at least 10 repetitions and at least one measured
second when 10 repetitions were shorter than that. Reported times are the
average slowest-rank times.

## Initialization and cached query results

| MPI ranks | Initialization before (ms) | after (ms) | change | Query before (ns) | after (ns) | change |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 136.33 | 146.19 | +7.23% | 39.72 | 40.12 | +1.01% |
| 2 | 59.77 | 64.44 | +7.82% | 39.54 | 39.26 | -0.71% |
| 4 | 32.17 | 34.90 | +8.48% | 39.86 | 39.93 | +0.18% |
| 8 | 19.69 | 20.69 | +5.07% | 40.68 | 40.49 | -0.47% |

Initialization includes the required count writes for every atom/voxel overlap.
The maximum measured increase, 8.48%, is below the 15% investigation threshold.
Cached query latency changes by at most 1.01%, below its 3% threshold; the query
path does not read blocker counts.

## Existing adsorption-repair results

| MPI ranks | Detection before/after (us) | Best before/after (ms) | Medium before/after (ms) | Worst before/after (ms) |
|---:|---:|---:|---:|---:|
| 1 | 0.079 / 0.085 | 4.218 / 4.297 | 238.43 / 233.84 | 354.95 / 366.92 |
| 2 | 0.423 / 0.420 | 2.354 / 2.327 | 122.67 / 123.66 | 230.72 / 238.79 |
| 4 | 0.945 / 0.940 | 1.292 / 1.284* | 69.74 / 70.83 | 140.59 / 142.12 |
| 8 | 1.345 / 1.383 | 0.836 / 0.874 | 45.30 / 45.50 | 101.21 / 102.20 |

The largest stable increase is 7.59% for the one-rank detection-only path,
which is far below one microsecond. Best repair changes by -1.15% to +4.54%.
Medium and worst repair change by -1.93% to +3.50%. These remain below the 10%
end-to-end adsorption threshold, and substantial flood-fill cases remain below
the 5% repair-traversal threshold.

`*` The four-rank best fixture performs work on only one rank and showed a
bimodal 0.6/1.3 ms distribution under WSL scheduling for both binaries. Five
independent launches per build gave median launch averages of 1.292 ms before
and 1.284 ms after. The medium and worst all-rank cases did not show this mode.

## Memory result

The blocker array has a deterministic aggregate payload of 8 MiB for a 128^3
grid (`2,097,152 * 4` bytes), divided among MPI owners. The medium repair
fixture produced these summed process peak-RSS measurements:

| MPI ranks | Before (MiB) | After (MiB) | Increase (MiB) |
|---:|---:|---:|---:|
| 1 | 195.46 | 203.52 | 8.06 |
| 2 | 174.72 | 182.53 | 7.81 |
| 4 | 182.52 | 190.10 | 7.58 |
| 8 | 241.91 | 250.86 | 8.95 |

Peak RSS contains allocator and process-launch noise, but the measured
7.6-9.0 MiB range is consistent with the expected 8 MiB payload.

## Phase R1 conclusion

Phase R1 passes its completion gate. Reversible, overlap-safe occupancy is in
place; existing adsorption behavior remains correct; query latency is
effectively unchanged; initialization remains within budget; and meaningful
repair workloads remain within the planned performance thresholds. Phase R2
can build serial opening repair on the newly-gas transition output.
