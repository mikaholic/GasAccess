# Phase 8 affected-region repair verification

Date: 2026-08-14

Phase 8 replaces the Phase 7 full-classification fallback with serial searches
limited to components adjacent to newly solid voxels. This record focuses on
correctness and locality; statistically controlled latency measurements remain
part of Phase 10.

## Local pinch-off fixture

The deterministic locality test uses a `9x1x9` grid with a high-z reservoir,
two solid trench walls, and a three-voxel roof deposition.

| Measurement | Result |
|---|---:|
| Total voxels | 81 |
| Newly solid roof voxels | 3 |
| Newly closed cavity voxels | 21 |
| Repair BFS visits | 21 |
| Total changed voxels | 24 |
| Full-reference state differences | 0 |

The repair completely explores the 21-voxel source-disconnected cavity. The
exterior-side seeds are themselves reservoir voxels, so the repair does not
scan the remaining accessible volume.

## Other deterministic coverage

- Removing the only source from a five-voxel line visits and closes all four
  surviving gas voxels.
- Removing one of two sources preserves the connected gas component.
- A two-voxel cut in an eleven-voxel line closes only the three-voxel middle
  component.
- A periodic-seam cut produces the same states as full classification.
- The incremental and forced-reference modes return identical summaries,
  changed IDs, and voxel states.

## Large synthetic regression

The Phase 6 million-atom workload retains its exact state checksums:

```text
initial_state_checksum=0x34a6ae79e2f13ec8
final_state_checksum=0x32609a43f1cd6264
```

Its three sample updates are locally safe Phase 7 events, so they correctly
perform neither affected-region repair nor full reclassification. The driver
now reports those method counts and repair visit/closure counts explicitly.

## Remaining limitation

`DepositionUpdater` still takes a full state snapshot and scans it to discover
newly solid voxels and build the result. Phase 8 removes global connectivity
traversal for pinch-off repair; it does not yet make the entire deposition
pipeline local-complexity.
