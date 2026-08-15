# Phase 7 local-filter comparison

Date: 2026-08-14

This comparison repeats the Phase 6 million-atom synthetic workload after the
conservative local topology filter was integrated. It is a single-run sanity
measurement on the same environment recorded in
[`PHASE6_BASELINE.md`](PHASE6_BASELINE.md), not a statistically controlled
performance claim.

## Command

```sh
./build/gasaccess_reference_driver \
    --scenario bulk --nx 128 --ny 64 --nz 128 \
    --atom-count 1000000 --query-count 1000000 --update-count 3 \
    --seed 12345
```

## Result

| Measurement | Phase 6 reference | Phase 7 filter |
|---|---:|---:|
| Atom records | 1,000,000 | 1,000,000 |
| Voxels | 1,048,576 | 1,048,576 |
| Initial checksum | `0x34a6ae79e2f13ec8` | `0x34a6ae79e2f13ec8` |
| Final checksum | `0x32609a43f1cd6264` | `0x32609a43f1cd6264` |
| Geometry-changing updates | 3 | 3 |
| Full reclassifications | 3 | 0 |
| Total update time | 82.267 ms | 24.625 ms |
| Peak RSS | 43,835,392 B | 40,214,528 B |

The three updates were individually proven safe, reducing their measured total
time by about 70% (`3.34x`). The state checksums remained identical to the full
reference path.

The safe path still snapshots and scans the complete state array to discover
newly solid voxels and construct the public update result. Consequently, this
is not yet a fully local-complexity update. Later phases will replace those
transitional full-array operations; Phase 7 specifically establishes the
correct conservative topology decision.
