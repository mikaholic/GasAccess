# Phase 6 serial reference baseline

Date: 2026-08-12

Code baseline: commit `ac2bbf9` (before the Phase 7 incremental filter). Current
development builds may produce faster update timings while retaining these
state checksums. The Phase 7 comparison is recorded in
[`PHASE7_FILTER.md`](PHASE7_FILTER.md).

This record establishes a reproducible correctness and turnaround-time baseline
before incremental connectivity algorithms are introduced. It uses generated
structures because the production atomic-structure format and representative
input have not yet been supplied.

## Environment

- Build: CMake `Release`, C++17, warnings treated as errors
- Compiler: GCC 8.5.0
- CMake: 3.26.5
- OS: Linux 6.6.87.2 under WSL2, x86-64
- CPU: Intel Core i7-13700, 24 logical CPUs
- Execution: one process and one thread; no MPI

These are single-run reference numbers, not statistically controlled
performance claims. Timing changes should be compared on the same machine with
an otherwise idle system.

## Commands

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
    -DCMAKE_C_FLAGS=-Werror -DCMAKE_CXX_FLAGS=-Werror
cmake --build build --parallel

./build/gasaccess_reference_driver \
    --scenario open-trench --query-count 1000000 --update-count 3

./build/gasaccess_reference_driver \
    --scenario sealed-trench --query-count 1000000 --update-count 3

./build/gasaccess_reference_driver \
    --scenario bulk --nx 128 --ny 64 --nz 128 \
    --atom-count 1000000 --query-count 1000000 --update-count 3 \
    --seed 12345
```

All three runs used spacing `1.0`, atom radius `0.2`, precursor radius `0.25`,
and non-periodic boundaries with the high-z face as the reservoir. The generated
atom centers lie on voxel centers, so the combined radius `0.45` blocks the
intended occupied voxel without reaching its face-neighbor centers.

## Correctness results

| Scenario | Atom records | Voxels | Solid | Outside | Closed | Initial checksum |
|---|---:|---:|---:|---:|---:|---|
| Open trench | 4,032 | 131,072 | 4,032 | 127,040 | 0 | `0x58555ea38df5f565` |
| Sealed trench | 4,672 | 131,072 | 4,672 | 86,720 | 39,680 | `0x476e87e5f025d7e5` |
| Bulk, one million atoms | 1,000,000 | 1,048,576 | 644,130 | 348,817 | 55,629 | `0x34a6ae79e2f13ec8` |

The open trench has no closed gas. Adding the deterministic roof produces a
39,680-voxel closed region. The bulk generator permits repeated atom locations,
which is why one million atom records produce 644,130 unique solid voxels.

## Timing results

All times are milliseconds. The query column is the total for one million
cached position queries. The update column is the total for three one-atom
Phase 5 updates; each geometry-changing update intentionally performs a full
snapshot, global serial reclassification, and state diff.

| Scenario | Atom generation | Voxelization | Classification | Queries | Queries/s | 3 full updates |
|---|---:|---:|---:|---:|---:|---:|
| Open trench | 0.073 | 0.245 | 3.892 | 22.242 | 44.96 million | 12.417 |
| Sealed trench | 0.132 | 0.261 | 2.874 | 35.385 | 28.26 million | 9.830 |
| Bulk, one million atoms | 19.761 | 69.656 | 23.727 | 52.739 | 18.96 million | 82.267 |

## Memory results

| Scenario | Persistent state | Frontier payload bound | Process peak RSS |
|---|---:|---:|---:|
| Open trench | 131,072 B | 1,048,576 B | 5,267,456 B |
| Sealed trench | 131,072 B | 1,048,576 B | 4,956,160 B |
| Bulk, one million atoms | 1,048,576 B | 8,388,608 B | 43,835,392 B |

Persistent gas state is exactly one byte per voxel. The frontier number is a
conservative `voxel_count * sizeof(VoxelId)` bound on the logical identifier
payload of the current full flood-fill frontier, not its measured allocation or
live occupancy. A vector allocator may reserve more capacity. Linux `VmHWM`
supplies peak RSS; it includes atom storage, traversal and update temporaries,
allocator/runtime overhead, and the executable itself, so it must not be
interpreted as gas-state storage alone.

## Scope remaining

This baseline does not read the production atomic-structure format and does not
exercise MPI. A format adapter will be added when the structure format and a
representative sample are available. Incremental timing improvements begin in
Phase 7; Phase 6 deliberately measures the full-recomputation reference path.
