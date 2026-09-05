# Phase 15 lifecycle efficiency tests

This phase adds one repeated MPI driver for the three operations seen by a KMC
application: initial grid construction/classification, cached coordinate
queries, and collective deposition repair. The driver emits one `key=value`
field per line so benchmark automation can archive and compare runs.

## Build

```sh
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON \
    -DGASACCESS_ENABLE_MPI=ON \
    -DGASACCESS_MPI_TEST_RANK_COUNTS="1;2;4;8"
cmake --build build --parallel
```

`GASACCESS_MPI_TEST_RANK_COUNTS` controls all registered MPI rank matrices and
defaults to `1;2;4;8`.

## Timing method

- At least three untimed warmups run by default.
- Initialization and repair repeat until both the requested minimum repetition
  count and cumulative measured duration are reached.
- The reported operation time for each repetition is the maximum elapsed time
  across ranks. `average_seconds` is the arithmetic mean of those maximum-rank
  samples.
- Query timing uses a large allocation-free batch and reports the slowest-rank
  average query time and global throughput.
- Repair fixture construction, initial classification, correctness checks, and
  state restoration occur outside the timed `apply_deposition()` interval.
- Minimum, maximum, and standard deviation accompany the primary average.

The default minimum measured duration is one second. Short registered smoke
tests override it so routine CTest runs remain practical; release measurements
should retain the default or request a longer duration.

## Initialization

Initialization starts with mock SPPARKS domain and atom arrays but no GasAccess
grid. Each repetition constructs the distributed grid, adapts atoms, voxelizes
owned occupancy, performs full distributed classification, and synchronizes
face ghosts before it returns.

```sh
mpiexec -n 8 ./build/gasaccess_mpi_efficiency_driver \
    --operation initialization \
    --scenario million-slab \
    --min-measured-seconds 1 \
    --min-repetitions 10
```

The `million-slab` fixture uses a `128x128x128` grid with spacing `1.0`. Its
lower 64 layers contain 1,048,576 atoms and its upper 64 layers are gas.
Initialization reports elapsed statistics, visited voxels, and distributed
frontier traffic. Every owned voxel must finish classified and sent/received
frontier totals must balance.

## Coordinate queries

The query benchmark initializes once, precomputes owned atom coordinates, and
times only repeated calls to the cached seven-voxel coordinate query.

```sh
mpiexec -n 8 ./build/gasaccess_mpi_efficiency_driver \
    --operation query \
    --scenario million-slab \
    --query-count 1000000 \
    --min-measured-seconds 1
```

`--query-count` is the minimum number of calls per rank. The driver continues
whole coordinate batches until both that count and the requested measured
duration are reached.

## Deposition detection and repair

All ranks call the same collective `apply_deposition()` operation. A rank may
have no changed owned voxel, but remains in change detection, repair
termination, and halo synchronization. Four deterministic cases isolate the
main paths:

| Case | Geometry | Expected work |
|---|---|---|
| `baseline` | Deposit into an already-solid barrier voxel | Detection only; no geometry change or repair |
| `best` | Close an `8x8x8` pocket adjacent to the reservoir | 512 newly closed voxels on one rank |
| `medium` | Close the lower quarter of the global grid | Approximately 25% of voxels close across ranks |
| `worst` | Close the lower three quarters of the global grid | Approximately 75% of voxels close and every rank participates |

Example matrix:

```sh
for ranks in 1 2 4 8; do
    for case_name in baseline best medium worst; do
        mpiexec -n "$ranks" ./build/gasaccess_mpi_efficiency_driver \
            --operation repair \
            --case "$case_name" \
            --nx 128 --ny 128 --nz 128 \
            --min-measured-seconds 1 \
            --min-repetitions 10
    done
done
```

The x dimension must be divisible by the rank count. The `best` fixture also
requires at least 12 x voxels per rank so its pocket and shell remain local.

Repair output includes:

- end-to-end `apply_deposition()` statistics;
- expected, visited, closed, and changed voxel counts;
- count of ranks that visited and closed repair voxels;
- minimum and maximum per-rank visits and closures;
- repair communication rounds;
- sent and received frontier entries; and
- average seconds per visited voxel.

The worst case is an explicit all-rank acceptance test: it fails unless every
rank visits repair voxels, more than half the grid becomes closed, the upper
reservoir-connected control remains accessible, and communication counts
balance.

## Correctness coverage

`gasaccess_distributed_updater_tests` retains differential comparison among the
incremental MPI result, forced full distributed reclassification, and the
serial reference. Its dominant-cavity fixture closes 384 of 512 voxels and
asserts nonzero repair participation on every rank. The registered suite runs
that fixture at 1, 2, 4, and 8 ranks.

Efficiency timings are measurements rather than strict performance thresholds.
Correct state, path, ownership, participation, and communication invariants are
the pass/fail conditions.
