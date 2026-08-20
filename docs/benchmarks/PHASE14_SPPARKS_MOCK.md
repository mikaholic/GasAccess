# Phase 14 SPPARKS mock acceptance

Date: 2026-08-19

## Environment

- Linux 6.6.87.2 under WSL2
- Intel Core i7-13700, 12 physical cores / 24 logical CPUs
- GCC 8.5.0, release build (`-O3`)
- OpenMPI/OpenRTE 4.1.1
- one local machine; these numbers are not multi-node network measurements

## Workload

The `million-slab` fixture has a 128 x 128 x 128 gas grid. Its lower 64 voxel
layers contain one atom at each voxel center, giving exactly 1,048,576 atoms;
the upper 64 layers are free gas. The mock SPPARKS domain is periodic in all
axes. GasAccess is periodic in x/y, non-periodic in z, and uses only the top-z
face as a reservoir.

Each run performs owned-atom generation, neighbor-directed atom ghost
exchange, SPPARKS-to-GasAccess adaptation, distributed voxelization, full
distributed flood-fill, and one cached query for every owned atom. Times are
the maximum across ranks. Peak RSS is summed across ranks.

Command pattern:

```sh
mpiexec -n RANKS ./build/gasaccess_spparks_mock_driver \
    --scenario million-slab
```

## Results

| Ranks | Process grid | Total (s) | Voxelize (s) | Classify (s) | Query (s) | Queries/s | RSS sum (bytes) |
|---:|:---:|---:|---:|---:|---:|---:|---:|
| 1 | 1x1x1 | 0.321573 | 0.075007 | 0.078059 | 0.053428 | 19,625,993 | 197,046,272 |
| 2 | 2x1x1 | 0.167258 | 0.035720 | 0.036177 | 0.027468 | 38,174,025 | 233,545,728 |
| 4 | 2x2x1 | 0.093688 | 0.019398 | 0.019852 | 0.014442 | 72,603,939 | 279,031,808 |

Every rank count produced the same physical result:

- 1,048,576 solid voxels;
- 1,048,576 outside-accessible gas voxels;
- zero closed-void voxels;
- exactly 16,384 accessible atoms, corresponding to the top solid layer;
- `acceptance=pass`.

The observed end-to-end speedups relative to one rank were approximately
1.92x on two ranks and 3.43x on four ranks. These are synthetic single-machine
measurements, not a performance promise for the production KMC structure.

## Correctness and regression tests

The final release test command was:

```sh
ctest --test-dir build --output-on-failure
```

All 31 registered tests passed. The new tests include:

- open/sealed static integration on one, two, and four ranks;
- x-, y-, z-, and 2x2 process decompositions where applicable;
- distributed state and site-query comparison with the serial reference;
- owned and face-ghost state comparison;
- independent SPPARKS and GasAccess boundary settings;
- atom ghost-cutoff rejection;
- million-atom acceptance on one, two, and four ranks;
- direct compile compatibility with the available SPPARKS `Domain` and `App`
  headers.

A one-rank sealed-trench Valgrind run reported no invalid GasAccess memory
access. Full leak reporting showed only small process-lifetime allocations
whose stacks originate in this OpenMPI/OpenRTE component loader and `MPI_Init`,
consistent with the earlier MPI phase checks.
