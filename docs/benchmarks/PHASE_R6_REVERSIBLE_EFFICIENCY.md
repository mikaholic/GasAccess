# Phase R6 reversible atom-change efficiency and scaling

Date: 2026-09-05

Phase R6 extends the MPI lifecycle benchmark from adsorption-only repair to
adsorption, desorption, and atomic mixed batches. It covers detection-only,
best, medium, and worst cases at 1, 2, 4, and 8 ranks: 48 repeated benchmark
configurations in total.

The complete machine-readable results are archived in
[`PHASE_R6_RESULTS.csv`](PHASE_R6_RESULTS.csv).

## Benchmark contract

The repair driver accepts:

```text
--operation repair
--change-kind adsorption|desorption|mixed
--case baseline|best|medium|worst
```

Every configuration constructs identical incremental and forced-full grids,
applies the same atomic batch to each, and verifies that final owned states,
blocker counts, and changed-coordinate lists are identical. The two timing
series stop independently so a very fast incremental case does not force
thousands of unnecessary full classifications.

Each timing series uses:

- a Release build;
- a fixed `128x128x128` global grid (2,097,152 voxels);
- 3 untimed warmups;
- at least 10 measured repetitions; and
- at least 1 second of cumulative measured time.

Only the complete collective `apply_atom_changes()` call is timed. Fixture
construction, initial classification, state restoration, and correctness
checks are outside the interval. A sample is the maximum elapsed rank time;
the primary result is the arithmetic mean of those samples. The CSV also
records standard deviation and repetition count. The executable additionally
prints minimum, maximum, cumulative time, time per visited voxel, and separate
closing/opening communication metrics.

Environment:

- GCC 8.5.0, `-O3 -DNDEBUG`;
- Open MPI/OpenRTE 4.1.1;
- WSL2 Linux 6.6.87.2;
- Intel Core i7-13700, 12 cores/24 hardware threads; and
- x-only MPI decomposition with identical global dimensions at every rank
  count.

## Fixtures

| Change kind | Baseline | Best | Medium | Worst |
|---|---|---|---|---|
| Adsorption | Add a second blocker to an already-solid plug | Close one local `8x8x8` pocket | Close 25% of the grid | Close 75% of the grid across every rank |
| Desorption | Remove one of two overlapping blockers | Open one local `8x8x8` pocket | Open 25% of the grid | Open more than half the grid across every rank |
| Mixed | Add and remove the same footprint, giving zero net change | Close and open disjoint local `8x8x8` pockets | Close/open disjoint components totaling about 25% | Both passes cross every rank and affect about 75% of the grid |

The mixed fixtures commit both occupancy directions before either connectivity
pass runs. They therefore measure one atomic KMC batch, not sequentially
observable adsorption and desorption events.

## Repeated timing results

Each cell is `incremental ms / forced-full ms (full/incremental)`. A ratio
above 1 means incremental repair is faster; a ratio below 1 means the affected
work is large enough that full classification is faster.

### Adsorption

| Ranks | Baseline | Best | Medium | Worst |
|---:|---:|---:|---:|---:|
| 1 | 0.000172 / 8.041754 (46862.36x) | 0.131968 / 128.440214 (973.27x) | 206.300104 / 123.558158 (0.60x) | 375.840060 / 132.263748 (0.35x) |
| 2 | 0.001063 / 4.231130 (3981.92x) | 0.165791 / 66.016935 (398.19x) | 117.961470 / 63.250126 (0.54x) | 234.371558 / 62.974577 (0.27x) |
| 4 | 0.002218 / 3.327218 (1499.77x) | 0.162489 / 33.167354 (204.12x) | 63.060194 / 31.878836 (0.51x) | 137.080315 / 31.979370 (0.23x) |
| 8 | 0.003539 / 1.346500 (380.50x) | 0.162300 / 20.039976 (123.48x) | 40.495184 / 19.165228 (0.47x) | 94.240356 / 19.008211 (0.20x) |

### Desorption

| Ranks | Baseline | Best | Medium | Worst |
|---:|---:|---:|---:|---:|
| 1 | 0.000159 / 8.050103 (50522.03x) | 0.120808 / 127.986763 (1059.42x) | 89.753551 / 147.636404 (1.64x) | 327.877444 / 211.597059 (0.65x) |
| 2 | 0.001056 / 4.223125 (4000.25x) | 0.146782 / 65.965506 (449.41x) | 59.622236 / 85.489115 (1.43x) | 209.240711 / 134.350868 (0.64x) |
| 4 | 0.002294 / 2.192853 (955.76x) | 0.144656 / 33.886496 (234.26x) | 35.841053 / 48.264402 (1.35x) | 123.751516 / 83.666738 (0.68x) |
| 8 | 0.003620 / 1.321586 (365.06x) | 0.187110 / 20.889174 (111.64x) | 25.851777 / 34.164988 (1.32x) | 87.581191 / 61.100452 (0.70x) |

### Mixed adsorption and desorption

| Ranks | Baseline | Best | Medium | Worst |
|---:|---:|---:|---:|---:|
| 1 | 0.000199 / 8.200990 (41213.18x) | 0.239376 / 131.643603 (549.94x) | 232.601425 / 138.143278 (0.59x) | 452.658069 / 168.110200 (0.37x) |
| 2 | 0.001121 / 4.204522 (3750.48x) | 0.270229 / 65.885180 (243.81x) | 131.420392 / 74.159125 (0.56x) | 271.990222 / 97.887037 (0.36x) |
| 4 | 0.002366 / 2.211771 (934.62x) | 0.262851 / 33.926843 (129.07x) | 69.002566 / 39.986612 (0.58x) | 158.890293 / 57.538870 (0.36x) |
| 8 | 0.003580 / 1.372530 (383.38x) | 0.268492 / 19.789309 (73.71x) | 42.983441 / 26.209802 (0.61x) | 102.482327 / 42.307019 (0.41x) |

The forced-full baseline does not run classification because no geometry
transition occurs. It still snapshots all owned states before detection, as
required by forced-full mode. Baseline ratios therefore compare the normal
detection path with reference-mode setup, not incremental flood fill with a
completed full classification.

## Scaling and work

| Case | Adsorption 1-to-8-rank speedup | Desorption speedup | Mixed speedup |
|---|---:|---:|---:|
| Best, local fixed work | 0.81x | 0.65x | 0.89x |
| Medium, distributed | 5.09x | 3.47x | 5.41x |
| Worst, all-rank | 3.99x | 3.74x | 4.42x |

The local best case intentionally remains on one owner, so adding ranks cannot
divide its 512-voxel component and only adds collective overhead. Medium and
worst cases distribute work over x and show useful strong scaling.

The 8-rank worst-case traversal checks were:

| Change kind | Closing visited / relabeled | Opening visited / relabeled | Participating ranks | Rounds (close/open) | Sent entries (close/open) |
|---|---:|---:|---:|---:|---:|
| Adsorption | 1,644,288 / 1,572,864 | 0 / 0 | 8 / 0 | 6 / 0 | 179,968 / 0 |
| Desorption | 0 / 0 | 1,572,865 / 1,572,865 | 0 / 8 | 0 / 5 | 0 / 172,033 |
| Mixed | 857,856 / 786,432 | 774,145 / 774,145 | 8 / 8 | 6 / 7 | 93,952 / 84,674 |

For every row, globally sent and received frontier-entry totals are equal.
Worst desorption opens 1,572,865 voxels, more than half the grid. Worst mixed
closes 786,432 and opens 774,145 voxels; both passes visit every rank.

Summed peak RSS ranges from about 116 MB for the paired one-rank baseline to
493 MB for the paired eight-rank adsorption worst case. These values include
two simultaneously resident benchmark grids—incremental and forced-full—and
must not be interpreted as the memory requirement of one production grid.

## Interpretation

- Detection-only updates remain microsecond-scale or faster: 0.159-0.199 us
  at one rank and 3.539-3.620 us at eight ranks.
- A local 512-voxel repair completes in 0.121-0.268 ms and is 73.7x to
  1059x faster than forced full classification.
- Medium desorption is 1.32x to 1.64x faster than full classification while
  opening 524,289 voxels.
- Large distributed repairs scale, but incremental traversal is not always
  cheaper than a full scan. The medium closing/mixed cases and all worst cases
  expose the expected crossover. A future adaptive policy could select full
  reclassification when estimated affected work is a large fraction of the
  grid.

The relevant KMC regime is dominated by local events, where the benchmark
shows very low update cost and avoids a global scan. The all-rank fixtures are
safety and scalability tests, not an assertion that incremental repair should
win when most of the grid changes.

## Historical Phase R1 adsorption comparison

Phase R1 measured the adsorption-specific path immediately before and after
adding reversible blocker counts. Those same-grid results are retained here as
the adsorption regression reference:

| MPI ranks | Detection before/after (us) | Best before/after (ms) | Medium before/after (ms) | Worst before/after (ms) |
|---:|---:|---:|---:|---:|
| 1 | 0.079 / 0.085 | 4.218 / 4.297 | 238.43 / 233.84 | 354.95 / 366.92 |
| 2 | 0.423 / 0.420 | 2.354 / 2.327 | 122.67 / 123.66 | 230.72 / 238.79 |
| 4 | 0.945 / 0.940 | 1.292 / 1.284 | 69.74 / 70.83 | 140.59 / 142.12 |
| 8 | 1.345 / 1.383 | 0.836 / 0.874 | 45.30 / 45.50 | 101.21 / 102.20 |

The Phase R1 before/after result showed no material adsorption regression.
Current R6 timings use the unified R5 `apply_atom_changes()` path and a paired
incremental/full harness, so they are reported separately rather than treated
as a controlled continuation of that historical before/after experiment.

## Correctness and completion

All 48 release configurations completed with `acceptance=pass`. Each checked:

- exact incremental/full state and blocker-count equality;
- exact transition, relabel, final-state, and changed-coordinate counts;
- sorted, unique, identical incremental/full changed-coordinate lists;
- reservoir, closing-cavity, and opening-cavity accessibility probes;
- balanced closing and opening frontier traffic;
- one-rank locality for best cases; and
- all-rank participation for the relevant worst-case pass or passes.

CMake also registers short one-repetition smoke variants for all 48
configurations. Together with initialization and query efficiency tests, all
56 efficiency-labeled tests pass. A clean Release rebuild completed without
compiler warnings, and the complete project regression passed 102 of 102
tests, including every 1/2/4/8-rank MPI test.
