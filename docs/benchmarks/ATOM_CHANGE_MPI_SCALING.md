# Repeated atom-change MPI scaling benchmark

Date: 2026-09-05

## Contract

This benchmark measures the incremental collective `apply_atom_changes()`
operation for adsorption-only, desorption-only, and combined atom-change
batches. Every change kind is tested with best, medium, and worst repair
fixtures at 1, 2, 4, and 8 MPI ranks.

- Release build (`-O3 -DNDEBUG`)
- fixed `128x128x128` global grid (2,097,152 voxels)
- x-only strong-scaling decomposition
- 3 untimed warmups
- at least 10 measured repetitions and at least 1 second of cumulative
  incremental timing per configuration
- one sample is the maximum elapsed time across ranks
- reported time is the arithmetic mean of those collective samples
- speedup is the matching one-rank mean divided by the multi-rank mean
- setup, restoration, and correctness validation are outside the timed interval

Every configuration is compared with a forced full reclassification. Final
owned states, blocker counts, and changed-coordinate lists must match exactly.

## Fixtures

- **Best:** one local `8x8x8` component per repair direction; repair remains
  on one MPI owner.
- **Medium:** distributed components affect approximately 25% of the grid.
- **Worst, adsorption:** sealing one inlet closes 1,572,864 voxels across every
  rank.
- **Worst, desorption:** opening one inlet makes 1,572,865 voxels accessible
  across every rank.
- **Worst, combined:** one large component initially uses an old inlet. The
  atomic batch adsorbs at that inlet and desorbs at a distant new inlet. The
  closing pass relabels 1,572,864 voxels and the opening pass relabels
  1,572,865 voxels, with both traversals reaching every rank. The component
  ends accessible, so only the two inlet voxels differ between initial and
  final states.

## Results

### Adsorption

| MPI ranks | Best time / speedup | Medium time / speedup | Worst time / speedup |
|---:|---:|---:|---:|
| 1 | 0.135 ms / 1.00x | 207.147 ms / 1.00x | 377.844 ms / 1.00x |
| 2 | 0.160 ms / 0.84x | 117.593 ms / 1.76x | 232.256 ms / 1.63x |
| 4 | 0.162 ms / 0.83x | 61.882 ms / 3.35x | 135.348 ms / 2.79x |
| 8 | 0.152 ms / 0.89x | 39.350 ms / 5.26x | 93.914 ms / 4.02x |

### Desorption

| MPI ranks | Best time / speedup | Medium time / speedup | Worst time / speedup |
|---:|---:|---:|---:|
| 1 | 0.116 ms / 1.00x | 86.629 ms / 1.00x | 318.224 ms / 1.00x |
| 2 | 0.144 ms / 0.81x | 58.612 ms / 1.48x | 207.070 ms / 1.54x |
| 4 | 0.139 ms / 0.84x | 35.446 ms / 2.44x | 125.045 ms / 2.54x |
| 8 | 0.133 ms / 0.88x | 25.308 ms / 3.42x | 88.305 ms / 3.60x |

### Combined adsorption + desorption

| MPI ranks | Best time / speedup | Medium time / speedup | Worst time / speedup |
|---:|---:|---:|---:|
| 1 | 0.222 ms / 1.00x | 226.782 ms / 1.00x | 658.133 ms / 1.00x |
| 2 | 0.251 ms / 0.88x | 125.341 ms / 1.81x | 394.586 ms / 1.67x |
| 4 | 0.266 ms / 0.83x | 67.947 ms / 3.34x | 285.914 ms / 2.30x |
| 8 | 0.250 ms / 0.89x | 42.649 ms / 5.32x | 230.871 ms / 2.85x |

## Interpretation

The local best fixtures intentionally do not distribute work, so additional
ranks add collective overhead and produce speedups below 1. Medium repairs
show useful strong scaling: at 8 ranks, adsorption, desorption, and combined
updates reach 5.26x, 3.42x, and 5.32x speedup, respectively. Worst adsorption
and desorption reach 4.02x and 3.60x. The adversarial combined worst case
processes the same large component twice and still reaches 2.85x at 8 ranks.

Exact measurements and traversal counters are in
[`ATOM_CHANGE_MPI_SCALING_RESULTS.csv`](ATOM_CHANGE_MPI_SCALING_RESULTS.csv).
