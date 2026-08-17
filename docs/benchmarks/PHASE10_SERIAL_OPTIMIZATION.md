# Phase 10 serial scale and storage optimization

Date: 2026-08-17

Phase 10 measured the serial hot paths, removed the demonstrated production
update bottleneck, and retained exact reference results. The production atomic
structure was not yet available, so these measurements use the deterministic
generators retained from Phase 6.

## Environment

- Build: CMake `Release`, C++17, configured warnings plus `-Werror`
- Compiler: GCC 8.5.0
- CMake: 3.26.5
- OS: Linux 6.6.87.2 under WSL2, x86-64
- CPU: Intel Core i7-13700, 24 logical CPUs
- Execution: one process and one thread; no MPI

Timings below are local measurements, not portable performance guarantees.
The driver now reports per-update median, p95, and maximum latency so future
runs can retain distributions rather than only aggregate time.

## Measured bottleneck and implementation

Before Phase 10, every deposition update copied all voxel states, counted them,
and scanned the full array again to discover which voxels became solid. The
million-voxel workload spent 23.901399 ms on three locally safe one-voxel
updates even though connectivity traversal was already avoided.

Phase 10 makes the production affected-region mode local:

- `GasGrid` maintains four state counts while states change, making the current
  classification summary an O(1) lookup.
- `AtomVoxelizer` can record exact newly solid voxel IDs and their previous
  states during its already-required bounded candidate traversal.
- `DepositionUpdater` reuses this change-record storage and no longer snapshots,
  counts, or searches the full state array in production mode.
- `ConnectivityRepairMode::FullReclassification` intentionally retains the
  full snapshot and state diff as the correctness/debug reference.

The persistent state field remains one byte per voxel. The four 64-bit counters
add only 32 bytes per grid, independent of grid size.

## Million-atom comparison

Command:

```sh
./build/gasaccess_reference_driver \
    --scenario bulk --nx 128 --ny 64 --nz 128 \
    --atom-count 1000000 --query-count 1000000 --update-count 3 \
    --seed 12345
```

| Measurement | Before Phase 10 | After Phase 10 |
|---|---:|---:|
| Atom records | 1,000,000 | 1,000,000 |
| Voxels | 1,048,576 | 1,048,576 |
| Three locally safe updates | 23.901399 ms | 0.004028 ms |
| Safe-update median | not recorded | 1.163 us |
| Safe-update p95/max | not recorded | 1.260 us |
| One million queries | 50.973350 ms | 50.617192 ms |
| Query throughput | 19.62 million/s | 19.76 million/s |
| Peak RSS | 40,247,296 B | 40,321,024 B |

The three-update aggregate improved by approximately 5,934 times in this run.
The small RSS difference is normal process-level measurement variation and does
not indicate per-voxel growth.

Correctness checksums are unchanged:

```text
initial_state_checksum=0x34a6ae79e2f13ec8
final_state_checksum=0x32609a43f1cd6264
```

Voxelization and initial classification remained close to their prior scale:
67.691128 ms and 22.311115 ms respectively in the recorded optimized run.

## Common updates and repair tail

A 1,000-update run on the same million-atom structure separated 872 locally
safe events from 128 affected-region repairs:

```sh
./build/gasaccess_reference_driver \
    --scenario bulk --nx 128 --ny 64 --nz 128 \
    --atom-count 1000000 --query-count 1000000 --update-count 1000 \
    --seed 12345
```

| Path | Samples | Median | p95 | Maximum |
|---|---:|---:|---:|---:|
| Locally safe | 872 | 0.127 us | 1.630 us | 32.676 us |
| Affected-region repair | 128 | 16.772 ms | 18.555 ms | 22.392 ms |

These synthetic bulk repairs collectively visited 35,310,024 voxels. They are
deliberately reported separately: their cost reflects the searched components,
whereas the common proven-safe path no longer scales with total grid volume.

## Deterministic pinch-off

The new `pinch-off` scenario deposits a 640-voxel roof over an open trench. The
first 639 deposits are locally safe; the final deposit closes and repairs the
39,680-voxel cavity.

```sh
./build/gasaccess_reference_driver \
    --scenario pinch-off --nx 64 --ny 32 --nz 64 \
    --query-count 1000000 --update-count 640 --seed 12345
```

| Measurement | Result |
|---|---:|
| Safe-update median | 1.265 us |
| Safe-update p95 | 1.319 us |
| Safe-update maximum | 10.461 us |
| Pinch-off repair latency | 3.184 ms |
| Repair visits/newly closed | 39,680 / 39,680 |
| Total time for all 640 updates | 3.998 ms |

The final checksum is `0x476e87e5f025d7e5`, identical to constructing the
equivalent sealed trench directly.

## Storage decisions

The dense backend remains appropriate for the tested scale:

- one million voxels require exactly 1 MiB of persistent gas states;
- the one-million-atom input occupies most of the roughly 40 MiB peak RSS in
  the reference driver, but the core library receives it as a non-owning view;
- affected-region epoch arrays and frontier storage appear only when repair is
  required and remain proportional to the local dense grid.

Packing the one-byte state into two bits would save only 0.75 MiB at this grid
size while complicating and potentially slowing the hot query. A compile-time
float ABI would not reduce persistent voxel storage and is not justified without
the production KMC atom layout. Therefore Phase 10 adds neither packed state nor
float input, and no tiled/chunked backend phase is currently required.

If the production structure uses far more voxels, or profiling on the target
KMC hardware changes these conclusions, tiled storage and explicit float input
entry points should be reconsidered as separately scoped work.
