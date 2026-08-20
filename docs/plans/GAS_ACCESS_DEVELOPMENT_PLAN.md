# GasAccess Development Plan

Status: in development — Phase 13 complete
Last updated: 2026-08-19

This is the working plan for developing GasAccess as an independent C++/C
library and later integrating it with an MPI-parallel kinetic Monte Carlo
simulator. It should be updated as interfaces, measurements, and physical
requirements become known.

## 1. Objective

GasAccess will determine whether a reaction site is connected to an external
gas reservoir through space that can accommodate a configured finite-radius
precursor. Connectivity is evaluated through an auxiliary voxel field, not
through solid-atom connectivity.

The runtime architecture must be:

```text
rare geometry change
        |
        v
update cached gas connectivity
        |
        v
many reaction queries
        |
        v
fixed-cost cached lookup
```

The library has two eventual acceptance goals:

1. Process a supplied structure containing millions of atoms and report
   initialization time, update time, query throughput, and memory consumption.
2. Integrate with the existing KMC domain decomposition and propensity-update
   workflow without adding MPI communication to ordinary accessibility queries.

## 2. Project-wide requirements

### 2.1 Naming convention

All production code and tests will use these conventions:

| Construct | Convention | Example |
|---|---|---|
| Classes and public types | `PascalCase` | `GasGrid`, `GridSpec` |
| Functions and methods | `snake_case` | `apply_deposition()` |
| Variables | `snake_case` | `voxel_count` |
| Private data members | trailing-underscore `snake_case_` | `states_` |
| C API functions | prefixed `snake_case` | `ga_is_accessible()` |

Enum values will follow the established scoped style, for example
`GasState::Solid`, unless a later project-wide convention replaces it.

### 2.2 Boundary conditions

Periodicity is independently configurable for x, y, and z. All eight periodic
axis combinations must be supported.

The GasAccess boundary configuration is independent of the host KMC boundary
configuration. In particular, a production KMC domain may remain periodic in
all three axes while its GasAccess view uses periodic x/y, non-periodic z, a
top-z reservoir face, and a closed/non-source bottom-z face. The adapter must
use the caller-supplied GasAccess settings rather than copying the SPPARKS
periodicity blindly.

Periodic behavior applies consistently to:

- voxel neighbor traversal;
- atoms whose exclusion spheres cross a periodic seam;
- reaction-site lookup near a seam;
- initial flood-fill;
- local topology tests;
- incremental connectivity repair;
- MPI ownership and ghost exchange.

For each axis, a periodic setting identifies the low and high faces. A periodic
face is not an external reservoir boundary. Each non-periodic face may be
configured as either a gas-reservoir face or a closed/non-source face.

An explicit set or mask of reservoir seed voxels will also be supported by the
connectivity interface. This permits specialized configurations, including
fully periodic domains. A domain with no reservoir faces and no explicit seeds
has no external gas source; all empty voxels are therefore closed with respect
to external gas.

### 2.3 Initial physical model

- Axis-aligned Cartesian grid with uniform spacing along each individual axis.
  The x, y, and z spacings may differ, so voxels may be rectangular cuboids
  rather than cubes.
- Six-face connectivity between gas voxels.
- No gas leakage through edge-only or corner-only contact.
- Connectivity is position-only geometric connectivity. Molecular orientation,
  rotational motion, and orientation-dependent transport are not represented.
- The KMC simulator remains responsible for orientation-dependent reaction
  physics at a reaction site.
- One isotropic effective precursor radius is a required physical parameter for
  each production grid instance. A zero value is permitted only for point-probe
  reference tests or models that explicitly require it.
- The initial finite-size model is a hard-sphere steric-exclusion approximation.
  It represents whether the precursor can geometrically fit through free space;
  energetic barriers and orientation-dependent reaction rates remain the KMC
  simulator's responsibility.
- A voxel is blocked when its center is inside the excluded radius
  `R_atom + R_precursor` of a solid atom.
- Expanding every atom by `R_precursor` must close a throat whose physical
  clearance is too small for the precursor, even if the underlying atom surfaces
  do not touch. Consequently, accessibility is specific to the configured
  precursor radius.
- Monotonic deposition is implemented first: gas voxels may become solid, but
  solid voxels do not become gas.
- Gas-to-gas connectivity and reaction-site-to-gas adjacency are separate
  policies. Gas connectivity remains six-neighbor even if a different fixed
  site-contact stencil is later justified.

The graph algorithms do not require cubic voxels: six-neighbor flood-fill,
cached site queries, the index-space `3x3x3` topology filter, and affected-region
repair depend only on voxel indices and face adjacency. Geometry operations do
require per-axis handling. Coordinate lookup, voxel centers, periodic lengths,
and atom candidate bounds must use the corresponding x/y/z spacing. Spherical
steric exclusion remains a sphere in physical coordinates: the final test is
the Euclidean distance to the voxel center, not a scaled index distance that
would incorrectly turn the sphere into an ellipsoid.

### 2.4 Performance requirements

- Accessibility queries must not flood-fill, allocate memory, or communicate.
- Query work must be bounded by a small, fixed site-neighbor stencil.
- Normal queries use only owned or ghost state in MPI builds.
- Full flood-fill is allowed for initialization, testing, debugging, and
  exceptional fallback, but not as the normal optimized update path.
- Production MPI code must not replicate the global gas graph or routinely
  all-gather all voxels.

### 2.5 Voxel and state representation

Voxels will be implicit grid cells rather than individually allocated objects.
Their coordinates are derived from a linear identifier, and their states are
stored contiguously:

```cpp
using VoxelId = std::uint64_t;

struct VoxelCoord {
    std::int64_t x;
    std::int64_t y;
    std::int64_t z;
};

enum class GasState : std::uint8_t {
    Unclassified = 0,
    Solid = 1,
    OutsideAccessible = 2,
    ClosedVoid = 3
};
```

The initial dense backend will use `std::vector<GasState> states_`, giving one
persistent state byte per owned or ghost voxel. It will not use
`std::vector<bool>`. `Unclassified` is a construction/debug state and must not
remain after a successful classification.

Public/global voxel identifiers are 64-bit. Integer voxel coordinates are
signed so neighbor offsets can be calculated safely before boundary rejection
or periodic wrapping. Local container indexing may use `std::size_t` after
checked conversion. Linearization, dimension multiplication, and conversions
must be checked for integer overflow.

The C API will represent stored state with `uint8_t` and named `GA_GAS_*`
constants rather than relying on the compiler-dependent storage size of a C
enum. MPI state exchange will likewise send the fixed-width byte representation.

`GridSpec` will use a three-component spacing value in Phase 11. Each component
is a finite positive `double`; the C representation will likewise use three
components. An isotropic grid is represented by equal components. This is a
geometry/API generalization only and does not change voxel identifiers, state
storage, or the six-neighbor graph.

Bit-packing is deliberately deferred. If Phase 10 demonstrates that state
memory is limiting, two measured alternatives are a two-bit packed state or
separate `occupied` and `outside_connected` bitplanes. Any alternative must be
benchmarked against the one-byte representation because bit extraction may
slow the query hot path.

### 2.6 Numerical precision policy

The first implementation will use `double` for atom coordinates, atom radii,
grid origin, voxel spacing, periodic-distance calculations, and atom-to-voxel
geometry tests. Voxel topology and accessibility state remain integer data and
are unaffected by floating-point precision.

There will be no global float-versus-double compilation option initially. Such
an option would change public layouts or ABI if exposed carelessly, double the
required numerical test matrix, and provide little benefit to cached reaction
queries. Atom input does not need to remain resident after voxelization, so
double-precision input also need not create a permanent per-voxel memory cost.

If Phase 10 measurements show a material input-memory or voxelization-bandwidth
benefit, explicit float and double input entry points may be added while keeping
geometry calculations in double precision. This is preferred over a macro that
silently changes the installed library ABI. Supporting float input would not
recover precision absent from the source data, but conversion to double would
keep subsequent calculations consistent.

## 3. Proposed library boundary

The implementation will use C++17 internally and expose both a C++ API and a
thin C API with opaque handles. Its initial geometry scalar type will be
`double`. Exact signatures will be finalized in Phase 0, but the responsibilities
are expected to resemble:

```cpp
class GasAccess {
public:
    void build_from_atoms(AtomView atoms);
    void classify_all();
    UpdateResult apply_deposition(AtomView new_atoms);
};

class GasAccessibilityQuery {
public:
    bool is_site_accessible(const Point3& atom_position) const;
};
```

`UpdateResult` exposes changed voxel identifiers for diagnostics and optional
future optimizations. They are not part of the required KMC query path. The
permanent reference API also provides an explicit full recomputation for
validation.

The core API accepts in-memory array views. Atomic-structure file readers are
adapters, preventing a particular simulation file format from becoming a core
library dependency.

## 4. Development phases

Code estimates are guardrails for new production code and exclude tests,
generated code, and documentation. If a phase grows substantially beyond its
range, it should be split before continuing.

### Milestone A: serial correctness baseline

#### Phase 0: contracts and build scaffold

Scope:

- Establish the CMake project, public/private header layout, tests, and examples.
- Define `GridSpec`, atom input, site input, boundary/source configuration,
  error handling, index types, and units.
- Establish the fixed-width voxel identifier/state types and double-precision
  geometry contract.
- Record the naming convention in contributor-facing configuration/docs.
- Establish a C++ API and opaque-handle C ABI skeleton.

Tests:

- Build and link one C++ client and one pure-C client.
- Verify C and C++ state values and fixed-width representations agree.
- Verify `GasState` occupies one byte and `VoxelId` occupies eight bytes on every
  supported build.
- Reject zero spacing, zero-sized axes, negative or non-finite radii,
  inconsistent periodic faces, invalid source faces, integer-overflowing
  dimensions, and out-of-domain inputs where wrapping is not permitted.

Exit gate:

- Public types and configuration behavior are reviewed before algorithm code
  depends on them.

Expected production code: 100-250 lines.

#### Phase 1: voxel grid geometry and topology

Scope:

- Dense voxel storage and stable local voxel identifiers.
- Coordinate/index conversion and overflow-safe size calculations.
- Checked 64-bit linear identifiers with safe conversion to local indices.
- Six-neighbor iteration with independent x/y/z periodic wrapping.
- Reservoir-face and explicit-source identification.

Tests:

- Coordinate/index round trips.
- Corners, edges, degenerate small dimensions, and out-of-range coordinates.
- Exact neighbor sets for all eight periodic-axis combinations.
- No duplicated neighbors when periodic dimensions are very small.
- Correct identification of non-periodic source faces and explicit seeds.

Exit gate:

- Grid topology is correct without atoms or connectivity classification.

Expected production code: 250-450 lines.

#### Phase 2: static atom voxelization

Scope:

- Convert atom positions and radii into blocked voxels.
- Apply the single precursor radius as spherical steric exclusion.
- Use atom-centered voxel bounding boxes rather than comparing every atom with
  every voxel.
- Wrap exclusion geometry across each enabled periodic seam.

Tests:

- Analytically known blocked cells for a single atom.
- Overlapping atoms and heterogeneous atom radii.
- Pores just below, at, and above the effective precursor diameter.
- Atom surfaces that remain separated while their precursor-exclusion regions
  overlap and seal the passage.
- Atoms touching non-periodic boundaries.
- Atoms crossing x, y, z, and combined periodic seams.
- Atom centers exactly on and immediately around exclusion-distance thresholds.
- Large coordinate origins combined with small voxel spacing.
- Repeatability independent of atom input order.

Exit gate:

- Occupancy matches a slow brute-force test oracle on small randomized grids.

Expected production code: 250-500 lines.

#### Phase 3: initial exterior classification

Scope:

- Seed all empty reservoir voxels.
- Run serial breadth-first flood-fill through six-connected empty voxels.
- Classify voxels as solid, outside-accessible, or closed void.
- Retain this implementation permanently as the reference classifier.

Tests:

- Empty domain, completely solid domain, and domain without a source.
- Open trench and trench with a sealed roof.
- Enclosed cavity and cavity with a one-voxel face-connected channel.
- Edge-only and corner-only contacts remain disconnected.
- Periodic paths crossing each enabled seam.

Exit gate:

- Every outside voxel has a path to a source and every closed voxel lacks one.

Expected production code: 180-350 lines.

#### Phase 4: cached reaction-site queries and C interface

Scope:

- Query gas state by voxel identifier and coordinate.
- Query site accessibility through a bounded site-contact stencil.
- Add the functional C wrapper for construction, classification, and queries.
- Keep site-contact policy independent from six-neighbor gas connectivity.

Tests:

- Surface sites next to outside gas, closed gas, and solid-only neighborhoods.
- Sites adjacent to periodic seams.
- Equivalent results through C and C++ APIs.
- Query instrumentation confirms no graph traversal, allocation, or MPI call.

Exit gate:

- Ordinary reaction checks are fixed-cost cached lookups.

Expected production code: 200-400 lines.

#### Phase 5: deposition with full recomputation

Scope:

- Accept one atom or a small deposition batch.
- Update only potentially affected occupancy voxels.
- Re-run the full reference classification after geometry changes.
- Return the voxels whose accessibility state changed.
- Enforce or clearly report the initial deposition-only update contract.

Tests:

- Staged sidewall growth followed by trench pinch-off.
- Deposition that changes no occupancy.
- Overlapping deposition and multi-voxel deposition.
- Closure across a periodic seam.
- Removal of an empty reservoir seed by deposition.

Exit gate:

- The serial baseline correctly detects sealed voids after arbitrary supported
  deposition sequences.

Expected production code: 250-450 lines.

#### Phase 6: standalone reference driver and baseline measurements

Scope:

- Provide a small standalone executable using the in-memory library API.
- Generate deterministic synthetic open-trench, closed-trench, and bulk cases.
- Add a file-format adapter after the target structure format is supplied.
- Record voxelization, classification, query, and full-update timings.
- Record atom count, voxel count, grid spacing, and peak memory with every run.
- Report persistent bytes per voxel separately from transient traversal memory.

Tests:

- Deterministic state checksums and summary counts.
- End-to-end reproduction of all topology fixtures.
- Successful execution on progressively larger generated structures.

Exit gate:

- Correctness and baseline timing are reviewed before incremental complexity is
  introduced.

Expected production/tool code: 250-500 lines, excluding a complex format reader.

### Milestone B: incremental serial updates

#### Phase 7: conservative local topology filter

Scope:

- Examine accessible neighbors of every newly blocked voxel.
- Prove harmless removal using a fixed local neighborhood, initially 3x3x3.
- Treat inconclusive cases as possible pinch-offs.
- Handle periodic wrapping and source voxels conservatively.
- Continue using full recomputation for suspected pinch-offs.

Tests:

- Known safe deletions, bridges, narrow necks, and periodic-seam bridges.
- Exhaustive small-grid patterns around the removed voxel where practical.
- Differential test: whenever the filter declares an event safe, skipping the
  full recomputation must produce exactly the reference result.

Exit gate:

- The filter has no false-safe results in exhaustive and randomized testing.

Expected production code: 150-300 lines.

#### Phase 8: serial affected-region repair

Scope:

- Start component searches from surviving accessible neighbors of removed gas.
- Detect which resulting regions retain a path to a source.
- Relabel fully explored source-disconnected regions as closed.
- Preserve a configurable full-recompute fallback/debug path.
- Use reusable visitation storage or epochs to avoid per-update large allocations.

Tests:

- Multiple pinch-off shapes and nested cavities.
- Pinch-off involving source-adjacent cells.
- Connectivity through periodic seams.
- Random deposition sequences compared voxel-for-voxel with full flood-fill
  after every event.

Exit gate:

- Incremental state is identical to the reference state for every tested event.

Expected production code: 400-650 lines.

#### Phase 9: single-site KMC query contract

Scope:

- Make the existing call below the only required GasAccess operation while the
  KMC iterates over atoms or reaction sites:

  ```cpp
  bool accessible = query.is_site_accessible(atom_position);
  ```

- Define `atom_position` as the atom's current position supplied by KMC,
  including any position change produced by MD relaxation.
- Document that the default contact stencil is the voxel containing the atom
  plus its six face-neighbor voxels. The result is true when at least one of
  those voxels is `OutsideAccessible`.
- Confirm that a query object can be constructed once and reused: it references
  the current grid state and therefore observes connectivity updates without
  rebuilding or registering atoms.
- Keep reaction creation, rate calculation, atom storage, MD synchronization,
  and KMC iteration outside GasAccess.
- Add no atom identifiers, atom registry, callbacks, affected-atom lists,
  batch-query requirement, or changed-voxel handoff requirement.
- Retain existing auxiliary APIs, including changed voxel IDs and custom voxel
  stencils, but do not require the KMC integration to use them.

Tests:

- A minimal mock KMC loop calls only
  `query.is_site_accessible(atom_position)` for each atom.
- Atoms adjacent to exterior-connected gas return true, while atoms adjacent
  only to solid or closed-void gas return false.
- The same query object returns the new result after a deposition update creates
  a pinch-off or otherwise changes cached connectivity.
- Queries use the atom's latest supplied position after a simulated MD move.
- Periodic seams and non-periodic out-of-domain positions obey the documented
  boundary rules.
- The hot query remains allocation-free and has a fixed upper bound on voxel
  state inspections.

Exit gate:

- The mock KMC integration needs only
  `query.is_site_accessible(atom_position)` to decide whether a gas-dependent
  reaction may be considered.
- The answer always reflects the latest completed GasAccess connectivity update
  and no KMC-owned atom or reaction data is stored by GasAccess.

Expected production/documentation code: 25-100 lines because the core query
already exists; most work is integration testing and contract documentation.

#### Phase 10: serial scale and storage optimization

Scope:

- Benchmark a structure containing millions of atoms.
- Profile voxelization, initial classification, harmless deposition, real
  pinch-off repair, and query throughput separately.
- Measure common-update median and tail latency, fallback frequency, and peak RSS.
- Improve data layout and traversal only where measurements justify it.
- Evaluate float input or packed state only if profiling identifies the
  corresponding input-bandwidth or state-memory bottleneck.
- Decide whether the dense backend is sufficient. If not, introduce a separate
  planned phase for tiled/chunked storage rather than silently expanding scope.

Tests:

- Results and checksums remain identical before and after each optimization.
- If float input is added, its conversion and documented accuracy behavior pass
  the same boundary and periodic-seam cases as double input.
- Smaller benchmark cases remain usable in continuous testing.
- Large performance runs produce reports rather than fragile machine-specific
  pass/fail assertions.

Exit gate:

- Serial turnaround time and memory are documented on the target hardware and
  accepted, or a measured storage/performance issue receives its own phase.

Expected production code: 200-500 lines unless a new storage backend is approved.

### Milestone C: MPI and KMC integration

Milestone C reuses the SPPARKS-style off-lattice decomposition; GasAccess will
not implement an independent domain decomposition. The host supplies a small
decomposition descriptor rather than exposing SPPARKS private bin internals to
the core library. The reviewed SPPARKS data provide global/subdomain bounds,
`procgrid`, `myloc`, per-axis periodicity, and the existing `world` MPI
communicator.

The production off-lattice application currently uses periodic boundaries on
all three axes, while the generic GasAccess MPI implementation and its tests
must also support non-periodic axes.

#### Phase 11: anisotropic grid, decomposition adapter, and gas ghost exchange

Scope:

- Replace the current scalar grid spacing in the C++ and C APIs with explicit
  x/y/z spacing and update coordinate lookup, voxel centers, periodic lengths,
  atom candidate ranges, validation, drivers, and tests.
- Preserve spherical excluded-volume voxelization in physical coordinates when
  the per-axis spacings differ. Do not scale the radius independently by axis
  in the final distance test.
- Add an alignment helper that takes a requested nominal resolution and, for
  each axis, chooses a global voxel dimension divisible by that axis's process
  count, then reports the exact spacing as `box_length / dimension`. Do not
  silently modify a directly supplied `GridSpec`.
- Add an optional MPI decomposition descriptor containing the existing `world`
  communicator, global/local box bounds, `procgrid`, and `myloc`. Gas periodic
  axes and reservoir faces are supplied separately from the KMC boundary
  settings.
- Require exact voxel/decomposition alignment for the initial implementation:
  grid origin equals the global box lower bound, the dimension times spacing
  equals the corresponding box length on each axis, and every global voxel
  dimension is divisible by the process-grid dimension on that axis. Equal
  x/y/z spacing is not required.
- Derive owned global voxel-index ranges by integer arithmetic so no voxel is
  split between ranks and every SPPARKS subdomain face is a voxel face.
- Add one layer of face-connected gas ghost voxels, sufficient for the
  six-neighbor connectivity model and the default position query.
- Derive the six face-neighbor ranks from `procgrid`, `myloc`, and periodicity.
  Use local copies for periodic self-neighbors and no neighbor at a
  non-periodic global edge.
- Exchange the fixed-width one-byte gas states over `world` separately from
  `CommOffLattice`, whose messages contain atom/site arrays rather than gas
  voxels.
- Consume already synchronized owned and ghost atoms when building owned gas
  occupancy. Ignore periodic image atoms across a global face that GasAccess
  treats as non-periodic. Validate again inside GasAccess that the supplied
  host atom ghost distance (the minimum relevant off-lattice bin size) is at
  least `max(R_atom + R_precursor)`, even though KMC is expected to guarantee
  this already.
- Ensure a completed gas halo exchange precedes the KMC query loop; ordinary
  `query.is_site_accessible(atom_position)` calls perform no MPI operation.
- Do not implement distributed flood-fill yet.

Tests:

- Per-axis coordinate/index round trips, voxel centers, periodic wrapping, and
  invalid spacing values on non-cubic grids.
- Non-cubic atom voxelization against a brute-force physical-distance oracle,
  including periodic seams and different spacing orderings.
- Serial classification, local topology filtering, affected-region repair, and
  position queries on non-cubic grids match their reference results, confirming
  that topology behavior is independent of voxel aspect ratio.
- One-rank behavior matches the serial library.
- The alignment helper returns exactly covered boxes, near-target per-axis
  spacings, and divisible dimensions. Direct validation accepts compatible
  grids and rejects misaligned origin, extent, spacing, and process-grid
  divisibility.
- Owned index ranges cover the global grid exactly once on one, two, and four
  ranks.
- Owned/ghost states match on x-, y-, and z-split decompositions.
- Tests cover periodic self-copies, cross-rank periodic seams, non-periodic
  global edges, and queries whose face stencil crosses a rank boundary.
- Ghost-distance validation accepts adequate SPPARKS atom halos and rejects an
  excluded radius larger than the supplied distance.
- A fully periodic KMC descriptor can be paired with a GasAccess configuration
  having periodic x/y, non-periodic z, only the top-z face as a reservoir, and
  no gas-state wrap or periodic atom image across z.

Exit gate:

- Non-cubic geometry is correct, every rank has correct owned and face-ghost
  gas state after exchange, and a normal site query requires only local/ghost
  reads.

Expected production code: 450-750 lines, depending on the supplied adapter API.

#### Phase 12: distributed initial flood-fill

Scope:

- Maintain a local frontier per rank.
- Send frontier entries only across face boundaries to the owning neighbor
  rank using `world`.
- Use a small collective only to detect global termination.
- Avoid global gas-graph replication and routine all-gather.
- Seed non-periodic reservoir faces only on ranks owning those global faces.
  A KMC-periodic axis that the caller configures as non-periodic for GasAccess
  follows this face-source rule; a GasAccess domain that remains fully periodic
  uses configured explicit reservoir voxels.
- Synchronize final face halos before returning control to the KMC query loop.

Tests:

- Compare distributed output with the serial reference for 1, 2, and 4 or more
  ranks using several decompositions.
- Include cavities and access paths crossing rank and periodic boundaries.
- Verify termination when some ranks own no active frontier.

Exit gate:

- Distributed initial classification is decomposition-independent and matches
  the serial voxel state exactly.

Expected production code: 350-600 lines.

#### Phase 13: distributed incremental repair

Scope:

- Keep the local topology filter conservative at rank boundaries. If its
  `3x3x3` proof would require unavailable edge/corner halo data, escalate to
  distributed repair rather than declaring the deletion safe.
- Exchange connectivity-repair frontiers only when required.
- Detect distributed termination.
- Relabel affected components, synchronize face ghosts, and retain changed
  owned voxel IDs as an auxiliary diagnostic result.
- Retain distributed full recomputation as a debug/fallback mode.
- Do not introduce atom/site registration, affected-atom lists, or selective
  propensity invalidation.

Tests:

- Pinch-offs occurring within a rank and exactly on rank boundaries.
- Cavities spanning several ranks.
- Random distributed deposition sequences compared with the serial reference.
- Results are identical under multiple decompositions of the same domain.

Exit gate:

- Incremental MPI results match full distributed and serial recomputation.

Expected production code: 450-750 lines.

#### Phase 14: KMC integration and production acceptance

Scope:

- Adapt the real atom, deposition, and SPPARKS `Domain` interfaces into the
  neutral GasAccess descriptors.
- Supply the effective GasAccess periodic axes and reservoir faces independently
  from the SPPARKS boundary flags; the expected production configuration is
  periodic x/y, non-periodic z, with the top-z face connected to the reservoir.
- Build/initialize the grid from the KMC structure.
- Schedule work as: synchronize SPPARKS atoms, update gas occupancy and
  connectivity, exchange gas-state halos, then enter the KMC site loop.
- Use only `query.is_site_accessible(atom_position)` while KMC iterates sites;
  KMC retains responsibility for its normal full rate/event recalculation.
- Use the incremental deposition updater only for monotonic gas-to-solid
  changes. If MD relaxation can unblock or move occupancy between voxels,
  rebuild occupancy from the synchronized current atoms and run distributed
  reclassification as the correctness baseline.
- Run the supplied million-atom structure and representative event sequence.
- Measure strong scaling, communication volume, update latency, query throughput,
  and end-to-end KMC overhead.

Tests:

- A mock KMC integration test precedes changes to the real simulator.
- Integrated results match the standalone library on the same geometry/events.
- Normal site queries perform no collectives and require only local/ghost state.
- Atom ghost cutoff validation is exercised against the maximum excluded
  radius used by the integration.
- MD-relaxed geometries match a clean rebuild from the same current atoms.
- Debug full recomputation periodically checks incremental state during long runs.

Exit gate:

- Physical results are accepted and measured performance meets targets agreed
  from Phase 10 and the KMC baseline.

Code size depends on the KMC interfaces and will be estimated after those
interfaces are reviewed.

## 5. Testing strategy retained throughout development

### Deterministic unit tests

Small grids will make the exact expected state visible. Tests will cover open
paths, sealed cavities, diagonal contacts, source rules, periodic seams, and
degenerate dimensions.

### Reference differential tests

The full flood-fill remains the correctness oracle. Incremental state will be
compared against it after each event in randomized deposition sequences.

Required invariants are:

- Solid voxels are never outside-connected.
- Every outside-connected voxel has a six-neighbor path to a source.
- No closed-void voxel has a path to a source.
- Every empty source voxel is outside-connected.
- Incremental and full-recompute classifications are identical.

### MPI equivalence tests

The same global geometry and event stream will be executed with different rank
counts and decompositions. Global checksums and selected gathered test results
must match the serial reference. Gathering is permitted in tests but not in the
normal production algorithm.

### Performance measurements

Performance reports will separate:

- structure input and voxelization;
- initial flood-fill;
- harmless local deposition;
- suspected and confirmed pinch-off repair;
- accessibility query throughput;
- peak memory and bytes per voxel;
- MPI communication volume and scaling.

Atom count, voxel count, voxel spacing, geometry, compiler configuration, rank
count, and hardware must accompany every reported turnaround time.

## 6. Dependency policy

- The initial production core should require only C++17 and CMake/CTest.
- No explicit graph library is needed because the regular grid has implicit
  six-neighbor topology.
- Test and benchmark packages may be added as pinned, test-only dependencies if
  they improve reliability without entering the production library.
- MPI will be an optional build feature beginning in Phase 11.
- HDF5, ADIOS2, or another file library will be considered only after the target
  atomic-structure format and I/O scale are known.
- Any proposed dependency will be reviewed for license, maintenance, portability,
  memory cost, and hot-path overhead before adoption.

## 7. Deliberately deferred functionality

These features are outside the first implementation and require separate phases
if they become necessary:

- etching or other solid-to-gas transitions;
- species-dependent precursor radii or multiple simultaneous gas grids;
- molecular-orientation states, orientation transitions, or
  orientation-dependent transport connectivity;
- adaptive, sparse, or multiresolution voxels unless Phase 10 requires them;
- GPU acceleration;
- diagonal gas connectivity;
- a generic replacement for the KMC application's domain decomposition.

The internal separation between occupancy, connectivity, and querying should
allow these additions without changing the basic reaction-query contract.

## 8. Information checkpoints

Before Phase 0 is finalized:

- atom record layout, units, and radius source;
- physical domain extents;
- expected voxel spacing or desired resolution study;
- effective precursor radius and how it is calibrated for the target molecule;
- periodic axes and reservoir/source definition;
- deposition event payload;
- reaction-site location and desired site-contact rule.

Before Phase 6 uses the supplied structure:

- atomic-structure file format and representative sample;
- target hardware and preferred compiler;
- required output or visualization format.

Before Phase 9:

- confirm that KMC supplies the current atom position at query time;
- confirm that the containing voxel plus six face neighbors is the desired
  atom-to-gas contact stencil.

Before Phase 11:

- Confirmed: reuse SPPARKS `Domain` global/subdomain bounds, `procgrid`, and
  `myloc`; do not create another domain decomposition.
- Confirmed: use exact voxel/subdomain face alignment and integer owned ranges.
- Confirmed: allow per-axis spacing so rectangular cuboid voxels can align
  exactly with each decomposed box axis; equal x/y/z spacing is not required.
- Confirmed: the production off-lattice application is fully periodic, while
  the generic library must also support non-periodic MPI boundaries. The
  GasAccess boundary configuration may override the KMC setting; the expected
  production call makes z non-periodic and uses the top-z reservoir face.
- Confirmed: use the existing `world` communicator for a separate gas-state
  face-halo exchange; atom ghosts are already synchronized before queries.
- Confirmed: implement a defensive GasAccess check that the supplied atom ghost
  distance covers `max(R_atom + R_precursor)` even though KMC should already
  enforce this condition.
- Confirmed: the KMC integration requires only the position-based accessibility
  query and does not need selective invalidation.
- Still required during integration: provide the requested nominal voxel
  resolution; the alignment helper will derive divisible dimensions and report
  the resulting per-axis spacing.
- Still required during integration: expose the minimum off-lattice bin/ghost
  distance and maximum `R_atom + R_precursor` to the defensive validator.
- Still required during integration: confirm whether the bottom-z GasAccess
  face remains closed/non-source; the current planned default is top-z source
  only.
- Confirmed: OpenMPI 4.1.1 is available through `mpicxx` and `mpiexec`; CMake
  discovers MPI when `GASACCESS_ENABLE_MPI=ON`.

## 9. Progress record

- [x] Phase 0: contracts and build scaffold
- [x] Phase 1: voxel grid geometry and topology
- [x] Phase 2: static atom voxelization
- [x] Phase 3: initial exterior classification
- [x] Phase 4: cached reaction-site queries and C interface
- [x] Phase 5: deposition with full recomputation
- [x] Phase 6: standalone reference driver and baseline measurements
- [x] Phase 7: conservative local topology filter
- [x] Phase 8: serial affected-region repair
- [x] Phase 9: single-site KMC query contract
- [x] Phase 10: serial scale and storage optimization
- [x] Phase 11: anisotropic grid, decomposition adapter, and gas ghost exchange
- [x] Phase 12: distributed initial flood-fill
- [x] Phase 13: distributed incremental repair
- [ ] Phase 14: KMC integration and production acceptance

Update this checklist only when a phase's tests and exit gate have passed.

### Completed phase notes

#### Phase 0 completion — 2026-08-12

- Finalized the current C++17 in-memory contracts for grid geometry, atoms,
  sources, states, classification, and bounded site queries.
- Added fixed-width public types, checked configuration validation, documented
  units and naming rules, and status-based error handling for C callers.
- Completed the CMake scaffold and verified both C++ and pure-C clients compile,
  link, and agree on fixed-width state and voxel-ID representations.
- The APIs may still gain backward-compatible adapters when the real KMC
  interfaces are reviewed before Milestone C.

#### Phase 1 completion — 2026-08-11

- Added the CMake C++17 library/test scaffold required to build Phase 1.
- Implemented dense one-byte state storage, checked 64-bit indexing,
  coordinate conversion, six-face neighbors, independent x/y/z wrapping,
  reservoir faces, and explicit source voxels.
- Passed seven deterministic test groups with GCC 8.5 and Clang 20.1 using the
  configured warning set plus `-Werror`.
- Passed Valgrind Memcheck with zero errors and no memory leaks.

#### Phase 2 completion — 2026-08-11

- Added non-owning `AtomView` input with double-precision atom positions and
  radii in the same units as `GridSpec`.
- Implemented additive hard-sphere steric voxelization using
  `R_atom + R_precursor`, bounded per-atom candidate ranges, minimum-image
  distances, and independent periodic wrapping on all axes.
- Added validation for precursor/atom radii, positions, atom-view pointer/count
  consistency, and floating-point grid resolution. All batch inputs are
  validated before the grid is mutated.
- Passed seven Phase 2 test groups, including 96 deterministic randomized
  systems covering all eight periodic-axis combinations and checked against a
  full brute-force atom-versus-voxel oracle.
- Passed the full two-executable CTest suite under GCC 8.5 with the configured
  warnings treated as errors.
- Passed Valgrind Memcheck for the voxelizer suite with zero errors and no
  memory leaks.

#### Phase 3 completion — 2026-08-12

- Added `ExteriorClassifier::classify()` as the permanent full serial reference
  classifier. It preserves `Solid`, resets all cached empty states, seeds every
  non-solid reservoir voxel, and flood-fills through six-face neighbors.
- Added `ClassificationSummary` with solid, outside-accessible, and closed-void
  voxel counts.
- Verified empty, completely solid, sourceless, blocked-source, explicit-source,
  open-trench, sealed-trench, enclosed-cavity, and one-voxel-channel cases.
- Verified edge-only and corner-only contacts remain disconnected and verified
  paths crossing the periodic x, y, and z seams independently.
- Passed state/count/path invariants for 128 deterministic randomized grids
  spanning all eight periodic-axis configurations.
- Passed the full three-executable CTest suite under GCC 8.5 with the configured
  warnings treated as errors.
- Passed Valgrind Memcheck for the classifier suite with zero errors and no
  memory leaks.

#### Phase 4 completion — 2026-08-12

- Added `GasAccessibilityQuery` for direct voxel-state queries, a default
  position-based stencil containing the site voxel plus six face neighbors, and
  a caller-provided stencil capped at 27 voxel IDs.
- Kept gas-to-gas connectivity at six neighbors while allowing the site-contact
  policy to be supplied independently by a future KMC adapter.
- Verified site queries across periodic seams, outside non-periodic boundaries,
  containing voxels, face neighbors, diagonal-only contacts, invalid IDs, empty
  stencils, and oversized stencils.
- Instrumented 3,000 representative hot queries and verified zero dynamic
  allocations and no cached-state mutation.
- Added the opaque-handle C99 API for grid creation/destruction, state access,
  atom voxelization, full exterior classification, and all accessibility query
  forms. C++ exceptions are converted to fixed status codes and a thread-local
  error message.
- Avoided a full atom-array copy in the C wrapper by converting validated atom
  input in fixed-size stack chunks.
- Passed the full five-executable CTest suite under GCC 8.5 with both C and C++
  warnings treated as errors, including a client compiled as pure C99.
- Passed Valgrind Memcheck for both the query suite and pure-C end-to-end API
  test with zero errors and no memory leaks.

#### Phase 5 completion — 2026-08-12

- Added `DepositionUpdater::apply_deposition()` for atoms whose placement and
  deposition physics are determined by the caller/KMC simulator.
- Enforced the monotonic deposition-only contract and required a fully
  classified input grid. Each event batch performs additive steric voxelization
  and, when occupancy changes, exactly one full Phase 3 reclassification.
- Added `DepositionUpdateResult` with the newly solid count, post-update
  classification counts, and sorted unique IDs for every voxel whose state
  changed. Fully overlapping and empty batches skip reclassification.
- Added opaque C update-result handles and accessors without changing ownership
  of the caller's deposited-atom data.
- Verified staged trench sidewall growth and pinch-off, overlapping/no-op
  events, batch validation before mutation, closure across a periodic seam, and
  deposition that blocks the only explicit reservoir source.
- Compared 160 deterministic one- and two-atom event batches across all eight
  periodic-axis configurations against explicit voxelization plus full
  reference classification after every event.
- Passed the full six-executable CTest suite under GCC 8.5 with C and C++
  warnings treated as errors.
- Passed Valgrind Memcheck for the C++ deposition suite and pure-C end-to-end
  update path with zero errors and no memory leaks.

#### Phase 6 completion — 2026-08-12

- Added `gasaccess_reference_driver`, a standalone C++17 executable that uses
  only the public in-memory library API and emits machine-readable `key=value`
  output.
- Added deterministic open-trench, sealed-trench, and seeded bulk generators;
  all accept independent x/y/z periodic settings, grid size, spacing, atom and
  precursor radii, query count, and update count from the command line.
- Reported state summary counts and FNV-1a state checksums before and after
  updates, plus separate timings for grid creation, atom generation,
  voxelization, classification, cached queries, and full deposition updates.
- Reported exact persistent state storage, a conservative traversal-frontier
  logical-payload bound, and Linux process peak RSS with clear accounting
  limitations.
- Added five end-to-end CTest cases covering open, sealed, z-periodic sealed,
  seeded bulk, and fully periodic bulk configurations. Each checks an exact
  deterministic initial-state checksum; the generated sizes progress from 960
  to 12,288 voxels in routine testing.
- Completed and recorded a Release baseline with one million atom records,
  1,048,576 voxels, one million cached queries, and three full-recomputation
  deposition updates. See `docs/benchmarks/PHASE6_BASELINE.md` for commands,
  platform details, timings, counts, checksums, and memory results.
- Passed Valgrind Memcheck for the sealed-trench end-to-end driver path with
  zero errors and no memory leaks.
- Deferred the production structure-file adapter until the file format and a
  representative input are supplied, as required by the information
  checkpoint before this phase uses external structures.

#### Phase 7 completion — 2026-08-14

- Added `LocalTopologyFilter` with a fixed-capacity `3x3x3` six-neighbor search
  that performs no dynamic allocation. It proves a single removed gas voxel
  safe only when all surviving accessible neighbors reconnect locally.
- Treats removal of a reservoir source, a locally disconnected neighbor set,
  and every multi-voxel change as requiring full reference reclassification.
  Blocking an already closed-void voxel and accessible removals with zero or
  one surviving accessible neighbor are safe under the documented classified
  input contract.
- Integrated the filter into `DepositionUpdater`. Proven-safe updates preserve
  cached accessibility and adjust summary counts locally; inconclusive updates
  retain the exact Phase 5 classifier and changed-state diff as a fallback.
- Added `full_reclassification_performed` reporting to the C++ and C deposition
  results and added a matching counter to the standalone reference driver.
- Verified the safe/fallback decision against full reclassification for all
  128 occupancy patterns around a center voxel in a `3x3x1` fixture and for 160
  deterministic event batches spanning all eight periodic-axis combinations.
- Added direct tests for open-space removal, a narrow bridge, explicit-source
  removal, closed-void removal, multi-voxel fallback, a periodic local bypass,
  and a periodic-seam bridge.
- Retained identical initial and final checksums on the million-atom Phase 6
  workload. All three sample updates were proven safe; total update time fell
  from 82.267 ms to 24.625 ms in the single-run comparison recorded in
  `docs/benchmarks/PHASE7_FILTER.md`.
- The updater still performs transitional full-array snapshot and discovery
  scans. Removing those operations is intentionally left to later serial
  storage/update optimization rather than expanding the topology-filter phase.
- Passed the full 12-case CTest suite with GCC 8.5 and C/C++ warnings treated
  as errors. Valgrind Memcheck reported zero errors and no leaks for the local
  filter, deposition updater, and pure-C API suites.

#### Phase 8 completion — 2026-08-14

- Added `AffectedRegionRepair`, which seeds searches only from surviving
  `OutsideAccessible` neighbors of newly solid voxels. A fully explored
  source-free component is relabeled `ClosedVoid`; a search stops early when
  it reaches a reservoir source or a component already proven external.
- Added reusable 32-bit search and repair epochs, a retained BFS frontier, and
  retained seed storage. Ordinary calls avoid clearing visitation arrays, and
  the frontier itself also serves as the discovered-region list.
- Integrated affected-region repair after inconclusive Phase 7 decisions.
  Locally safe events still skip repair, while suspected single- and
  multi-voxel cuts no longer invoke the full classifier in production mode.
- Added `ConnectivityRepairMode::FullReclassification` as a configurable C++
  reference/debug path. The C API exposes the equivalent
  `ga_apply_deposition_with_mode()` operation and reuses updater workspace on
  its grid handle when the radius and mode remain unchanged.
- Added update-result metrics for repair use, visited voxels, and newly closed
  voxels in C++, C, and the standalone reference driver.
- Verified a 21-voxel trench cavity is repaired with exactly 21 BFS visits in
  an 81-voxel grid. Also verified removal of the only source, preservation by a
  second source, a two-voxel cut with one closed middle component, and a
  periodic-seam cut.
- Compared incremental states, summaries, and sorted changed IDs against forced
  full reclassification after every event in the existing 160-event randomized
  suite spanning all eight periodic-axis combinations. The 128-pattern
  exhaustive local suite also remains reference-identical.
- Retained the million-atom initial and final Phase 6 checksums. See
  `docs/benchmarks/PHASE8_REPAIR.md` for the locality and regression record.
- The full state snapshot and newly-solid discovery scan remain transitional;
  Phase 8 eliminates global connectivity traversal for affected repairs, not
  all full-array work in the deposition pipeline.
- Passed the full 13-case CTest suite with GCC 8.5 and C/C++ warnings treated
  as errors. Valgrind Memcheck reported zero errors and no leaks for the repair,
  deposition-updater, and pure-C API suites.

#### Phase 9 completion — 2026-08-17

- Finalized `query.is_site_accessible(atom_position)` as the sole operation
  required from GasAccess inside the KMC atom or reaction-site loop.
- Added a focused mock-KMC integration suite that performs this single call for
  every site and leaves atom storage, reaction state, and iteration under KMC
  ownership.
- Verified that one query object observes an affected-region pinch-off after a
  deposition update without being reconstructed or receiving changed-voxel
  identifiers.
- Verified that KMC can pass a new atom position directly after a simulated MD
  move; GasAccess stores no atom identifier, registry, or position cache.
- Verified periodic wrapping and rejection outside a non-periodic boundary
  through the same position-only query contract.
- Documented the default containing-plus-six-face-neighbor stencil, required
  update ordering, bounded seven-state lookup, and lack of flood-fill, MPI,
  callback, or allocation on the hot query path.
- Passed the full 14-case CTest suite with GCC 8.5 and C/C++ warnings treated as
  errors. Valgrind Memcheck reported zero errors and no leaks for both the KMC
  integration and accessibility-query suites; the existing allocation wrapper
  continued to report zero dynamic allocations across hot queries.

#### Phase 10 completion — 2026-08-17

- Measured the existing million-atom, 1,048,576-voxel workload before changing
  storage or update behavior. Three harmless one-voxel updates required
  23.901399 ms because each update copied, counted, and searched the full grid.
- Added incrementally maintained counts for all four gas states. Classification
  summaries are now available in O(1), with 32 bytes of fixed metadata per grid
  and no change to the one-byte-per-voxel state field. Added matching C++ and C
  state-count accessors with validation tests.
- Added exact newly-solid change capture to `AtomVoxelizer`, including each
  voxel's state immediately before removal. `DepositionUpdater` reuses this
  storage and no longer performs production full-grid snapshots or discovery
  scans.
- Kept full state snapshots and diffs exclusively in the explicitly selected
  `FullReclassification` reference/debug mode.
- Extended the reference driver with per-path update latency distributions and
  a deterministic `pinch-off` scenario. The CI-sized pinch-off reaches the
  exact sealed-trench checksum after one affected-region repair.
- Reduced the same three million-scale harmless updates from 23.901399 ms to
  0.004028 ms in the recorded run while preserving the initial and final state
  checksums. One million cached queries remained approximately 19.76 million
  queries/s.
- Measured 872 locally safe and 128 repair events on a 1,000-update million-atom
  run. Also measured a controlled 39,680-voxel trench closure: 639 safe updates
  had 1.265 us median latency and the final repair took 3.184 ms.
- Retained the dense backend. Persistent state is 1 MiB for the million-voxel
  case; packed state, a compile-time float ABI, and tiled storage were not
  justified by the measured data. Detailed commands, timings, memory, and
  rationale are recorded in
  `docs/benchmarks/PHASE10_SERIAL_OPTIMIZATION.md`.
- Passed the full 15-case CTest suite under GCC 8.5 with C and C++ warnings
  treated as errors. Valgrind Memcheck reported zero errors and no leaks for
  the grid, change-capturing voxelizer, deposition updater, C API, and
  end-to-end pinch-off driver suites.

#### Phase 11 completion — 2026-08-19

- Replaced scalar C++ and C grid spacing with explicit x/y/z components.
  Coordinate mapping, voxel centers, periodic lengths, candidate bounds, the C
  adapter, and reference-driver inputs now handle rectangular cuboid voxels.
- Preserved spherical steric exclusion in physical coordinates and extended
  the randomized brute-force voxelization differential test to non-cubic
  spacing across all eight periodic-axis combinations.
- Added `make_aligned_grid_geometry()` to select near-requested per-axis
  spacing with global dimensions divisible by the supplied MPI process grid.
  Direct specifications remain validation-only and are never silently changed.
- Added the optional `GasAccess::gasaccess_mpi` target with
  `MpiDecompositionSpec`, exact SPPARKS-compatible x-fastest ownership,
  global/local bound validation, integer owned ranges, and independently
  configured GasAccess periodic/reservoir boundaries.
- Added `DistributedGasGrid` with owned-voxel occupancy construction from
  synchronized atom views, one layer of face ghost state, reusable message
  buffers, nonblocking point-to-point halo exchange, periodic self-copy, and no
  communication across a non-periodic global edge.
- Added defensive validation of the declared and per-atom excluded radius
  against the supplied atom ghost distance. Periodic image atoms outside an
  effective non-periodic GasAccess boundary are ignored.
- Added `DistributedGasAccessibilityQuery`; its normal
  `query.is_site_accessible(atom_position)` call reads only owned/face-ghost
  state and performs no MPI operation.
- Added focused one-, two-, and four-rank tests covering aligned ownership,
  x/y/z decompositions, mixed and fully periodic gas boundaries, periodic
  seams, self-neighbors, non-periodic z edges, cross-rank queries, non-cubic
  distributed voxelization, and ghost-distance rejection.
- Passed all 18 registered CTest cases with GCC 8.5 and C/C++ warnings treated
  as errors. Serial grid, voxelizer, and C API Valgrind checks reported no
  errors or definite leaks. The MPI suite passed an invalid-memory Valgrind
  check; full leak reporting showed only process-lifetime allocations rooted
  in this OpenMPI build's `MPI_Init` and component loader.
- Documented the build, SPPARKS field mapping, boundary override, ownership,
  halo, atom coverage, and query-ordering contract in
  `docs/integration/MPI_GRID.md`. Distributed exterior classification remains
  Phase 12 work.

#### Phase 12 completion — 2026-08-19

- Added `DistributedExteriorClassifier`, which preserves owned solid occupancy,
  resets every owned free voxel, seeds configured global reservoir faces and
  explicit source voxels, and performs six-face distributed exterior
  classification.
- Added reusable local BFS storage and compact face-frontier buffers. Messages
  contain only tangential offsets on the receiving face; the implementation
  never gathers or replicates the global gas grid.
- Added count/payload point-to-point exchange with periodic self-neighbor
  handling and an `MPI_Allreduce` used only for global frontier termination.
  Ranks without active local work remain in the collective protocol until all
  ranks are idle.
- Added `DistributedClassificationSummary` with local state counts, visited
  voxel count, communication rounds, and sent/received frontier-entry counts.
- Added final face-halo synchronization inside `classify()` so the existing
  `query.is_site_accessible(atom_position)` contract is immediately valid when
  classification returns.
- Added differential MPI tests against the serial `ExteriorClassifier` and
  serial query for one, two, and four ranks. The fixtures cover x/y/z and 2x2
  decompositions, non-cubic spacing, top-face sources, no-source termination,
  sealed and open barriers, fully periodic explicit-source traversal,
  deterministic obstacles, periodic seams, inactive ranks, state summaries,
  final ghost states, and query results.
- Passed all 21 registered serial, reference-driver, and MPI tests, including
  the new classifier suite on one, two, and four ranks, with strict GCC warnings
  treated as errors. A two-rank Valgrind run using OpenMPI's TCP/self transport
  completed with no GasAccess invalid-memory errors.

#### Phase 13 completion — 2026-08-19

- Added `DistributedDepositionUpdater` with the same affected-region versus
  forced-full mode selection as the serial updater. Calls are collective, while
  the existing position-based KMC query remains communication-free.
- Extended `DistributedGasGrid` with incrementally maintained owned-state
  counts and change-capturing owned atom voxelization. This avoids a local
  full-grid validation/count scan during ordinary incremental events.
- Added a conservative fixed-capacity local topology proof for one globally
  removed accessible voxel. Zero/one-neighbor deletions remain constant-cost;
  larger proofs at rank boundaries escalate when edge/corner ghost data would
  be required.
- Added the safe fast paths for no geometry change and removal entirely from an
  existing closed void. Multi-voxel removals involving accessible gas
  conservatively invoke distributed repair.
- Added rare-path synchronization of only removed-voxel coordinates and prior
  states. No global gas-state field or gas graph is gathered.
- Added distributed affected-component searches with reusable epoch arrays,
  local BFS storage, compact tangential face offsets, periodic self-neighbor
  handling, source/termination reductions, and owned relabelling of newly
  closed components.
- Added local changed-coordinate diagnostics and traversal/communication
  metrics. Final face ghosts are synchronized on every geometry-changing
  incremental path.
- Added forced Phase 12 reclassification as a distributed debug/reference mode
  with exact changed-owned-coordinate comparison.
- Added one-, two-, and four-rank differential tests that maintain incremental,
  forced-full, and serial grids together. Coverage includes no-change and
  closed-void fast paths, proven-safe local removal, local and rank-boundary
  pinch-offs, periodic-seam closure, multi-voxel closure, spanning cavities,
  final ghosts/queries/counts, and a deterministically shuffled 64-event
  deposition sequence.
- Passed all 24 registered serial, reference-driver, and MPI tests with strict
  GCC warnings treated as errors. A two-rank Valgrind run of the complete
  distributed updater suite using OpenMPI's TCP/self transport completed with
  no GasAccess invalid-memory errors.
