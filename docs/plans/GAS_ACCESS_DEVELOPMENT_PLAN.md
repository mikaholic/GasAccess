# GasAccess Development Plan

Status: planning  
Last updated: 2026-08-11

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

- Uniform Cartesian voxel grid.
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
    bool is_site_accessible(const SiteProbe& site) const;
    UpdateResult apply_deposition(AtomView new_atoms);
};
```

`UpdateResult` will eventually expose changed voxel identifiers and, when site
registration is enabled, affected reaction-site identifiers. The permanent
reference API will also provide an explicit full recomputation for validation.

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

#### Phase 9: selective reaction-site invalidation

Scope:

- Return stable identifiers for voxels whose accessibility changed.
- Add optional registration of reaction sites and their nearby gas voxels.
- Maintain a reverse voxel-to-site mapping or an integration callback so only
  affected propensities need recomputation.
- Avoid forcing KMC-specific types into the core library.

Tests:

- Changed sites match a full scan of all registered sites.
- Unaffected sites are not reported.
- Multiple changed voxels do not duplicate site notifications.
- Periodic-boundary sites map to the correct wrapped voxels.

Exit gate:

- The library can drive selective KMC propensity updates without global scans.

Expected production code: 250-450 lines.

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

Milestone C does not begin until the existing decomposition interfaces are
provided. Required information includes ownership rules, ghost width and
exchange API, Cartesian rank neighbors, periodic rank topology, event ownership,
global/local coordinate conversion, and propensity invalidation interfaces.

#### Phase 11: decomposition adapter and gas ghost exchange

Scope:

- Map global gas voxels to owned and ghost storage using the KMC decomposition.
- Exchange occupancy and accessibility flags with spatial neighbor ranks.
- Wrap MPI-global edges consistently on every periodic axis.
- Do not implement distributed flood-fill yet.

Tests:

- One-rank behavior matches the serial library.
- Owned/ghost states match on two- and four-rank decompositions.
- Tests cover internal rank interfaces and global periodic seams.

Exit gate:

- Every rank has correct local and ghost gas state after exchange.

Expected production code: 300-550 lines, depending on the supplied adapter API.

#### Phase 12: distributed initial flood-fill

Scope:

- Maintain a local frontier per rank.
- Send frontier entries to the owning neighbor rank.
- Use a small collective only to detect global termination.
- Avoid global gas-graph replication and routine all-gather.

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

- Extend topology checks across owned/ghost boundaries.
- Exchange connectivity-repair frontiers only when required.
- Detect distributed termination.
- Relabel affected components, synchronize ghosts, and report changed local
  sites.
- Retain distributed full recomputation as a debug/fallback mode.

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

- Adapt real atom, site, deposition, domain, and propensity interfaces.
- Build/initialize the grid from the KMC structure.
- Invoke updates only after geometry-changing events.
- Recalculate only affected reaction propensities.
- Run the supplied million-atom structure and representative event sequence.
- Measure strong scaling, communication volume, update latency, query throughput,
  and end-to-end KMC overhead.

Tests:

- A mock KMC integration test precedes changes to the real simulator.
- Integrated results match the standalone library on the same geometry/events.
- Normal site queries perform no collectives and require only local/ghost state.
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

Before Phase 11:

- KMC decomposition and ghost-exchange interfaces;
- MPI topology and ownership rules;
- local/global atom and site identifiers;
- existing propensity-update mechanism.

## 9. Progress record

- [ ] Phase 0: contracts and build scaffold
- [x] Phase 1: voxel grid geometry and topology
- [ ] Phase 2: static atom voxelization
- [ ] Phase 3: initial exterior classification
- [ ] Phase 4: cached reaction-site queries and C interface
- [ ] Phase 5: deposition with full recomputation
- [ ] Phase 6: standalone reference driver and baseline measurements
- [ ] Phase 7: conservative local topology filter
- [ ] Phase 8: serial affected-region repair
- [ ] Phase 9: selective reaction-site invalidation
- [ ] Phase 10: serial scale and storage optimization
- [ ] Phase 11: decomposition adapter and gas ghost exchange
- [ ] Phase 12: distributed initial flood-fill
- [ ] Phase 13: distributed incremental repair
- [ ] Phase 14: KMC integration and production acceptance

Update this checklist only when a phase's tests and exit gate have passed.

### Completed phase notes

#### Phase 1 completion — 2026-08-11

- Added the CMake C++17 library/test scaffold required to build Phase 1.
- Implemented dense one-byte state storage, checked 64-bit indexing,
  coordinate conversion, six-face neighbors, independent x/y/z wrapping,
  reservoir faces, and explicit source voxels.
- Passed seven deterministic test groups with GCC 8.5 and Clang 20.1 using the
  configured warning set plus `-Werror`.
- Passed Valgrind Memcheck with zero errors and no memory leaks.
- Phase 0 remains open because its complete C ABI and final input contracts are
  intentionally scheduled for review before later functionality depends on
  them.
