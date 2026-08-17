# KMC single-site query contract

GasAccess requires one operation inside the KMC atom or reaction-site loop:

```cpp
const bool accessible = query.is_site_accessible(atom_position);
```

`atom_position` is the current position supplied by KMC. If MD relaxation moves
an atom, KMC passes the updated position on its next iteration; GasAccess does
not register atoms or cache their positions.

The position query checks the voxel containing the atom and its six
face-neighbor voxels. It returns `true` when at least one checked voxel has
`GasState::OutsideAccessible`. It returns `false` when none do, including when
the neighboring free space is a closed void. A position outside a
non-periodic boundary is inaccessible. A position outside a periodic boundary
is wrapped onto the corresponding axis.

A query may be constructed once and reused:

```cpp
gasaccess::GasAccessibilityQuery query(gas_grid);

for (auto& atom : atoms) {
    const bool accessible = query.is_site_accessible(atom.position);
    atom.gas_reaction_rate = accessible
        ? calculate_gas_reaction_rate(atom)
        : 0.0;
}
```

`GasAccessibilityQuery` holds a reference to `gas_grid`, so a reused query sees
the latest completed connectivity update. KMC must finish any geometry and gas
connectivity update before starting the site loop. Concurrent writes to the
grid while it is being queried are outside this serial contract.

The query performs no flood-fill, MPI communication, atom lookup, callback, or
dynamic allocation. Its default stencil inspects at most seven cached voxel
states, so the cost is bounded and independent of the total grid or atom count.

Changed-voxel output and caller-defined stencil queries remain available as
auxiliary APIs, but the KMC integration does not need them.
