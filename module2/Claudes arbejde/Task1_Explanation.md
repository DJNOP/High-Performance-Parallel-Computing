# Task 1: Struct-of-Arrays (SoA) Transformation

## Overview
This document explains the transformation from Array-of-Structures (AoS) to Struct-of-Arrays (SoA) layout in the molecular dynamics water simulation code.

## Key Changes Made

### 1. **Atom → Atoms Class Transformation**

**Original (AoS):**
```cpp
class Atom {
public:
    double mass, ep, sigma, charge;
    std::string name;
    Vec3 p, v, f;  // Single atom's position, velocity, force
};
```

**New (SoA):**
```cpp
class Atoms {
public:
    double mass, ep, sigma, charge;  // Same for all atoms of this type
    std::string name;
    std::vector<Vec3> p, v, f;  // Arrays for all molecules' atoms of this type
};
```

**Rationale:** Instead of one Atom per molecule, we now have one Atoms block containing data for ALL molecules' atoms of the same type (e.g., all oxygen atoms, all hydrogen1 atoms, all hydrogen2 atoms). This creates contiguous memory layout where p[0], p[1], p[2]... are sequential in memory.

### 2. **Molecule → Molecules Class Transformation**

**Original (AoS):**
```cpp
class Molecule {
    std::vector<Atom> atoms;       // 3 atoms per molecule
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<size_t> neighbours;
};

class System {
    std::vector<Molecule> molecules;  // N molecules
};
```

**New (SoA):**
```cpp
class Molecules {
    std::vector<Atoms> atoms;  // 3 Atoms blocks (O, H1, H2)
    std::vector<Bond> bonds;   // Shared topology
    std::vector<Angle> angles; // Shared topology
    std::vector<std::vector<size_t>> neighbours;  // One list per molecule
    size_t no_mol;
};

class System {
    Molecules molecules;  // Single SoA container
};
```

**Rationale:** Since all water molecules have identical topology (same bonds, same angles), we can share the topology description and just store N copies of the atomic data in arrays indexed by molecule number.

### 3. **Loop Restructuring**

The critical change is reversing the loop order to make the innermost loop iterate over molecule indices.

**Original Pattern (AoS):**
```cpp
for (Molecule& molecule : sys.molecules)          // Outer: molecules
    for (auto& atom : molecule.atoms) {           // Inner: atoms
        atom.v += dt/atom.mass * atom.f;
        atom.p += dt * atom.v;
    }
```

**New Pattern (SoA):**
```cpp
Molecules& mol = sys.molecules;
for (auto& atom : mol.atoms)                      // Outer: atom types
    for (size_t i = 0; i < mol.no_mol; i++) {    // Inner: molecules
        atom.v[i] += dt/atom.mass * atom.f[i];
        atom.p[i] += dt * atom.v[i];
    }
```

**Rationale:** The innermost loop now accesses consecutive array elements (atom.v[0], atom.v[1], atom.v[2]...), which:
- Improves cache locality (entire cache line is used)
- Enables SIMD vectorization (process multiple molecules simultaneously)
- Reduces cache misses dramatically

### 4. **Detailed Function Transformations**

#### **BuildNeighborList**
- **Change:** Use `mol.atoms[0].p[i]` instead of `sys.molecules[i].atoms[0].p`
- **Reason:** Oxygen positions now stored in contiguous array; molecule i is the i-th element

#### **UpdateBondForces**
```cpp
// Old: Loop over molecules, then bonds
for (Molecule& molecule : sys.molecules)
    for (Bond& bond : molecule.bonds) { ... }

// New: Loop over bonds, then molecules
for (Bond& bond : mol.bonds) {
    Atoms& A1 = mol.atoms[bond.a1];
    Atoms& A2 = mol.atoms[bond.a2];
    for (size_t i = 0; i < mol.no_mol; i++) {
        Vec3 dp = A1.p[i] - A2.p[i];
        Vec3 f = -bond.K*(1 - bond.L0/dp.mag())*dp;
        A1.f[i] += f;
        A2.f[i] -= f;
    }
}
```
**Key insight:** Since bond topology is identical for all molecules, we loop over the 2 bond types once, then vectorize over all molecules.

#### **UpdateAngleForces**
```cpp
// Similar pattern: loop angle types, then molecules
for (Angle& angle : mol.angles) {
    Atoms& A1 = mol.atoms[angle.a1];
    Atoms& A2 = mol.atoms[angle.a2];
    Atoms& A3 = mol.atoms[angle.a3];
    for (size_t i = 0; i < mol.no_mol; i++) {
        // Angle force calculation using A1.p[i], A2.p[i], A3.p[i]
    }
}
```

#### **UpdateNonBondedForces**
```cpp
// Old: 4-level nesting over molecules and atoms
for (size_t i : molecules)
    for (size_t j : neighbours[i])
        for (Atom& atom1 : molecules[i].atoms)
            for (Atom& atom2 : molecules[j].atoms)

// New: 4-level nesting but with SoA indexing
for (size_t i = 0; i < mol.no_mol; i++)
    for (auto& j : mol.neighbours[i])
        for (size_t a = 0; a < mol.atoms.size(); a++)
            for (size_t b = 0; b < mol.atoms.size(); b++) {
                Vec3 dp = mol.atoms[a].p[i] - mol.atoms[b].p[j];
                // ... force calculation
                mol.atoms[a].f[i] += f;
                mol.atoms[b].f[j] -= f;
            }
```
**Key change:** Access pattern is now `atoms[type].p[molecule_index]` instead of `molecules[index].atoms[type].p`

#### **UpdateKDK (Leap-Frog Integrator)**
```cpp
for (auto& atom : mol.atoms)
    for (size_t i = 0; i < mol.no_mol; i++) {
        atom.v[i] += sc.dt/atom.mass*atom.f[i];
        atom.f[i]  = {0,0,0};
        atom.p[i] += sc.dt*atom.v[i];
    }
```
**Vectorization opportunity:** Inner loop processes `mol.no_mol` independent operations on contiguous arrays.

#### **MakeWater (Initialization)**
```cpp
// Create 3 Atoms blocks sized N_molecules
Atoms O(16, 0.65, 0.31, -0.82, "O", N_molecules);
Atoms H1(1, 0.18828, 0.238, 0.41, "H", N_molecules);
Atoms H2(1, 0.18828, 0.238, 0.41, "H", N_molecules);

// Fill positions by molecule index
for (size_t i = 0; i < N_molecules; i++) {
    Vec3 P0 = calculate_position(i);
    mols.atoms[0].p[i] = P0;                     // Oxygen
    mols.atoms[1].p[i] = P0 + offset1;           // Hydrogen 1
    mols.atoms[2].p[i] = P0 + offset2;           // Hydrogen 2
}
```

#### **WriteOutput**
```cpp
// Old: Loop molecules, then atoms
for (Molecule& molecule : sys.molecules)
    for (auto& atom : molecule.atoms)
        file << atom.p.x << " " << atom.p.y << " " << atom.p.z;

// New: Loop molecules, then atom types
for (size_t i = 0; i < mol.no_mol; i++)
    for (size_t a = 0; a < mol.atoms.size(); a++)
        file << mol.atoms[a].p[i].x << " " 
             << mol.atoms[a].p[i].y << " " 
             << mol.atoms[a].p[i].z;
```

## Memory Layout Comparison

### Array-of-Structures (Original):
```
Memory: [Mol0: O,H1,H2][Mol1: O,H1,H2][Mol2: O,H1,H2]...
```
- When processing oxygen atoms, must skip over H1 and H2 data
- Cache line contains mixed data (O, H1, H2)
- Cannot vectorize across molecules easily

### Struct-of-Arrays (New):
```
Memory: [O0,O1,O2,...][H1_0,H1_1,H1_2,...][H2_0,H2_1,H2_2,...]
```
- All oxygen positions are consecutive
- Cache line fully utilized when processing one atom type
- Perfect for SIMD: can process 4-8 molecules at once

## Performance Benefits

1. **Cache Efficiency:**
   - AoS: Loading molecule i brings unwanted atom types into cache
   - SoA: Loading atom[i] brings neighboring molecules' data (which we need next)

2. **Vectorization Potential:**
   - Inner loops now operate on contiguous arrays
   - Compiler can auto-vectorize with `-fopenmp-simd`
   - Can process multiple molecules per CPU instruction

3. **Memory Bandwidth:**
   - Better utilization of cache lines
   - Fewer cache misses
   - Prefetching more effective

## Verification

Both versions produce identical checksums:
```
Accumulated forces Bonds   : 7.4899e+07
Accumulated forces Angles  : 4.7449e+07
Accumulated forces Non-bond: 9.8363e+08
```

This confirms the transformation preserves correctness while enabling better performance.

## Summary

The AoS → SoA transformation involved:
1. Converting `Atom` (single) → `Atoms` (array of N)
2. Converting `vector<Molecule>` → single `Molecules` container
3. Reversing loop nesting: atom types outer, molecule index inner
4. Changing access pattern: `molecules[i].atoms[type]` → `atoms[type].data[i]`
5. Maintaining identical topology (bonds, angles) shared across molecules

This enables the compiler to vectorize the inner loops, dramatically improving performance for large systems while maintaining identical numerical results.
