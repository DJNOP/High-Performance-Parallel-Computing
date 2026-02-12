# Water Molecular Dynamics - SoA Implementation
## Code Walkthrough for Video Presentation

---

## 1. NEW DATA STRUCTURES (Lines 115-144)

### Atoms Class - Holds N atoms of same type (SoA)
```cpp
class Atoms {
public:
    double mass, ep, sigma, charge;  // Properties shared by all atoms of this type
    std::string name;
    std::vector<Vec3> p, v, f;       // Arrays sized N_molecules
    
    // Constructor allocates arrays for N_molecules
    Atoms(double mass, double ep, double sigma, double charge, 
          std::string name, size_t N_identical) 
    : mass{mass}, ep{ep}, sigma{sigma}, charge{charge}, name{name}, 
      p{N_identical, {0,0,0}},  // All positions
      v{N_identical, {0,0,0}},  // All velocities
      f{N_identical, {0,0,0}}   // All forces
    {}
};
```
**KEY POINT:** Instead of 1 atom, we store N atoms' data in parallel arrays

### Molecules Class - Container for all molecules (SoA)
```cpp
class Molecules {
public:
    std::vector<Atoms> atoms;                      // 3 blocks: O, H1, H2
    std::vector<Bond> bonds;                       // Shared topology
    std::vector<Angle> angles;                     // Shared topology
    std::vector<std::vector<size_t>> neighbours;   // Per-molecule neighbors
    size_t no_mol;                                 // Number of molecules
};
```
**KEY POINT:** All water molecules share same topology, so we only need 3 Atoms blocks

---

## 2. NEIGHBOR LIST (Lines 194-234)

### Critical Access Pattern Change
```cpp
void BuildNeighborList(System& sys){
    Molecules& mol = sys.molecules;
    
    for(size_t i = 0; i < mol.no_mol; i++){
        for (size_t j = 0; j < mol.no_mol; j++) {
            // OLD: sys.molecules[i].atoms[0].p
            // NEW: mol.atoms[0].p[i]
            Vec3 dp = mol.atoms[0].p[i] - mol.atoms[0].p[j];
            distances2[j] = dp.mag2();
        }
        // ... sorting and neighbor assignment ...
    }
}
```
**KEY POINT:** Oxygen positions are now in a contiguous array: atoms[0].p[0], atoms[0].p[1], ...

---

## 3. BOND FORCES (Lines 236-254)

### Loop Reversal - The Heart of SoA
```cpp
void UpdateBondForces(System& sys){
    Molecules& mol = sys.molecules;
    
    // OLD PATTERN (AoS):
    // for (Molecule& molecule : sys.molecules)     // Outer: N molecules
    //     for (Bond& bond : molecule.bonds)        // Inner: 2 bonds
    
    // NEW PATTERN (SoA):
    for (Bond& bond : mol.bonds){                   // Outer: 2 bond types
        Atoms& A1 = mol.atoms[bond.a1];            // Get atom type 1 block
        Atoms& A2 = mol.atoms[bond.a2];            // Get atom type 2 block
        
        for (size_t i = 0; i < mol.no_mol; i++){   // Inner: N molecules
            Vec3 dp  = A1.p[i] - A2.p[i];          // Access consecutive memory
            Vec3 f   = -bond.K*(1 - bond.L0/dp.mag())*dp;
            A1.f[i] += f;                           // Update consecutive memory
            A2.f[i] -= f;
        }
    }
}
```
**KEY POINT:** Inner loop now accesses p[0], p[1], p[2]... which are consecutive in memory
**VECTORIZATION:** Compiler can process 4-8 molecules simultaneously with SIMD

---

## 4. ANGLE FORCES (Lines 256-292)

### Same Pattern with 3 Atom Types
```cpp
void UpdateAngleForces(System& sys){
    Molecules& mol = sys.molecules;
    
    for (Angle& angle : mol.angles){               // Just 1 angle type (H-O-H)
        Atoms& A1 = mol.atoms[angle.a1];          // H1 atoms
        Atoms& A2 = mol.atoms[angle.a2];          // O atoms
        Atoms& A3 = mol.atoms[angle.a3];          // H2 atoms
        
        for (size_t i = 0; i < mol.no_mol; i++){   // All molecules
            // Calculate angle between A2.p[i]-A1.p[i] and A2.p[i]-A3.p[i]
            Vec3 d21 = A2.p[i] - A1.p[i];
            Vec3 d23 = A2.p[i] - A3.p[i];
            // ... angle force calculation ...
            A1.f[i] += f1;
            A2.f[i] -= (f1 + f3);
            A3.f[i] += f3;
        }
    }
}
```
**KEY POINT:** Loop over 1 angle type, then vectorize over N molecules

---

## 5. NON-BONDED FORCES (Lines 294-325)

### SoA Indexing in Nested Loops
```cpp
void UpdateNonBondedForces(System& sys){
    Molecules& mol = sys.molecules;
    
    for (size_t i = 0; i < mol.no_mol; i++)          // Molecule i
    for (auto& j : mol.neighbours[i])                 // Neighbor j
    for (size_t a = 0; a < mol.atoms.size(); a++)     // Atom type a (0=O, 1=H1, 2=H2)
    for (size_t b = 0; b < mol.atoms.size(); b++){    // Atom type b
        
        Atoms& atom1 = mol.atoms[a];
        Atoms& atom2 = mol.atoms[b];
        
        // OLD: atom1.p - atom2.p
        // NEW: atom1.p[i] - atom2.p[j]
        Vec3 dp = atom1.p[i] - atom2.p[j];
        
        // Calculate Lennard-Jones + Coulomb forces
        Vec3 f = (LJ_term + Coulomb_term) * dp;
        
        atom1.f[i] += f;    // Force on molecule i, atom type a
        atom2.f[j] -= f;    // Equal/opposite on molecule j, atom type b
    }
}
```
**KEY POINT:** Access pattern is atoms[type].p[molecule] instead of molecules[i].atoms[type].p

---

## 6. LEAP-FROG INTEGRATOR (Lines 327-341)

### Perfect Vectorization Target
```cpp
void UpdateKDK(System &sys, Sim_Configuration &sc){
    Molecules& mol = sys.molecules;
    
    // OLD PATTERN (AoS):
    // for (Molecule& molecule : sys.molecules)
    //     for (auto& atom : molecule.atoms)
    
    // NEW PATTERN (SoA):
    for (auto& atom : mol.atoms)                      // 3 atom types
    for (size_t i = 0; i < mol.no_mol; i++){         // N molecules
        
        atom.v[i] += sc.dt/atom.mass*atom.f[i];      // v[0], v[1], v[2]... consecutive
        atom.f[i]  = {0,0,0};                         // f[0], f[1], f[2]... consecutive
        atom.p[i] += sc.dt*atom.v[i];                // p[0], p[1], p[2]... consecutive
    }
    
    sys.time += sc.dt;
}
```
**KEY POINT:** This is the most important loop for vectorization!
**PERFORMANCE:** Inner loop operates on 3 * N_molecules contiguous array elements

---

## 7. INITIALIZATION (Lines 343-384)

### Creating SoA Structure
```cpp
System MakeWater(size_t N_molecules){
    // Create 3 Atoms blocks, each sized N_molecules
    Atoms O (16, 0.65,    0.31,  -0.82, "O",  N_molecules);  // All oxygens
    Atoms H1( 1, 0.18828, 0.238,  0.41, "H",  N_molecules);  // All H1 atoms
    Atoms H2( 1, 0.18828, 0.238,  0.41, "H",  N_molecules);  // All H2 atoms
    
    // Topology is shared (same for all molecules)
    std::vector<Bond> waterbonds = { bond1, bond2 };
    std::vector<Angle> waterangle = { angle1 };
    
    Molecules mols({O, H1, H2}, waterbonds, waterangle, N_molecules);
    
    // Initialize positions by molecule index
    for (size_t i = 0; i < N_molecules; i++){
        Vec3 P0 = calculate_sphere_position(i);
        
        // OLD: Oatom.p = P0; sys.molecules.push_back(...)
        // NEW: Direct array assignment
        mols.atoms[0].p[i] = P0;                       // Oxygen
        mols.atoms[1].p[i] = P0 + hydrogen1_offset;    // H1
        mols.atoms[2].p[i] = P0 + hydrogen2_offset;    // H2
    }
    
    return System(std::move(mols));
}
```
**KEY POINT:** Allocate all arrays upfront, then fill by index

---

## 8. OUTPUT (Lines 386-400)

### Preserve Output Order
```cpp
void WriteOutput(System& sys, std::ofstream& file){
    Molecules& mol = sys.molecules;
    
    // Must write in same order as AoS version: molecule-by-molecule
    for (size_t i = 0; i < mol.no_mol; i++)          // Molecule i
    for (size_t a = 0; a < mol.atoms.size(); a++){   // Atom types O, H1, H2
        auto& atom = mol.atoms[a];
        file << sys.time << " " << atom.name << " "
             << atom.p[i].x << " "    // Position of atom type a in molecule i
             << atom.p[i].y << " "
             << atom.p[i].z << '\n';
    }
}
```
**KEY POINT:** Output order unchanged (molecule-by-molecule), just different access pattern

---

## SUMMARY OF TRANSFORMATION

### Memory Layout Change:
```
AoS: [Mol0: O,H1,H2][Mol1: O,H1,H2][Mol2: O,H1,H2]...
     ^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^
     molecule 0      molecule 1      molecule 2

SoA: [O: mol0,mol1,mol2,...][H1: mol0,mol1,mol2,...][H2: mol0,mol1,mol2,...]
     ^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^
     all oxygens              all hydrogen-1           all hydrogen-2
```

### Access Pattern Change:
```
AoS: sys.molecules[i].atoms[type].property
                   ^        ^
                   mol      atom

SoA: sys.molecules.atoms[type].property[i]
                        ^               ^
                        atom            mol
```

### Loop Nesting Change:
```
AoS: for molecule → for atom_type     (outer: molecules)
SoA: for atom_type → for molecule     (outer: atom types, VECTORIZABLE inner loop)
```

### Performance Benefits:
1. **Cache efficiency:** Consecutive access to p[i], p[i+1], p[i+2]...
2. **Vectorization:** Inner loops can use SIMD (4-8x speedup potential)
3. **Memory bandwidth:** Better utilization of cache lines

### Correctness Verification:
Both versions produce identical checksums, confirming the transformation is correct!
