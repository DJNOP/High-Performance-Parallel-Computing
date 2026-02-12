#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cassert>
#include <math.h>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <string> // needed for std::string (some compilers won't pull it in via iostream)

const double deg2rad = acos(-1)/180.0;    // pi/180 for changing degs to radians
double accumulated_forces_bond  = 0.;     // Checksum: accumulated size of forces
double accumulated_forces_angle = 0.;     // Checksum: accumulated size of forces
double accumulated_forces_non_bond = 0.;  // Checksum: accumulated size of forces
constexpr size_t nClosest = 8;            // Number of closest neighbors to consider in neighbor list. 
// constexpr just helps the compiler to optimize the code better

class Vec3 {
public:
    double x, y, z;
    Vec3(double x, double y, double z): x(x), y(y), z(z) {} // initialization of vector
    double mag2() const{                                    // squared size of vector (slightly faster)
        return x*x+y*y+z*z;
    }
    double mag() const{                                     // size of vector
    return sqrt(mag2());
    }
    Vec3 operator-(const Vec3& other) const{                // subtraction of two vectors
        return {x - other.x, y - other.y, z - other.z};
    }
    Vec3 operator+(const Vec3& other) const{                // addition of two vectors
        return {x + other.x, y + other.y, z + other.z};
    }
    Vec3 operator*(double scalar) const{                   // multiplication of vector by scalar (vec x scalar)
        return {scalar*x, scalar*y, scalar*z};
    }
    Vec3 operator/(double scalar) const{                  // division of vector by scalar
        return {x/scalar, y/scalar, z/scalar};
    }
    Vec3& operator+=(const Vec3& other){                 // add and assign to vector
        x += other.x; y += other.y; z += other.z;
        return *this;
    }
    Vec3& operator-=(const Vec3& other){                // subtract and assign to vector
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }
    Vec3& operator*=(double scalar){                    // multiply and assign to vector
        x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }
    Vec3& operator/=(double scalar){                    // divide and assign to vector
        x /= scalar; y /= scalar; z /= scalar;
        return *this;
    }
};
Vec3 operator*(double scalar, const Vec3& y){           // multiplication of scalar by vector (scalar x vec)
    return y*scalar;
}
Vec3 cross(const Vec3& a, const Vec3& b){
    return { a.y*b.z-a.z*b.y,
             a.z*b.x-a.x*b.z,
             a.x*b.y-a.y*b.x };
}
double dot(const Vec3& a, const Vec3& b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/* atom class (AoS version kept for reference / original layout) */
class Atom {
public:
    double mass;      // The mass of the atom in (U)
    double ep;        // epsilon for LJ potential
    double sigma;     // Sigma, somehow the size of the atom
    double charge;    // charge of the atom (partial charge)
    std::string name; // Name of the atom
    // the position in (nm), velocity (nm/ps) and forces (k_BT/nm) of the atom
    Vec3 p,v,f;
    // constructor, takes parameters and allocates p, v and f properly
    Atom(double mass, double ep, double sigma, double charge, std::string name) 
    : mass{mass}, ep{ep}, sigma{sigma}, charge{charge}, name{name}, p{0,0,0}, v{0,0,0}, f{0,0,0}
    {}
};

/* class for the covalent bond between two atoms U=0.5k(r12-L0)^2 */
class Bond {
public:
    double K;      // force constant
    double L0;     // relaxed length
    size_t a1, a2; // the indexes of the atoms at either end
};

/* class for the angle between three atoms U=0.5K(phi123-phi0)^2 */
class Angle {
public:
    double K;
    double Phi0;
    size_t a1, a2, a3; // the indexes of the three atoms, with a2 being the centre atom
};

/* molecule class (AoS version kept for reference / original layout) */
class Molecule {
public:
    std::vector<Atom> atoms;          // list of atoms in the molecule
    std::vector<Bond> bonds;          // the bond potentials, eg for water the left and right bonds
    std::vector<Angle> angles;        // the angle potentials, for water just the single one, but keep it a list for generality
    std::vector<size_t> neighbours;   // indices of the neighbours
};

// ===============================================================================
// Task 1: SoA layout (Structure-of-Arrays) so data for same atom-type is contiguous
// ===============================================================================

/* atoms class, representing N instances of identical atoms (SoA block) */
class Atoms {
public:
    double mass;      // The mass of the atom in (U)
    double ep;        // epsilon for LJ potential
    double sigma;     // Sigma, somehow the size of the atom
    double charge;    // charge of the atom (partial charge)
    std::string name; // Name of the atom
    // SoA: p[i], v[i], f[i] are for molecule i (same atom-type)
    std::vector<Vec3> p,v,f;

    Atoms(double mass, double ep, double sigma, double charge, std::string name, size_t N_identical) 
    : mass{mass}, ep{ep}, sigma{sigma}, charge{charge}, name{name}, 
      p{N_identical, {0,0,0}}, v{N_identical, {0,0,0}}, f{N_identical, {0,0,0}}
    {}
};

/* molecules container (SoA): holds atom blocks + topology + neighbour lists */
class Molecules {
public:
    std::vector<Atoms> atoms;                     // atom blocks (O, H1, H2)
    std::vector<Bond> bonds;                      // same bonds for each molecule
    std::vector<Angle> angles;                    // same angles for each molecule
    std::vector<std::vector<size_t>> neighbours;  // neighbours[i] for molecule i
    size_t no_mol;

    Molecules(std::vector<Atoms> atoms, std::vector<Bond> bonds, std::vector<Angle> angles, size_t no_mol)
    : atoms{atoms}, bonds{bonds}, angles{angles}, neighbours{no_mol}, no_mol{no_mol}
    {}
};

// ===============================================================================

/* system class */
class System {
public:
    Molecules molecules; 
    double time = 0;

    // changed: System stores one SoA "Molecules" container instead of vector<Molecule>
    // reason: forces/integrator will index into contiguous arrays (better cache/SIMD)
    explicit System(Molecules mols) : molecules(std::move(mols)) {}
};

class Sim_Configuration {
public:
    size_t steps = 10000;     // number of steps
    size_t no_mol = 100;      // number of molecules
    double dt = 0.0005;       // integrator time step
    size_t data_period = 100; // how often to save coordinate to trajectory
    std::string filename = "trajectory.txt";   // name of the output file with trajectory

    Sim_Configuration(std::vector <std::string> argument){
        for (size_t i = 1; i<argument.size() ; i += 2){
            std::string arg = argument.at(i);
            if(arg=="-h"){ // Write help
                std::cout << "MD -steps <number of steps> -no_mol <number of molecules>"
                          << " -fwrite <io frequency> -dt <size of timestep> -ofile <filename> \n";
                exit(0);
                break;
            } else if(arg=="-steps"){
                steps = std::stoi(argument[i+1]);
            } else if(arg=="-no_mol"){
                no_mol = std::stoi(argument[i+1]);
            } else if(arg=="-fwrite"){
                data_period = std::stoi(argument[i+1]);
            } else if(arg=="-dt"){
                dt = std::stof(argument[i+1]);
            } else if(arg=="-ofile"){
                filename = argument[i+1];
            } else{
                std::cout << "---> error: the argument type is not recognized \n";
            }
        }

        dt /= 1.57350; /// convert to ps based on having energy in k_BT, and length in nm
    }
};

// Build neighbour list (SoA version)
// changed: use mol.no_mol and mol.atoms[0].p[i] (oxygen positions) instead of sys.molecules[i].atoms[0].p
// reason: molecules are identified by index i; oxygen array is contiguous and represents molecule position well
void BuildNeighborList(System& sys){
    Molecules& mol = sys.molecules;
    size_t N = mol.no_mol;

    std::vector<double> distances2(N);
    std::vector<size_t> index(N);

    for(size_t i = 0; i < N; i++){
        mol.neighbours[i].clear();

        for (size_t j = 0; j < N; j++) {
            Vec3 dp = mol.atoms[0].p[i] - mol.atoms[0].p[j]; // SoA: O position in molecule i/j
            distances2[j] = dp.mag2();
            index[j] = j;
        }
        distances2[i] = 1e99;

        size_t target_num = std::min(nClosest, N - 1);

        auto lambda_compare = [&](size_t &a, size_t &b) { return distances2[a] < distances2[b]; };

        std::partial_sort(index.begin(),
                          index.begin() + target_num,
                          index.end(),
                          lambda_compare);

        for (size_t jj = 0; jj < target_num; jj++) {
            size_t k = index[jj];
            if (k < i) {
                auto& neighbours_k = mol.neighbours[k];
                if (std::find(neighbours_k.begin(), neighbours_k.end(), i) == neighbours_k.end())
                    mol.neighbours[i].push_back(k);
            } else {
                mol.neighbours[i].push_back(k);
            }
        }
    }
}

// Bond forces (SoA version)
// changed: loop bonds once, then loop molecule index i and update A1.f[i], A2.f[i]
// reason: bond topology is identical for every water molecule; SoA keeps A1/A2 arrays contiguous
void UpdateBondForces(System& sys){
    Molecules& mol = sys.molecules;

    for (Bond& bond : mol.bonds){
        Atoms& A1 = mol.atoms[bond.a1];
        Atoms& A2 = mol.atoms[bond.a2];

        for (size_t i = 0; i < mol.no_mol; i++){
            Vec3 dp  = A1.p[i] - A2.p[i];
            Vec3 f   = -bond.K*(1 - bond.L0/dp.mag())*dp;
            A1.f[i] += f;
            A2.f[i] -= f;
            accumulated_forces_bond += f.mag();
        }
    }
}

// Angle forces (SoA version)
// changed: same idea as bonds, but with 3 atom blocks A1/A2/A3 and index i
// reason: identical angle topology per molecule; SoA gives tight loops over i
void UpdateAngleForces(System& sys){
    Molecules& mol = sys.molecules;

    for (Angle& angle : mol.angles){
        Atoms& A1 = mol.atoms[angle.a1];
        Atoms& A2 = mol.atoms[angle.a2];
        Atoms& A3 = mol.atoms[angle.a3];

        for (size_t i = 0; i < mol.no_mol; i++){
            Vec3 d21 = A2.p[i] - A1.p[i];
            Vec3 d23 = A2.p[i] - A3.p[i];

            double norm_d21 = d21.mag();
            double norm_d23 = d23.mag();
            double phi = acos(dot(d21, d23) / (norm_d21*norm_d23));

            Vec3 c21_23 = cross(d21, d23);
            Vec3 Ta = cross(d21, c21_23);
            Ta /= Ta.mag();

            Vec3 Tc = cross(c21_23, d23);
            Tc /= Tc.mag();

            Vec3 f1 = Ta*(angle.K*(phi-angle.Phi0)/norm_d21);
            Vec3 f3 = Tc*(angle.K*(phi-angle.Phi0)/norm_d23);

            A1.f[i] += f1;
            A2.f[i] -= (f1 + f3);
            A3.f[i] += f3;

            accumulated_forces_angle += f1.mag() + f3.mag();
        }
    }
}

// Non-bonded forces (SoA version)
// changed: loop over atom-types a,b instead of per-molecule atom vectors; use p[i] and p[j]
// reason: there is no "molecule.atoms" container anymore; atom blocks store all positions/forces in arrays
void UpdateNonBondedForces(System& sys){
    Molecules& mol = sys.molecules;
    double KC = 80*0.7;

    for (size_t i = 0; i < mol.no_mol; i++)
    for (auto& j : mol.neighbours[i])
    for (size_t a = 0; a < mol.atoms.size(); a++)
    for (size_t b = 0; b < mol.atoms.size(); b++){
        Atoms& atom1 = mol.atoms[a];
        Atoms& atom2 = mol.atoms[b];

        double ep = sqrt(atom1.ep*atom2.ep);
        double sigma2 = pow(0.5*(atom1.sigma+atom2.sigma),2);
        double q = KC*atom1.charge * atom2.charge;

        Vec3 dp = atom1.p[i] - atom2.p[j];   // SoA indexing: molecule i vs neighbour molecule j
        double r2 = dp.mag2();
        double r  = sqrt(r2);

        double sir = sigma2/r2;
        double sir3 = sir*sir*sir;

        Vec3 f = (ep*(12*sir3*sir3-6*sir3)*sir + q/(r*r2))*dp;

        atom1.f[i] += f;   // force on molecule i, atom-type a
        atom2.f[j] -= f;   // equal/opposite on molecule j, atom-type b
        accumulated_forces_non_bond += f.mag();
    }
}

// Integrator (SoA version)
// changed: loop atom blocks then loop i, updating v[i], p[i], resetting f[i]
// reason: makes the inner loop contiguous in memory (good cache / vectorization)
void UpdateKDK(System &sys, Sim_Configuration &sc){
    Molecules& mol = sys.molecules;

    for (auto& atom : mol.atoms)
    for (size_t i = 0; i < mol.no_mol; i++){
        atom.v[i] += sc.dt/atom.mass*atom.f[i];
        atom.f[i]  = {0,0,0};
        atom.p[i] += sc.dt*atom.v[i];
    }

    sys.time += sc.dt;
}

// Setup water system (SoA version)
// changed: allocate 3 Atoms blocks sized N_molecules, then fill mols.atoms[type].p[i]
// reason: SoA expects "all O", "all H1", "all H2" stored contiguously instead of per-molecule objects
System MakeWater(size_t N_molecules){
    const double L0 = 0.09584;
    const double angle = 104.45*deg2rad;

    Atoms O(16, 0.65,    0.31,  -0.82, "O", N_molecules);
    Atoms H1( 1, 0.18828, 0.238, 0.41, "H", N_molecules);
    Atoms H2( 1, 0.18828, 0.238, 0.41, "H", N_molecules);

    std::vector<Bond> waterbonds = {
        { .K = 20000, .L0 = L0, .a1 = 0, .a2 = 1},
        { .K = 20000, .L0 = L0, .a1 = 0, .a2 = 2}
    };

    std::vector<Angle> waterangle = {
        { .K = 1000, .Phi0 = angle, .a1 = 1, .a2 = 0, .a3 = 2 }
    };

    Molecules mols({O, H1, H2}, waterbonds, waterangle, N_molecules);

    double phi = acos(-1) * (sqrt(5.) - 1.);
    double radius = sqrt(N_molecules)*0.15;

    for (size_t i = 0; i < N_molecules; i++){
        double y = 1 - (i / (N_molecules - 1.));
        double r = sqrt(1 - y * y);
        double theta = phi * i;

        double x = cos(theta) * r;
        double z = sin(theta) * r;

        Vec3 P0{x*radius, y*radius, z*radius};

        mols.atoms[0].p[i] = {P0.x, P0.y, P0.z};
        mols.atoms[1].p[i] = {P0.x+L0*sin(angle/2), P0.y+L0*cos(angle/2), P0.z};
        mols.atoms[2].p[i] = {P0.x-L0*sin(angle/2), P0.y+L0*cos(angle/2), P0.z};
    }

    return System(std::move(mols));
}

// Output (SoA version)
// changed: loop molecule index i, then atom block a, write atom.p[i]
// reason: atoms are stored by type; molecule i is the i-th entry in each block
void WriteOutput(System& sys, std::ofstream& file){
    Molecules& mol = sys.molecules;

    for (size_t i = 0; i < mol.no_mol; i++)
    for (size_t a = 0; a < mol.atoms.size(); a++){
        auto& atom = mol.atoms[a];
        file << sys.time << " " << atom.name << " "
             << atom.p[i].x << " "
             << atom.p[i].y << " "
             << atom.p[i].z << '\n';
    }
}

//======================================================================================================
//======================== Main function ===============================================================
//======================================================================================================
int main(int argc, char* argv[]){    
    Sim_Configuration sc({argv, argv+argc}); // Load the system configuration from command line data
    
    System sys  = MakeWater(sc.no_mol);   // SoA system with sc.no_mol molecules
    std::ofstream file(sc.filename); // open file

    WriteOutput(sys, file);    // writing the initial configuration in the trajectory file
    
    auto tstart = std::chrono::high_resolution_clock::now(); // start time (nano-seconds)

    for (size_t step = 0; step < sc.steps; step++){
        // BuildNeighborList every 100th step
        if (step % 100 == 0)
            BuildNeighborList(sys);

        // Always evolve the system
        UpdateBondForces(sys);
        UpdateAngleForces(sys);
        UpdateNonBondedForces(sys);
        UpdateKDK(sys, sc);

        // Write output every data_period steps
        if (step % sc.data_period == 0)
        {
            WriteOutput(sys, file);
        }
    }

    auto tend = std::chrono::high_resolution_clock::now(); // end time (nano-seconds)

    std::cout <<  "Accumulated forces Bonds   : "  << std::setw(9) << std::setprecision(5) 
              << accumulated_forces_bond << "\n";
    std::cout <<  "Accumulated forces Angles  : "  << std::setw(9) << std::setprecision(5)
              << accumulated_forces_angle << "\n";
    std::cout <<  "Accumulated forces Non-bond: "  << std::setw(9) << std::setprecision(5)
              << accumulated_forces_non_bond << "\n";
    std::cout << "Elapsed total time:       " << std::fixed << std::setw(9) << std::setprecision(4)
              << (tend - tstart).count()*1e-9 << "\n";
}
