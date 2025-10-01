#pragma once
#include <vector>
#include <string>

namespace mdn {

struct FeatureSet {
    std::vector<double> distances; // flattened upper-triangular distances
    std::vector<double> hbond_scores; // simplistic proxy
    std::vector<double> pocket_center; // x,y,z
};

// Parse minimal PDB/SDF (toy). In practice, replace with robust parser or RDKit in Python.
FeatureSet compute_features(const std::vector<double>& coords_xyz, int n_atoms, double hbond_cutoff_angstrom=3.5);

// TODO: Hooks for MD steps (OpenMM/GROMACS). Stubs for now.
void run_mini_md(const std::string& pdb_path, int steps, double tempK);

}
