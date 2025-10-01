#include "backend.hpp"
#include "utils.hpp"
#include <stdexcept>
#include <cmath>

namespace mdn {

FeatureSet compute_features(const std::vector<double>& coords_xyz, int n_atoms, double hbond_cutoff_angstrom) {
    if ((int)coords_xyz.size() != n_atoms*3) throw std::runtime_error("coords size mismatch");
    FeatureSet f;
    // distances (upper triangular)
    for (int i=0;i<n_atoms;i++){
        for (int j=i+1;j<n_atoms;j++){
            double dx = coords_xyz[3*i+0]-coords_xyz[3*j+0];
            double dy = coords_xyz[3*i+1]-coords_xyz[3*j+1];
            double dz = coords_xyz[3*i+2]-coords_xyz[3*j+2];
            double d = std::sqrt(dx*dx+dy*dy+dz*dz);
            f.distances.push_back(d);
            // silly "hbond score" proxy: close and within cutoff
            if (d < hbond_cutoff_angstrom) f.hbond_scores.push_back(1.0 - d/hbond_cutoff_angstrom);
        }
    }
    // pocket center: mean of coords
    double cx=0, cy=0, cz=0;
    for (int i=0;i<n_atoms;i++){
        cx += coords_xyz[3*i+0];
        cy += coords_xyz[3*i+1];
        cz += coords_xyz[3*i+2];
    }
    cx/=n_atoms; cy/=n_atoms; cz/=n_atoms;
    f.pocket_center = {cx, cy, cz};
    return f;
}

void run_mini_md(const std::string& pdb_path, int steps, double tempK) {
    // Stub: integrate OpenMM here if available. For the portfolio, leaving a clear hook is enough.
    (void)pdb_path; (void)steps; (void)tempK;
    // no-op
}

}
