#pragma once

#include <cstddef>

#include "moldyn/core/system.hpp"
#include "moldyn/io/nist_config.hpp"

namespace moldyn {

// SPC/E rigid three-site water model parameters, matching NIST's Ewald/SPC-E
// reference calculation exactly (see
// docs/superpowers/plans/2026-09-24-phase2-ewald-water.md). rOH and
// angleHOH describe the rigid geometry for later tasks (SETTLE); they are
// not used by the energy calculation here, which takes atom positions as
// given.
struct SpceParameters {
    double qH = 0.42380;             // e; qO = -2 * qH
    double sigmaO = 3.16555789;      // Angstrom, oxygen LJ sigma
    double epsilonO = 78.19743111;   // Kelvin (epsilon/kB), oxygen LJ epsilon
    double rOH = 1.0;                // Angstrom, rigid O-H bond length
    double angleHOH = 109.47;        // degrees, rigid H-O-H angle
    double massOxygen = 15.9994;     // amu
    double massHydrogen = 1.008;     // amu
};

// A parsed NIST SPC/E configuration turned into a charged, molecule-tagged
// System (`full`) plus an oxygen-only System in the same box (`oxygens`),
// since dispersion is defined on oxygen sites only.
struct SpceSystem {
    System full;
    System oxygens;
    std::size_t moleculeCount = 0;
};

// Builds `full` and `oxygens` from a config whose `types` are "O"/"H",
// stored as molecules of three consecutive atoms (O, H, H repeating, as
// readNistConfig documents for the SPC/E files). Molecule id = atom index
// / 3. Throws std::invalid_argument if the atom count isn't a multiple of
// 3, `types` is empty, or a type other than "O"/"H" is found.
SpceSystem buildSpceSystem(const NistConfig& cfg, const SpceParameters& params = SpceParameters());

}  // namespace moldyn
