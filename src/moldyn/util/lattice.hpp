#pragma once

#include <cstdint>

#include "moldyn/core/system.hpp"

namespace moldyn {

// Builds N atoms on a face-centred cubic lattice filling a cubic box at the
// given number density, with velocities drawn from a Maxwell-Boltzmann
// distribution at `temperature` (kB = 1) and centre-of-mass motion removed.
// `cells` is the number of fcc unit cells per side, so N = 4 * cells^3.
// All atoms are given mass 1, matching the project's reduced-unit
// convention (see io/nist_config.hpp).
System fccLattice(int cells, double density, double temperature, uint64_t seed);

}  // namespace moldyn
