#pragma once

#include <string>
#include <vector>

#include "moldyn/core/system.hpp"

namespace moldyn {

// Reads the NIST Standard Reference Simulation configuration format:
//   line 1: Lx Ly Lz
//   line 2: N
//   lines 3..N+2: index x y z [type]
// All atoms are given mass 1 and zero velocity. The optional 5th column
// (atom type, present in the SPC/E files) is returned in `types`; for the
// Lennard-Jones files it is empty.
//
// Line 2 is not trusted as an atom count: for the SPC/E files it gives the
// molecule count, while the file holds three atom lines per molecule. The
// reader reads coordinate lines until end of file instead.
struct NistConfig {
    System system;
    std::vector<std::string> types;
};

// Throws std::runtime_error on a missing file or a malformed line.
NistConfig readNistConfig(const std::string& path);

}  // namespace moldyn
