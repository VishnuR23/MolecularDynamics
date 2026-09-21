#pragma once

#include "moldyn/core/system.hpp"

namespace moldyn {

// Common result of a pairwise energy/virial evaluation over a System.
struct EnergyVirial {
    double energy = 0.0;  // pair sum, no tail correction
    double virial = 0.0;  // raw sum of w(r), NOT divided by 3
};

// Minimal interface shared by force fields that can report a pairwise
// energy/virial and an analytic long-range correction. Phase 1 has exactly
// one implementation (LennardJones); this exists so later phases can add
// force fields behind the same interface without touching call sites.
class ForceField {
public:
    virtual ~ForceField() = default;

    virtual EnergyVirial computeEnergyVirial(const System& sys) const = 0;
    virtual double longRangeCorrection(const System& sys) const = 0;
};

}  // namespace moldyn
