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

    // Analytic long-range (tail) correction to the energy, beyond whatever
    // cutoff the implementation uses. Zero for a force field with no
    // truncated dispersion term, or for a scheme (e.g. LennardJones's
    // LinearForceShift) that shifts the force to zero at the cutoff --
    // stated explicitly as 0.0 rather than left for the caller to assume.
    virtual double longRangeCorrection(const System& sys) const = 0;

    // Analytic long-range (tail) correction to the pressure -- already
    // folded into a pressure (i.e. already divided by volume): add it
    // directly to a pressure computed as rho*T + virial/(3*V), do not
    // divide by volume again. Same zero cases as longRangeCorrection().
    virtual double longRangeCorrectionPressure(const System& sys) const = 0;
};

}  // namespace moldyn
