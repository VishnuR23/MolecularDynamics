#pragma once

#include "moldyn/core/system.hpp"
#include "moldyn/forcefield/ewald.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/forcefield/spce.hpp"

namespace moldyn {

// The complete SPC/E interaction: Lennard-Jones between oxygens plus Ewald
// electrostatics over every charged site, accumulated onto one System.
//
// buildSpceSystem hands back the oxygens as a separate System, which is the
// right shape for reproducing NIST's energy breakdown -- dispersion is defined
// on oxygen sites alone, so it is natural to sum it over a system containing
// only those. It is the wrong shape for dynamics: forces computed on that copy
// never reach the full system, and the two drift apart the moment anything
// moves. Electrostatics without the dispersion repulsion is also catastrophic,
// since nothing then stops the charges collapsing together.
//
// So the gather and scatter live here, in one place, rather than at each call
// site where they would eventually be forgotten.
class SpceForceField {
public:
    explicit SpceForceField(const SpceParameters& params = SpceParameters(),
                            double cutoff = 10.0,
                            int kmax = 5);

    // Zeroes the forces on `sys`, then accumulates dispersion and
    // electrostatics. Returns the total potential energy in Kelvin, including
    // the analytic dispersion tail correction.
    //
    // `sys` must be the `full` system from buildSpceSystem: charged, with
    // molecules as consecutive O, H, H triples.
    double computeForces(System& sys) const;

    // The same energy without touching forces.
    double computeEnergy(const System& sys) const;

    const Ewald& ewald() const { return ewald_; }
    const LennardJones& dispersion() const { return lj_; }

private:
    // An oxygen-only System mirroring the oxygens of `sys`, for the dispersion
    // term. Positions are copied in; forces are copied back by the caller.
    System gatherOxygens(const System& sys) const;

    SpceParameters params_;
    LennardJones lj_;
    Ewald ewald_;
};

}  // namespace moldyn
