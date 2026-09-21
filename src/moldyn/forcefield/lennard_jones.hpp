#pragma once

#include "moldyn/core/system.hpp"
#include "moldyn/forcefield/forcefield.hpp"

namespace moldyn {

enum class Truncation { Truncated, LinearForceShift };

// Lennard-Jones pair potential, evaluated by brute-force summation over all
// i<j pairs under the minimum image convention.
//
//   u_LJ(r) = 4*eps*[ (sig/r)^12 - (sig/r)^6 ]
//   w_LJ(r) = -r * du_LJ/dr = 24*eps*[ 2*(sig/r)^12 - (sig/r)^6 ]
//
// Truncated ("LRC" in the NIST tables): the potential is cut off at rc with
// no shift; longRangeCorrection() supplies the analytic tail separately.
//
// LinearForceShift ("LFS"): both the potential and its derivative are
// shifted to zero at rc, so longRangeCorrection() is exactly zero.
class LennardJones : public ForceField {
public:
    LennardJones(double sigma, double epsilon, double cutoff, Truncation truncation);

    // Pair sum by brute force over all i<j with the minimum image convention.
    // Does not touch forces.
    EnergyVirial computeEnergyVirial(const System& sys) const override;

    // Accumulates forces into sys and returns the same energy/virial as
    // computeEnergyVirial. Calls sys.zeroForces() first.
    //   f_vec_ij = (w(r) / r^2) * r_vec_ij,  r_vec_ij = r_i - r_j (minimum-imaged)
    //   atom i gains +f_vec_ij, atom j gains -f_vec_ij
    EnergyVirial computeForces(System& sys) const;

    // Analytic long-range correction to the energy. Zero for LinearForceShift.
    double longRangeCorrection(const System& sys) const override;

    double sigma() const;
    double epsilon() const;
    double cutoff() const;
    Truncation truncation() const;

private:
    // Shared pair kernel for both computeEnergyVirial and computeForces:
    // the i<j iteration, minimum image, cutoff test, and uLJ/wLJ arithmetic
    // live here exactly once. AccumulateForces is a compile-time switch so
    // the force-accumulation branch compiles away entirely for the
    // energy-only path. forceSink is written to iff AccumulateForces is
    // true, and must be non-null in that case (it aliases sys).
    template <bool AccumulateForces>
    EnergyVirial computePairSum(const System& sys, System* forceSink) const;

    double sigma_;
    double epsilon_;
    double cutoff_;
    Truncation truncation_;

    double uAtCutoff_;   // u_LJ(rc), precomputed
    double dudrAtCutoff_;  // du_LJ/dr at rc, precomputed
};

}  // namespace moldyn
