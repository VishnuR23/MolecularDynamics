#pragma once

#include "moldyn/core/system.hpp"
#include "moldyn/forcefield/forcefield.hpp"

namespace moldyn {

class CellList;
class VerletList;

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

    // Same computation as computeForces(sys), but only visiting the pairs the
    // neighbour list offers instead of every i<j pair. Results are identical
    // to the brute-force path to machine precision -- these are a pure
    // optimisation of how pairs are found, not a different calculation.
    // The CellList overload throws std::invalid_argument when
    // cells.usable() is false -- that list cannot enumerate pairs
    // correctly, and returning a plausible-looking wrong number instead
    // of failing would be worse than useless. Use the brute-force
    // overload for boxes too small to divide.
    EnergyVirial computeForces(System& sys, const CellList& cells) const;
    EnergyVirial computeForces(System& sys, const VerletList& verletList) const;

    // Analytic long-range correction to the energy. Zero for LinearForceShift.
    double longRangeCorrection(const System& sys) const override;

    // Analytic long-range correction to the pressure -- the virial tail
    // beyond the cutoff, already folded into a pressure (i.e. already
    // divided by volume): add it directly to a pressure computed as
    // rho*T + virial/(3*V), do not divide by volume again. Zero for
    // LinearForceShift, for the same reason longRangeCorrection() is:
    // that scheme shifts the force to zero at rc, so there is no
    // uncounted tail to correct.
    //
    // Intensive (depends only on density, not on N or V separately),
    // unlike longRangeCorrection() which is extensive. Like the energy
    // correction, this assumes the radial distribution function is
    // uniform (g(r) = 1) beyond rc -- accurate for a dense, homogeneous
    // fluid far from a phase boundary, which is the regime Phase 1
    // targets; it is not a substitute for a real g(r) tail near a
    // critical point or an interface.
    double longRangeCorrectionPressure(const System& sys) const override;

    double sigma() const;
    double epsilon() const;
    double cutoff() const;
    Truncation truncation() const;

private:
    // Shared pair kernel for computeEnergyVirial and every computeForces
    // overload: the minimum image, cutoff test, and uLJ/wLJ arithmetic live
    // here exactly once. AccumulateForces is a compile-time switch so the
    // force-accumulation branch compiles away entirely for the energy-only
    // path. forceSink is written to iff AccumulateForces is true, and must
    // be non-null in that case (it aliases sys).
    //
    // PairSource supplies the candidate (i, j) pairs via a
    // `forEachPair(fn)` template method with i < j -- the brute-force i<j
    // double loop, CellList, and VerletList all satisfy this, so this is
    // the only place that ever computes uLJ/wLJ from r2.
    template <bool AccumulateForces, class PairSource>
    EnergyVirial computePairSum(const System& sys, System* forceSink,
                                 const PairSource& pairs) const;

    double sigma_;
    double epsilon_;
    double cutoff_;
    Truncation truncation_;

    double uAtCutoff_;   // u_LJ(rc), precomputed
    double dudrAtCutoff_;  // du_LJ/dr at rc, precomputed
};

}  // namespace moldyn
