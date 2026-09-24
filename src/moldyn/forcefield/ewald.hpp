#pragma once

#include "moldyn/core/system.hpp"

namespace moldyn {

// SI physical constants used only to derive the Coulomb prefactor below.
// 2010 CODATA values, matching the constants NIST's own SPC/E reference
// page lists (docs/superpowers/plans/2026-09-24-phase2-ewald-water.md).
constexpr double kElementaryChargeCoulombs = 1.602176565e-19;       // C
constexpr double kVacuumPermittivityFaradPerMeter = 8.854187817e-12;  // F/m
constexpr double kBoltzmannJoulesPerKelvin = 1.3806488e-23;         // J/K
constexpr double kAngstromsPerMeter = 1e10;

// e^2 / (4 pi eps0 kB), in Kelvin*Angstrom per e^2: multiplying by
// q_i*q_j/r (r in Angstrom, q in e) gives an energy in Kelvin, matching
// NIST's E/kB reporting convention. Derived from the named constants above
// rather than pasted; should equal 167100.956632 K*A/e^2.
double coulombPrefactorKelvinAngstrom();

// The four separately retrievable Ewald energy components, in Kelvin
// (E/kB). Kept as four numbers rather than folded into one, because NIST
// publishes the breakdown and a discrepancy needs to localise to a term.
struct EwaldTerms {
    double real = 0.0;         // intermolecular pairs, r < cutoff, erfc(alpha r)/r
    double reciprocal = 0.0;   // k-space sum
    double self = 0.0;         // -C*(alpha/sqrt(pi)) * sum q_i^2
    double intra = 0.0;        // intramolecular correction (minimum-imaged)

    double total() const { return real + reciprocal + self + intra; }
};

// Ewald summation for point charges under tin-foil (conducting) boundary
// conditions, following the exact conventions verified in
// tools/spce_reference.py and docs/superpowers/plans/2026-09-24-phase2-ewald-water.md:
//
//   alpha  = 5.6 / min(box lengths)     -- recomputed per system, since it
//                                           depends on box size
//   k-space: k = 2*pi*n / L, n ranges over -kmax..kmax per axis (n != 0),
//            included only when |n|^2 < kmax^2 + 2 -- for kmax = 5 this is
//            NIST's stated |n|^2 < 27, the specific convention this Task
//            targets; it is not exercised here for any other kmax.
//   real space: brute-force i<j pairs, minimum image, intermolecular only
//               (same moleculeId is skipped -- handled by the intra term)
//   intramolecular correction: minimum-imaged distances within each
//               molecule (System::moleculeId groups them) -- required
//               because some NIST sample configurations store molecules
//               split across the periodic boundary.
class Ewald {
public:
    Ewald(double cutoff, int kmax);

    // Energy only; does not touch forces.
    EwaldTerms computeEnergies(const System& sys) const;

    // Zeroes sys's forces, then accumulates the real-space, reciprocal-space
    // and intramolecular-correction forces (the self term does not depend
    // on position, so it contributes none). Returns the same breakdown as
    // computeEnergies().
    EwaldTerms computeForces(System& sys) const;

    double cutoff() const;
    int kmax() const;

private:
    template <bool AccumulateForces>
    EwaldTerms compute(const System& sys, System* forceSink) const;

    double cutoff_;
    int kmax_;
};

}  // namespace moldyn
