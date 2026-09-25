#pragma once

namespace moldyn {
namespace units {

// The unit system the water side of this engine works in.
//
// Phase 1 used Lennard-Jones reduced units, where sigma = epsilon = mass =
// kB = 1 and nothing needed converting. SPC/E is parameterised in real units,
// so the choices have to be stated:
//
//     length   Angstrom
//     mass     atomic mass unit
//     charge   elementary charge
//     energy   Kelvin (that is, E/kB -- NIST reports the SPC/E reference
//              energies this way, so taking the same convention means the
//              reference values need no conversion at all)
//
// Those four fix the time unit. Requiring that 0.5*m*v^2 comes out in Kelvin
// with mass in amu and length in Angstrom gives
//
//     tau = Angstrom * sqrt(amu / kB) = 1.096688 ps
//
// which is not a round number, so the conversions below exist to keep call
// sites from open-coding it. A 2 fs timestep -- the usual choice for rigid
// water -- is 0.00182367 tau.
//
// Getting this wrong is silent: the dynamics stay self-consistent and only the
// physical timescale is off, so temperatures and structure still look right
// while every rate is wrong by a constant factor.

// One internal time unit, in picoseconds.
inline constexpr double kTauInPicoseconds = 1.0966875353;

// Convert a timestep to internal units.
inline constexpr double picoseconds(double ps) { return ps / kTauInPicoseconds; }
inline constexpr double femtoseconds(double fs) { return picoseconds(fs * 1.0e-3); }

// And back, for reporting.
inline constexpr double toPicoseconds(double tau) { return tau * kTauInPicoseconds; }

}  // namespace units
}  // namespace moldyn
