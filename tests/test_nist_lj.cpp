#include "doctest/doctest.h"
#include "moldyn/core/box.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/io/nist_config.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>

using moldyn::Box;
using moldyn::LennardJones;
using moldyn::System;
using moldyn::Truncation;
using moldyn::Vec3;

namespace {
struct Ref {
    int config; Truncation trunc; double rc;
    double upair; double wpair; double ulrc;
};

// NIST Standard Reference Simulation Website, Lennard-Jones fluid,
// cuboid cell. Values are reported to 5 significant figures.
const Ref kRefs[] = {
    {1, Truncation::Truncated,        3.0, -4.3515e3, -5.6867e2, -1.9849e2},
    {1, Truncation::Truncated,        4.0, -4.4675e3, -1.2639e3, -8.3769e1},
    {1, Truncation::LinearForceShift, 3.0, -3.8709e3,  3.1754e2,  0.0     },
    {2, Truncation::Truncated,        3.0, -6.9000e2, -5.6846e2, -2.4230e1},
    {2, Truncation::Truncated,        4.0, -7.0460e2, -6.5599e2, -1.0226e1},
    {2, Truncation::LinearForceShift, 3.0, -6.2012e2, -4.4533e2,  0.0     },
    {3, Truncation::Truncated,        3.0, -1.1467e3, -1.1649e3, -4.9622e1},
    {3, Truncation::Truncated,        4.0, -1.1754e3, -1.3371e3, -2.0942e1},
    {3, Truncation::LinearForceShift, 3.0, -1.0210e3, -9.3578e2,  0.0     },
    {4, Truncation::Truncated,        3.0, -1.6790e1, -4.6249e1, -5.4517e-1},
    {4, Truncation::Truncated,        4.0, -1.7060e1, -4.7869e1, -2.3008e-1},
    {4, Truncation::LinearForceShift, 3.0, -1.5001e1, -4.3096e1,  0.0     },
};

// Builds a System with exactly n atoms in a cubic box whose volume gives
// number density rho = n/V. longRangeCorrection[Pressure] depend only on
// N and V (not on where the atoms actually sit), so placing every atom at
// the origin is fine -- these tests never touch positions or forces.
System makeSystemAtDensity(std::size_t n, double rho) {
    const double volume = static_cast<double>(n) / rho;
    const double side = std::cbrt(volume);
    System sys(Box(Vec3{side, side, side}));
    sys.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        sys.addAtom(Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 0.0}, 1.0);
    }
    return sys;
}

// Reduced-unit state points (sigma = epsilon = 1) with independently
// hand-computed tail corrections, checked against
// P_LRC = (16/3)*pi*rho^2*eps*sig^3*[(2/3)*(sig/rc)^9 - (sig/rc)^3]
// alongside U_LRC/N from the existing longRangeCorrection().
struct PressureTailPoint {
    double rho;
    double rc;
    double uTailPerN;
    double pTail;
};

const PressureTailPoint kPressureTailPoints[] = {
    {0.8442, 2.5, -0.4520, -0.7621},
    {0.8442, 3.0, -0.2618, -0.4419},
    {0.7000, 2.5, -0.3748, -0.5240},
    {0.5000, 2.5, -0.2677, -0.2674},
};
// Half a unit in the last significant figure NIST printed.
//
// Every value in kRefs is published to exactly 5 significant figures, so
// the true value lies within half a unit of the last printed digit --
// e.g. -4.3515e3 stands for a true value in [-4351.55, -4351.45], a
// half-width of 0.05. That per-value half-width is exactly how precisely
// the reference is known, so it is both the tightest bound this
// comparison can honestly assert and the loosest one worth asserting.
//
// This replaces a flat 5e-5 relative tolerance, which was up to 9.4x
// looser than NIST's actual precision -- worst at W = -935.78, whose
// relative half-width is 0.005/935.78 = 5.34e-6 -- so a regression that
// broke the fifth published digit would have passed it. Do not loosen this back to a flat
// number. Every value here agrees with NIST to every digit NIST prints,
// and that is the claim under test.
double roundingHalfWidth(double referenceValue) {
    const double magnitude = std::abs(referenceValue);
    if (magnitude == 0.0) {
        return 0.0;
    }
    const double decadeExponent = std::floor(std::log10(magnitude));
    return 0.5 * std::pow(10.0, decadeExponent - 4.0);
}

// The LinearForceShift rows' U_LRC is an exact structural zero (that
// scheme shifts the potential to vanish at rc), not a rounded
// measurement, so there is no last printed digit to halve.
constexpr double kZeroLrcTolerance = 1e-12;
}  // namespace

TEST_CASE("roundingHalfWidth is half a unit in NIST's fifth significant figure") {
    // The helper the tolerance below rests on, checked against
    // hand-computed values spanning the four decades kRefs covers.
    CHECK(roundingHalfWidth(-4.3515e3) == doctest::Approx(0.05));
    CHECK(roundingHalfWidth(-5.6867e2) == doctest::Approx(0.005));
    CHECK(roundingHalfWidth(-1.6790e1) == doctest::Approx(0.0005));
    CHECK(roundingHalfWidth(-5.4517e-1) == doctest::Approx(5e-6));
    CHECK(roundingHalfWidth(0.0) == 0.0);
}

TEST_CASE("Lennard-Jones energy and virial match every digit NIST publishes") {
    double worstRatio = 0.0;

    for (const Ref& r : kRefs) {
        CAPTURE(r.config);
        CAPTURE(r.rc);
        const std::string path = std::string(MOLDYN_DATA_DIR) +
            "/nist/lj/lj_sample_config_periodic" + std::to_string(r.config) + ".txt";
        auto cfg = moldyn::readNistConfig(path);

        LennardJones lj(1.0, 1.0, r.rc, r.trunc);
        auto ev = lj.computeEnergyVirial(cfg.system);

        const double uTol = roundingHalfWidth(r.upair);
        const double wTol = roundingHalfWidth(r.wpair);
        const double uErr = std::abs(ev.energy - r.upair);
        const double wErr = std::abs(ev.virial - r.wpair);
        CAPTURE(uErr);
        CAPTURE(uTol);
        CAPTURE(wErr);
        CAPTURE(wTol);
        CHECK(uErr <= uTol);
        CHECK(wErr <= wTol);
        worstRatio = std::max({worstRatio, uErr / uTol, wErr / wTol});

        const double lrc = lj.longRangeCorrection(cfg.system);
        if (r.ulrc == 0.0) {
            CHECK(std::abs(lrc) <= kZeroLrcTolerance);
        } else {
            const double lrcTol = roundingHalfWidth(r.ulrc);
            const double lrcErr = std::abs(lrc - r.ulrc);
            CAPTURE(lrcErr);
            CAPTURE(lrcTol);
            CHECK(lrcErr <= lrcTol);
            worstRatio = std::max(worstRatio, lrcErr / lrcTol);
        }
    }

    // Reported, not asserted at a threshold: this is the headline number
    // ("our worst value uses this fraction of NIST's own precision"), and
    // the per-value CHECKs above are what actually gates it.
    MESSAGE("worst case: " << worstRatio << " of the per-value rounding half-width");
}

TEST_CASE("Lennard-Jones pressure long-range correction matches hand-computed values") {
    // N is arbitrary here: U_LRC is extensive (checked per-atom below) and
    // P_LRC is intensive, so any atom count at the target density works.
    const std::size_t n = 500;

    for (const PressureTailPoint& p : kPressureTailPoints) {
        CAPTURE(p.rho);
        CAPTURE(p.rc);

        System sys = makeSystemAtDensity(n, p.rho);
        LennardJones lj(1.0, 1.0, p.rc, Truncation::Truncated);

        const double uLrcPerN = lj.longRangeCorrection(sys) / static_cast<double>(n);
        CHECK(uLrcPerN == doctest::Approx(p.uTailPerN).epsilon(2e-3));

        const double pLrc = lj.longRangeCorrectionPressure(sys);
        CHECK(pLrc == doctest::Approx(p.pTail).epsilon(2e-3));
    }
}

TEST_CASE("Lennard-Jones pressure long-range correction is exactly zero for LinearForceShift") {
    System sys = makeSystemAtDensity(500, 0.8442);
    LennardJones lj(1.0, 1.0, 2.5, Truncation::LinearForceShift);

    CHECK(lj.longRangeCorrectionPressure(sys) == 0.0);
}

TEST_CASE("Pressure LRC is intensive; energy LRC is extensive") {
    // Fixed density (rho = 0.5, exactly representable), N and V doubled
    // together. rho = n/volume is bit-identical for both systems, so this
    // isolates whether the formulas themselves carry an extra factor of N.
    const double rho = 0.5;
    System small = makeSystemAtDensity(500, rho);
    System big = makeSystemAtDensity(1000, rho);

    LennardJones lj(1.0, 1.0, 2.5, Truncation::Truncated);

    const double pSmall = lj.longRangeCorrectionPressure(small);
    const double pBig = lj.longRangeCorrectionPressure(big);
    CHECK(pBig == doctest::Approx(pSmall).epsilon(1e-12));

    const double uSmall = lj.longRangeCorrection(small);
    const double uBig = lj.longRangeCorrection(big);
    CHECK(uBig == doctest::Approx(2.0 * uSmall).epsilon(1e-12));
}

TEST_CASE("Pressure LRC shrinks with cutoff and the (sigma/rc)^3 term dominates at large rc") {
    System sys = makeSystemAtDensity(500, 0.8442);

    const double p20 = LennardJones(1.0, 1.0, 2.0, Truncation::Truncated).longRangeCorrectionPressure(sys);
    const double p25 = LennardJones(1.0, 1.0, 2.5, Truncation::Truncated).longRangeCorrectionPressure(sys);
    const double p30 = LennardJones(1.0, 1.0, 3.0, Truncation::Truncated).longRangeCorrectionPressure(sys);
    const double p60 = LennardJones(1.0, 1.0, 6.0, Truncation::Truncated).longRangeCorrectionPressure(sys);

    // All negative (attractive tail dominates), shrinking in magnitude as
    // rc grows.
    CHECK(p20 < p25);
    CHECK(p25 < p30);
    CHECK(p30 < p60);
    CHECK(p60 < 0.0);

    // At large rc the (sig/rc)^9 term is negligible, so P_LRC should be
    // within a fraction of a percent of the (sig/rc)^3-only term.
    constexpr double kPi = 3.14159265358979323846;
    const double rho = 0.8442;
    const double sigOverRc3 = std::pow(1.0 / 6.0, 3);
    const double leadingTermOnly = (16.0 / 3.0) * kPi * rho * rho * 1.0 * 1.0 * (-sigOverRc3);
    CHECK(p60 == doctest::Approx(leadingTermOnly).epsilon(1e-3));
}
