#include "doctest/doctest.h"
#include "moldyn/core/box.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/forcefield/ewald.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/forcefield/spce.hpp"
#include "moldyn/io/nist_config.hpp"

#include <cmath>
#include <cstddef>
#include <string>

using moldyn::Box;
using moldyn::Ewald;
using moldyn::EwaldTerms;
using moldyn::LennardJones;
using moldyn::SpceSystem;
using moldyn::System;
using moldyn::Truncation;
using moldyn::Vec3;

namespace {

// Same idea as tests/test_nist_lj.cpp's roundingHalfWidth (half a unit in
// NIST's last published significant figure -- that's exactly how precisely
// each published value is known), reused rather than reimplemented for the
// concept, but NOT byte-identical: that helper is private to its
// translation unit (the header exposes no such function, so there is
// nothing to literally #include), and it is hardcoded to the LJ table's 5
// significant figures (exponent offset -4). The SPC/E table publishes 6
// (e.g. 9.95387e4, -2.84469e6 -- six digits each), so the offset here is
// -5: a value's last digit sits at 10^(decadeExponent - (sigFigs - 1)),
// and half a unit of that digit is 0.5*10^(decadeExponent - (sigFigs-1)).
// Getting this wrong (i.e. copying -4 unchanged) makes every bound here
// 10x too loose -- caught by checking the derived Etotal bounds below
// against the plan's table (2.16e-05, 1.46e-05, 9.04e-06, 3.29e-05), which
// only come out right with -5.
double roundingHalfWidth(double referenceValue) {
    const double magnitude = std::abs(referenceValue);
    if (magnitude == 0.0) {
        return 0.0;
    }
    const double decadeExponent = std::floor(std::log10(magnitude));
    return 0.5 * std::pow(10.0, decadeExponent - 5.0);
}

// NIST Standard Reference Simulation Website, SPC/E water, Ewald summation,
// cuboid cell, tin-foil boundary conditions. E/kB in Kelvin, published to 6
// significant figures.
//   https://www.nist.gov/mml/csd/informatics_research_group/spce_water
struct Ref {
    int config;
    double edisp;
    double elrc;
    double ereal;
    double efourier;
    double eself;
    double eintra;
    double etotal;
};

const Ref kRefs[] = {
    {1, 9.95387e4, -8.23715e2, -5.58889e5, 6.27009e3, -2.84469e6, 2.80999e6, -4.88604e5},
    {2, 1.93712e5, -3.29486e3, -1.19295e6, 6.03495e3, -5.68938e6, 5.61998e6, -1.06590e6},
    {3, 3.54344e5, -7.41343e3, -1.96297e6, 5.24461e3, -8.53407e6, 8.42998e6, -1.71488e6},
    {4, 4.48593e5, -1.37286e4, -3.57226e6, 7.58785e3, -1.42235e7, 1.41483e7, -3.20501e6},
};

std::string configPath(int config) {
    return std::string(MOLDYN_DATA_DIR) + "/nist/spce/spce_sample_config_periodic" +
           std::to_string(config) + ".txt";
}

}  // namespace

TEST_CASE("roundingHalfWidth is half a unit in NIST's sixth significant figure") {
    CHECK(roundingHalfWidth(9.95387e4) == doctest::Approx(0.05));
    CHECK(roundingHalfWidth(-2.84469e6) == doctest::Approx(5.0));
    CHECK(roundingHalfWidth(-8.23715e2) == doctest::Approx(0.0005));
    CHECK(roundingHalfWidth(0.0) == 0.0);
}

TEST_CASE("Coulomb prefactor matches NIST's own physical constants") {
    // e^2 / (4 pi eps0 kB), in K*Angstrom per e^2 -- derived in code from
    // named SI constants, not pasted, then checked against the number NIST
    // publishes (docs/superpowers/plans/2026-09-24-phase2-ewald-water.md).
    CHECK(moldyn::coulombPrefactorKelvinAngstrom() ==
          doctest::Approx(167100.956632).epsilon(1e-9));
}

TEST_CASE("SPC/E Ewald energy components match every digit NIST publishes, term by term") {
    double worstRatio = 0.0;

    for (const Ref& r : kRefs) {
        CAPTURE(r.config);
        auto nist = moldyn::readNistConfig(configPath(r.config));
        SpceSystem spce = moldyn::buildSpceSystem(nist);

        LennardJones oxygenLj(3.16555789, 78.19743111, 10.0, Truncation::Truncated);
        const double edisp = oxygenLj.computeEnergyVirial(spce.oxygens).energy;
        const double elrc = oxygenLj.longRangeCorrection(spce.oxygens);

        Ewald ewald(10.0, 5);
        const EwaldTerms terms = ewald.computeEnergies(spce.full);

        struct Component {
            const char* name;
            double got;
            double ref;
        };
        const Component components[] = {
            {"Edisp", edisp, r.edisp},
            {"ELRC", elrc, r.elrc},
            {"Ereal", terms.real, r.ereal},
            {"Efourier", terms.reciprocal, r.efourier},
            {"Eself", terms.self, r.eself},
            {"Eintra", terms.intra, r.eintra},
        };

        double sumOfHalfWidths = 0.0;
        for (const Component& c : components) {
            CAPTURE(c.name);
            const double tol = roundingHalfWidth(c.ref);
            const double err = std::abs(c.got - c.ref);
            CAPTURE(c.got);
            CAPTURE(c.ref);
            CAPTURE(err);
            CAPTURE(tol);
            CHECK(err <= tol);
            worstRatio = std::max(worstRatio, err / tol);
            sumOfHalfWidths += tol;
        }

        // Etotal's bound is derived, not chosen: the sum of the six
        // components' half-widths, divided by |Etotal| (it is a small
        // residual of large cancelling terms, so component rounding
        // amplifies). See the plan's table for the four values this
        // produces.
        const double etotal = edisp + elrc + terms.real + terms.reciprocal + terms.self + terms.intra;
        const double derivedBound = sumOfHalfWidths / std::abs(r.etotal);
        const double totalErr = std::abs(etotal - r.etotal) / std::abs(r.etotal);
        CAPTURE(etotal);
        CAPTURE(r.etotal);
        CAPTURE(derivedBound);
        CAPTURE(totalErr);
        CHECK(totalErr <= derivedBound);
        worstRatio = std::max(worstRatio, totalErr / derivedBound);
    }

    MESSAGE("worst case: " << worstRatio << " of the per-value bound");
}

TEST_CASE("Ewald forces equal minus the gradient of the Ewald energy") {
    // Configuration 1 (N=300) is small enough to keep this fast while still
    // exercising real space, reciprocal space, and the intramolecular
    // correction together.
    auto nist = moldyn::readNistConfig(configPath(1));
    SpceSystem spce = moldyn::buildSpceSystem(nist);

    Ewald ewald(10.0, 5);
    ewald.computeForces(spce.full);

    const double h = 1e-5;
    const double kGradientTolerance = 1e-7;
    double worstRelative = 0.0;

    // A handful of atoms spanning both elements and several molecules.
    for (std::size_t i : {std::size_t(0), std::size_t(1), std::size_t(2), std::size_t(150),
                           std::size_t(299)}) {
        Vec3 analytic = spce.full.force(i);
        Vec3 original = spce.full.position(i);
        double numeric[3];
        for (int d = 0; d < 3; ++d) {
            Vec3 plus = original;
            Vec3 minus = original;
            (d == 0 ? plus.x : d == 1 ? plus.y : plus.z) += h;
            (d == 0 ? minus.x : d == 1 ? minus.y : minus.z) -= h;

            spce.full.setPosition(i, plus);
            const double up = ewald.computeEnergies(spce.full).total();
            spce.full.setPosition(i, minus);
            const double um = ewald.computeEnergies(spce.full).total();
            spce.full.setPosition(i, original);

            numeric[d] = -(up - um) / (2.0 * h);
        }
        CAPTURE(i);
        CHECK(analytic.x == doctest::Approx(numeric[0]).epsilon(kGradientTolerance));
        CHECK(analytic.y == doctest::Approx(numeric[1]).epsilon(kGradientTolerance));
        CHECK(analytic.z == doctest::Approx(numeric[2]).epsilon(kGradientTolerance));

        const double components[3] = {analytic.x, analytic.y, analytic.z};
        for (int d = 0; d < 3; ++d) {
            const double scale = std::max(std::abs(components[d]), std::abs(numeric[d]));
            if (scale > 0.0) {
                worstRelative = std::max(worstRelative, std::abs(components[d] - numeric[d]) / scale);
            }
        }
    }

    MESSAGE("worst relative Ewald force/gradient disagreement: " << worstRelative);
}

TEST_CASE("Ewald forces sum to zero (Newton's third law)") {
    auto nist = moldyn::readNistConfig(configPath(2));
    SpceSystem spce = moldyn::buildSpceSystem(nist);

    Ewald ewald(10.0, 5);
    ewald.computeForces(spce.full);

    Vec3 total{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < spce.full.size(); ++i) {
        total += spce.full.force(i);
    }
    CHECK(total.x == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(total.y == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(total.z == doctest::Approx(0.0).epsilon(1e-8));
}
