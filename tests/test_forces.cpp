#include "doctest/doctest.h"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/io/nist_config.hpp"
#include <algorithm>
#include <cmath>
#include <string>

using moldyn::LennardJones; using moldyn::Truncation; using moldyn::Vec3;

TEST_CASE("analytic force equals minus the gradient of the energy") {
    const std::string path = std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic4.txt";   // N = 30, small and fast

    double worstRelative = 0.0;

    for (Truncation trunc : {Truncation::Truncated, Truncation::LinearForceShift}) {
        auto cfg = moldyn::readNistConfig(path);
        LennardJones lj(1.0, 1.0, 3.0, trunc);
        lj.computeForces(cfg.system);

        // 1e-8, which is what the spec asks for ("Analytic force equals
        // central finite difference of the energy to 1e-8") and what the
        // plan predicts for this step size: central differences reach
        // about eps_machine^(2/3) accuracy, so h = 1e-5 should agree near
        // 1e-8. It was 1e-6 -- two decades looser than the method's own
        // accuracy, loose enough that a genuine sign or factor error in a
        // small force component could hide inside it.
        //
        // The worst relative disagreement actually observed across these
        // nine components is reported by the MESSAGE below, so the margin
        // is a measured number in the test output, not a claim in a
        // comment. If this ever fails, suspect the force kernel before
        // suspecting the tolerance.
        const double h = 1e-5;
        const double kGradientTolerance = 1e-8;
        for (std::size_t i : {std::size_t(0), std::size_t(7), std::size_t(29)}) {
            Vec3 analytic = cfg.system.force(i);
            Vec3 original = cfg.system.position(i);
            double numeric[3];
            for (int d = 0; d < 3; ++d) {
                Vec3 plus = original, minus = original;
                (d == 0 ? plus.x : d == 1 ? plus.y : plus.z)   += h;
                (d == 0 ? minus.x : d == 1 ? minus.y : minus.z) -= h;

                cfg.system.setPosition(i, plus);
                double up = lj.computeEnergyVirial(cfg.system).energy;
                cfg.system.setPosition(i, minus);
                double um = lj.computeEnergyVirial(cfg.system).energy;
                cfg.system.setPosition(i, original);

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
                    worstRelative =
                        std::max(worstRelative, std::abs(components[d] - numeric[d]) / scale);
                }
            }
        }
    }

    MESSAGE("worst relative force/gradient disagreement: " << worstRelative
                                                            << " (tolerance 1e-8)");
}

TEST_CASE("forces sum to zero (Newton's third law)") {
    auto cfg = moldyn::readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic4.txt");
    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);
    lj.computeForces(cfg.system);
    Vec3 total{0, 0, 0};
    for (std::size_t i = 0; i < cfg.system.size(); ++i) total += cfg.system.force(i);
    CHECK(total.x == doctest::Approx(0.0).epsilon(1e-10));
    CHECK(total.y == doctest::Approx(0.0).epsilon(1e-10));
    CHECK(total.z == doctest::Approx(0.0).epsilon(1e-10));
}

TEST_CASE("computeForces reports the same energy as computeEnergyVirial") {
    auto cfg = moldyn::readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic2.txt");
    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);
    auto a = lj.computeEnergyVirial(cfg.system);
    auto b = lj.computeForces(cfg.system);
    CHECK(a.energy == doctest::Approx(b.energy).epsilon(1e-12));
    CHECK(a.virial == doctest::Approx(b.virial).epsilon(1e-12));
}
