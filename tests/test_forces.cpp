#include "doctest/doctest.h"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/io/nist_config.hpp"
#include <string>

using moldyn::LennardJones; using moldyn::Truncation; using moldyn::Vec3;

TEST_CASE("analytic force equals minus the gradient of the energy") {
    const std::string path = std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic4.txt";   // N = 30, small and fast

    for (Truncation trunc : {Truncation::Truncated, Truncation::LinearForceShift}) {
        auto cfg = moldyn::readNistConfig(path);
        LennardJones lj(1.0, 1.0, 3.0, trunc);
        lj.computeForces(cfg.system);

        const double h = 1e-5;
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
            CHECK(analytic.x == doctest::Approx(numeric[0]).epsilon(1e-6));
            CHECK(analytic.y == doctest::Approx(numeric[1]).epsilon(1e-6));
            CHECK(analytic.z == doctest::Approx(numeric[2]).epsilon(1e-6));
        }
    }
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
