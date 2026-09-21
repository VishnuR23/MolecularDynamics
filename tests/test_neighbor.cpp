#include "doctest/doctest.h"
#include "moldyn/neighbor/cell_list.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/io/nist_config.hpp"

using namespace moldyn;

TEST_CASE("cell list reproduces brute-force forces to machine precision") {
    // config 1 is N=800 in L=10 with rc=3 -> 3 cells per side, the tight case
    auto cfg = readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic1.txt");
    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);

    auto brute = lj.computeForces(cfg.system);
    std::vector<Vec3> reference;
    for (std::size_t i = 0; i < cfg.system.size(); ++i)
        reference.push_back(cfg.system.force(i));

    CellList cells(cfg.system.box(), 3.0);
    REQUIRE(cells.usable());
    cells.build(cfg.system);
    auto fast = lj.computeForces(cfg.system, cells);

    CHECK(fast.energy == doctest::Approx(brute.energy).epsilon(1e-12));
    CHECK(fast.virial == doctest::Approx(brute.virial).epsilon(1e-12));
    for (std::size_t i = 0; i < cfg.system.size(); ++i) {
        CAPTURE(i);
        CHECK(cfg.system.force(i).x == doctest::Approx(reference[i].x).epsilon(1e-12));
        CHECK(cfg.system.force(i).y == doctest::Approx(reference[i].y).epsilon(1e-12));
        CHECK(cfg.system.force(i).z == doctest::Approx(reference[i].z).epsilon(1e-12));
    }
}

TEST_CASE("cell list still matches NIST after the optimisation") {
    auto cfg = readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic1.txt");
    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);
    CellList cells(cfg.system.box(), 3.0);
    cells.build(cfg.system);
    auto ev = lj.computeForces(cfg.system, cells);
    CHECK(ev.energy == doctest::Approx(-4.3515e3).epsilon(5e-5));
    CHECK(ev.virial == doctest::Approx(-5.6867e2).epsilon(5e-5));
}

TEST_CASE("cell list reports itself unusable for a too-small box") {
    CellList cells(Box(Vec3{5.0, 5.0, 5.0}), 3.0);   // only 1 cell per side
    CHECK_FALSE(cells.usable());
}
