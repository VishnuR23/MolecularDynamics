#include "doctest/doctest.h"
#include "moldyn/neighbor/cell_list.hpp"
#include "moldyn/neighbor/verlet_list.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/io/nist_config.hpp"
#include <cstddef>
#include <vector>

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

namespace {
// config 1 (N=800, L=10, rc=3) leaves room for a 0.3 skin: cutoff+skin=3.3,
// floor(10/3.3)=3 cells/side, still usable.
constexpr double kSkin = 0.3;
}  // namespace

TEST_CASE("Verlet list reproduces brute-force forces to machine precision") {
    auto cfg = readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic1.txt");
    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);

    auto brute = lj.computeForces(cfg.system);
    std::vector<Vec3> reference;
    for (std::size_t i = 0; i < cfg.system.size(); ++i)
        reference.push_back(cfg.system.force(i));

    VerletList verlet(cfg.system.box(), 3.0, kSkin);
    verlet.build(cfg.system);
    auto fast = lj.computeForces(cfg.system, verlet);

    CHECK(fast.energy == doctest::Approx(brute.energy).epsilon(1e-12));
    CHECK(fast.virial == doctest::Approx(brute.virial).epsilon(1e-12));
    for (std::size_t i = 0; i < cfg.system.size(); ++i) {
        CAPTURE(i);
        CHECK(cfg.system.force(i).x == doctest::Approx(reference[i].x).epsilon(1e-12));
        CHECK(cfg.system.force(i).y == doctest::Approx(reference[i].y).epsilon(1e-12));
        CHECK(cfg.system.force(i).z == doctest::Approx(reference[i].z).epsilon(1e-12));
    }
}

TEST_CASE("a displacement under half the skin does not trigger a rebuild") {
    auto cfg = readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic1.txt");
    VerletList verlet(cfg.system.box(), 3.0, kSkin);
    verlet.build(cfg.system);

    // Move a single atom by 0.3*skin (< skin/2). The sum of the two largest
    // displacements is 0.3*skin + 0 = 0.3*skin, well under the skin.
    Vec3 r0 = cfg.system.position(0);
    cfg.system.setPosition(0, r0 + Vec3{0.3 * kSkin, 0.0, 0.0});

    CHECK_FALSE(verlet.needsRebuild(cfg.system));
}

TEST_CASE("the sum of the two largest displacements past the skin triggers a rebuild") {
    auto cfg = readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic1.txt");
    VerletList verlet(cfg.system.box(), 3.0, kSkin);
    verlet.build(cfg.system);

    // Move two atoms by 0.6*skin each: neither alone exceeds the skin, but
    // their sum (1.2*skin) does. A rebuild criterion that only looks at the
    // single largest displacement would wrongly say no rebuild is needed.
    Vec3 r0 = cfg.system.position(0);
    Vec3 r1 = cfg.system.position(1);
    cfg.system.setPosition(0, r0 + Vec3{0.6 * kSkin, 0.0, 0.0});
    cfg.system.setPosition(1, r1 + Vec3{0.0, 0.6 * kSkin, 0.0});

    CHECK(verlet.needsRebuild(cfg.system));
}

TEST_CASE("Verlet list forces still match brute force after 200 MD steps") {
    auto cfg = readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic1.txt");
    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);
    VerletList verlet(cfg.system.box(), 3.0, kSkin);
    verlet.build(cfg.system);

    // A minimal hand-rolled velocity-Verlet loop -- Task 9 hasn't landed a
    // general integrator yet, and none is needed here: this test only
    // exercises the neighbour list under realistic dynamics, not the
    // integrator's own correctness.
    const double dt = 2e-3;
    const std::size_t n = cfg.system.size();
    auto forces = lj.computeForces(cfg.system, verlet);

    for (int step = 0; step < 200; ++step) {
        for (std::size_t i = 0; i < n; ++i) {
            Vec3 v = cfg.system.velocity(i) + (0.5 * dt / cfg.system.mass(i)) * cfg.system.force(i);
            cfg.system.setVelocity(i, v);
            cfg.system.setPosition(i, cfg.system.position(i) + dt * v);
        }
        verlet.rebuildIfNeeded(cfg.system);
        forces = lj.computeForces(cfg.system, verlet);
        for (std::size_t i = 0; i < n; ++i) {
            Vec3 v = cfg.system.velocity(i) + (0.5 * dt / cfg.system.mass(i)) * cfg.system.force(i);
            cfg.system.setVelocity(i, v);
        }
    }
    (void)forces;

    std::vector<Vec3> fast;
    for (std::size_t i = 0; i < n; ++i) fast.push_back(cfg.system.force(i));

    auto brute = lj.computeForces(cfg.system);
    CHECK(brute.energy == doctest::Approx(forces.energy).epsilon(1e-12));
    CHECK(brute.virial == doctest::Approx(forces.virial).epsilon(1e-12));
    for (std::size_t i = 0; i < n; ++i) {
        CAPTURE(i);
        CHECK(cfg.system.force(i).x == doctest::Approx(fast[i].x).epsilon(1e-12));
        CHECK(cfg.system.force(i).y == doctest::Approx(fast[i].y).epsilon(1e-12));
        CHECK(cfg.system.force(i).z == doctest::Approx(fast[i].z).epsilon(1e-12));
    }

    // The rebuild trigger must actually have exercised its logic: 0 means it
    // never fired (the test proves nothing about the trigger), and 200 would
    // mean the skin is doing no work (rebuilding every step).
    CHECK(verlet.rebuildCount() > 0);
    CHECK(verlet.rebuildCount() < 200);
}
