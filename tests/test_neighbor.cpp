#include "doctest/doctest.h"
#include "moldyn/core/box.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/neighbor/cell_list.hpp"
#include "moldyn/neighbor/verlet_list.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/io/nist_config.hpp"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace moldyn;

// Why 1e-12 and not 0, throughout this file.
//
// The neighbour lists are a pure optimisation of *how pairs are found*,
// not a different calculation: every pair inside the cutoff is visited
// exactly once by either path, and each pair's contribution is computed by
// the one shared kernel (LennardJones::computePairSum). So the two paths
// sum an identical multiset of identical doubles -- but in a different
// order, because the cell/Verlet traversal visits pairs grouped by cell
// rather than in i<j order. Floating-point addition is not associative, so
// the sums can differ in the last bits, by roughly
// sqrt(n_pairs) * eps_machine ~ sqrt(1e5) * 2e-16 ~ 1e-13 for the N=800
// configuration used here. 1e-12 is that, with a decade of headroom; 0
// would be asserting an associativity that IEEE-754 does not provide.
//
// Phase 3 (SIMD, threading) will stress this deliberately: vectorised
// accumulation reorders and partially parallelises the same sum, and a
// threaded reduction reorders it again and non-deterministically. Expect
// this bound to need loosening then -- to something like 1e-10 -- and when
// it does, loosen it *with* a recomputed error estimate, not to whatever
// makes the run pass. Nothing about the physics changes; only the summation
// order does, and the tolerance should always state which.
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
    // Same bound tests/test_nist_lj.cpp uses: half a unit in NIST's own
    // fifth published significant figure, absolute. Going through the
    // cell list must not cost a single published digit.
    CHECK(std::abs(ev.energy - (-4.3515e3)) <= 0.05);
    CHECK(std::abs(ev.virial - (-5.6867e2)) <= 0.005);
}

TEST_CASE("cell list reports itself unusable for a too-small box") {
    CellList cells(Box(Vec3{5.0, 5.0, 5.0}), 3.0);   // only 1 cell per side
    CHECK_FALSE(cells.usable());
}

TEST_CASE("computeForces refuses an unusable cell list instead of answering wrongly") {
    // The dangerous case is not that this throws -- it is what would
    // happen if it did not. An unusable cell list still enumerates pairs;
    // they are just the wrong ones (a cell's neighbour set wraps onto
    // itself), so the call would return a plausible energy and plausible
    // forces that are silently wrong. Fail loudly instead.
    System sys(Box(Vec3{5.0, 5.0, 5.0}));
    sys.addAtom(Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 0.0}, 1.0);
    sys.addAtom(Vec3{1.1, 0.0, 0.0}, Vec3{0.0, 0.0, 0.0}, 1.0);

    CellList cells(sys.box(), 3.0);   // only 1 cell per side -> unusable
    REQUIRE_FALSE(cells.usable());
    cells.build(sys);

    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);
    CHECK_THROWS_AS(lj.computeForces(sys, cells), std::invalid_argument);
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
