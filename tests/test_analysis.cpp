#include "doctest/doctest.h"
#include "moldyn/analyze/block_average.hpp"
#include "moldyn/analyze/msd.hpp"
#include "moldyn/analyze/rdf.hpp"
#include "moldyn/analyze/vacf.hpp"
#include "moldyn/core/box.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/integrate/langevin.hpp"
#include "moldyn/util/random.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

using namespace moldyn;

TEST_CASE("g(r) of an ideal gas is 1 everywhere") {
    Box box(Vec3{20, 20, 20});
    Rng rng(4242);
    Rdf rdf(8.0, 80);
    for (int frame = 0; frame < 40; ++frame) {
        System sys(box);
        for (int i = 0; i < 1500; ++i)
            sys.addAtom(Vec3{(rng.uniform() - 0.5) * 20.0,
                             (rng.uniform() - 0.5) * 20.0,
                             (rng.uniform() - 0.5) * 20.0}, Vec3{0,0,0}, 1.0);
        rdf.accumulate(sys);
    }
    auto g = rdf.g();
    auto r = rdf.binCentres();
    // skip the first few bins, where the shell volume is tiny and noise is large
    for (std::size_t k = 5; k < g.size(); ++k) {
        CAPTURE(r[k]);
        CHECK(g[k] == doctest::Approx(1.0).epsilon(0.06));
    }
}

TEST_CASE("g(r) of two atoms puts all weight in one bin") {
    Box box(Vec3{20, 20, 20});
    Rdf rdf(8.0, 80);
    System sys(box);
    sys.addAtom(Vec3{0.0, 0.0, 0.0}, Vec3{0, 0, 0}, 1.0);
    sys.addAtom(Vec3{3.05, 0.0, 0.0}, Vec3{0, 0, 0}, 1.0);
    rdf.accumulate(sys);

    auto g = rdf.g();
    auto r = rdf.binCentres();
    const double binWidth = 8.0 / 80.0;
    std::size_t expectedBin = static_cast<std::size_t>(3.05 / binWidth);

    for (std::size_t k = 0; k < g.size(); ++k) {
        CAPTURE(r[k]);
        if (k == expectedBin) {
            CHECK(g[k] > 0.0);
        } else {
            CHECK(g[k] == doctest::Approx(0.0));
        }
    }
}

TEST_CASE("Msd ballistic limit: msd(t) = <v^2> t^2 for force-free atoms") {
    Box box(Vec3{1.0e6, 1.0e6, 1.0e6});
    const std::size_t natoms = 3;
    const double dt = 0.1;
    const std::size_t maxLag = 50;

    const std::vector<Vec3> v = {Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 2.0, 0.0}, Vec3{0.5, 0.5, 0.5}};

    System sys(box);
    for (std::size_t i = 0; i < natoms; ++i) {
        sys.addAtom(Vec3{0.0, 0.0, 0.0}, v[i], 1.0);
    }

    double vSqSum = 0.0;
    for (std::size_t i = 0; i < natoms; ++i) vSqSum += norm2(v[i]);
    const double vSqAvg = vSqSum / static_cast<double>(natoms);

    Msd msd(natoms, maxLag);
    for (std::size_t frame = 0; frame <= maxLag; ++frame) {
        const double t = static_cast<double>(frame) * dt;
        for (std::size_t i = 0; i < natoms; ++i) {
            sys.setPosition(i, v[i] * t);
        }
        msd.push(sys);
    }

    const auto m = msd.msd();
    for (std::size_t lag = 0; lag <= maxLag; ++lag) {
        const double t = static_cast<double>(lag) * dt;
        const double expected = vSqAvg * t * t;
        CAPTURE(lag);
        CHECK(m[lag] == doctest::Approx(expected).epsilon(1e-10));
    }
}

TEST_CASE("Msd unwrapping: an atom crossing the periodic boundary keeps accumulating") {
    const double boxLength = 4.0;
    Box box(Vec3{boxLength, boxLength, boxLength});
    const double dt = 0.05;
    const double vx = 1.7;
    const std::size_t maxLag = 1500;
    const std::size_t totalFrames = 2000;

    System sys(box);
    sys.addAtom(Vec3{0.0, 0.0, 0.0}, Vec3{vx, 0.0, 0.0}, 1.0);

    Msd msd(1, maxLag);
    for (std::size_t frame = 0; frame < totalFrames; ++frame) {
        const double t = static_cast<double>(frame) * dt;
        const Vec3 trueR{vx * t, 0.0, 0.0};
        sys.setPosition(0, box.wrap(trueR));  // simulate storage of wrapped coordinates
        msd.push(sys);
    }

    // Sanity: many boundary crossings actually occurred over the run.
    const double totalDistance = vx * static_cast<double>(totalFrames - 1) * dt;
    const int crossings = static_cast<int>(totalDistance / boxLength);
    REQUIRE(crossings > 10);

    const auto m = msd.msd();
    for (std::size_t lag = 0; lag <= maxLag; lag += 100) {
        const double t = static_cast<double>(lag) * dt;
        const double expected = vx * vx * t * t;
        CAPTURE(lag);
        CHECK(m[lag] == doctest::Approx(expected).epsilon(1e-8));
    }
    // The MSD keeps growing well past anything a wrapped (non-unwrapped)
    // difference could produce -- that failure mode is bounded by roughly
    // (boxLength/2)^2, since a wrapped coordinate can never differ from
    // another by more than half the box.
    CHECK(m[maxLag] > boxLength * boxLength);
}

TEST_CASE("Vacf of free particles is exactly 1 at every lag") {
    Box box(Vec3{100.0, 100.0, 100.0});
    const std::size_t natoms = 4;
    const std::size_t maxLag = 40;

    System sys(box);
    sys.addAtom(Vec3{0, 0, 0}, Vec3{1.0, 0.0, 0.0}, 1.0);
    sys.addAtom(Vec3{0, 0, 0}, Vec3{0.0, -2.0, 0.5}, 1.0);
    sys.addAtom(Vec3{0, 0, 0}, Vec3{0.3, 0.3, 0.3}, 1.0);
    sys.addAtom(Vec3{0, 0, 0}, Vec3{-1.5, 0.0, 2.0}, 1.0);

    Vacf vacf(natoms, maxLag);
    for (std::size_t frame = 0; frame <= maxLag + 5; ++frame) {
        vacf.push(sys);  // velocities never change: force-free
    }

    const auto c = vacf.vacf();
    for (std::size_t lag = 0; lag <= maxLag; ++lag) {
        CAPTURE(lag);
        CHECK(c[lag] == doctest::Approx(1.0).epsilon(1e-10));
    }
}

TEST_CASE("Block averaging: independent samples match sigma/sqrt(n) within 15%") {
    Rng rng(99);
    const std::size_t n = 8192;
    const double sigma = 2.0;
    std::vector<double> samples(n);
    for (std::size_t i = 0; i < n; ++i) {
        samples[i] = sigma * rng.normal();
    }

    const auto result = blockAverage(samples);
    const double naive = sigma / std::sqrt(static_cast<double>(n));
    CAPTURE(naive);
    CAPTURE(result.standardError);
    CHECK(result.standardError == doctest::Approx(naive).epsilon(0.15));
}

TEST_CASE("Block averaging: a strongly correlated series gives a larger error than naive") {
    Rng rng(1234);
    const std::size_t n = 16384;
    const double phi = 0.95;
    std::vector<double> samples(n);
    samples[0] = rng.normal();
    for (std::size_t i = 1; i < n; ++i) {
        samples[i] = phi * samples[i - 1] + std::sqrt(1.0 - phi * phi) * rng.normal();
    }

    // Naive standard error, computed directly from the raw (unblocked) samples.
    double mean = 0.0;
    for (double x : samples) mean += x;
    mean /= static_cast<double>(n);
    double sumSq = 0.0;
    for (double x : samples) sumSq += (x - mean) * (x - mean);
    const double variance = sumSq / static_cast<double>(n - 1);
    const double naive = std::sqrt(variance / static_cast<double>(n));

    const auto result = blockAverage(samples);
    CAPTURE(naive);
    CAPTURE(result.standardError);
    CHECK(result.standardError > naive * 2.0);
}

// Msd::diffusionCoefficient and Vacf::diffusionCoefficient had no coverage
// at all before this test: every check elsewhere validates a curve (g(r),
// msd(), vacf()), never a computed D. Neither formula being individually
// wrong is caught by comparing them to each other -- two implementations
// can share the same sign flip or the same missing factor. This test
// instead compares both against an exact analytic result.
//
// Non-interacting particles under the Langevin thermostat undergo an
// Ornstein-Uhlenbeck process in velocity, for which the diffusion
// coefficient is known exactly:
//
//   <v(0).v(t)> = 3*(kT/m) * exp(-friction*t)
//   D  =  (1/3) * integral_0^inf <v(0).v(t)> dt  =  kT / (m*friction)
//   MSD(t) -> 6*D*t   for t >> 1/friction
//
// With kT = m = friction = 1, D = 1 exactly. Forces are zeroed every step
// (LangevinBAOAB::step takes a force callback), so this is pure OU
// dynamics: no force field is evaluated.
//
// Note the O step of BAOAB (Leimkuhler-Matthews) is the *exact* OU
// transition kernel over one step of length dt, not an Euler
// approximation -- v_{n+1} = exp(-friction*dt)*v_n + noise is exact for
// any dt. So there is no discretisation bias in the velocity process
// itself; the only errors are (a) Green-Kubo truncation at maxLag, (b)
// trapezoid-rule quadrature error, (c) the MSD fit's ballistic-to-
// diffusive transition bias if the fit window starts too early, and (d)
// finite-sample statistical noise. Parameters below are chosen to make
// (a)-(c) negligible and (d) small, rather than widening the tolerance.
TEST_CASE("Diffusion coefficient of a free Ornstein-Uhlenbeck particle matches D = kT/(m*friction)") {
    const std::size_t natoms = 6000;
    const double dt = 0.2;
    const double friction = 1.0;
    const double kT = 1.0;
    const double exactD = kT / friction;  // mass = 1

    // maxLag*dt = 20 correlation times: Green-Kubo truncation error is
    // ~exp(-20) ~ 2e-9, not a percent.
    const std::size_t maxLag = 100;
    // 20 correlation times of equilibration before any frame is pushed.
    const std::size_t equilSteps = 100;
    // 150 time units (750 * dt) of production, giving each lag up to
    // maxLag hundreds of atoms times hundreds of time origins to average.
    const std::size_t totalSteps = 750;
    const std::uint64_t seed = 20260921;

    Box box(Vec3{1000.0, 1000.0, 1000.0});
    System sys(box);
    for (std::size_t i = 0; i < natoms; ++i) {
        sys.addAtom(Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 0.0}, 1.0);
    }

    // No force field: zero forces every step, so the thermostat evolves
    // pure free-particle OU dynamics.
    auto forceFn = [](System& s) {
        s.zeroForces();
        return 0.0;
    };
    sys.zeroForces();  // satisfy step()'s "forces current on entry" precondition

    LangevinBAOAB thermostat(dt, friction, kT, seed);
    for (std::size_t i = 0; i < equilSteps; ++i) {
        thermostat.step(sys, forceFn);
    }

    Msd msd(natoms, maxLag);
    Vacf vacf(natoms, maxLag);
    msd.push(sys);
    vacf.push(sys);
    for (std::size_t i = 0; i < totalSteps; ++i) {
        thermostat.step(sys, forceFn);
        msd.push(sys);
        vacf.push(sys);
    }

    // MSD fit window: start well past the ballistic-to-diffusive crossover
    // (t = 10/friction = 10 correlation times, so the exp(-friction*t)
    // curvature term contributes < 5e-5 of bias) and run to maxLag.
    const std::size_t fitFrom = static_cast<std::size_t>(10.0 / dt);
    const std::size_t fitTo = maxLag + 1;

    const double dMsd = msd.diffusionCoefficient(dt, fitFrom, fitTo);
    const double dVacf = vacf.diffusionCoefficient(dt);

    CAPTURE(dMsd);
    CAPTURE(dVacf);
    CAPTURE(exactD);

    // Tolerance: measured deviations at these parameters are ~0.01% (MSD)
    // and ~0.14% (VACF); 3% leaves more than a 20x margin over that noise
    // while still catching the bug classes this test exists for -- a sign
    // flip, a missing factor of 3 (normalised vs unnormalised VACF in the
    // Green-Kubo integral), or a wrong divisor in the Einstein slope --
    // every one of which mispredicts D by tens of percent or more.
    CHECK(dMsd == doctest::Approx(exactD).epsilon(0.03));
    CHECK(dVacf == doctest::Approx(exactD).epsilon(0.03));
    CHECK(std::fabs(dMsd - dVacf) < 0.03 * exactD);

    // Cheap extra check on the correlation machinery feeding the
    // Green-Kubo integral: the normalised VACF itself should follow
    // exp(-friction*t).
    const auto c = vacf.vacf();
    for (std::size_t lag = 0; lag <= maxLag; lag += 10) {
        const double t = static_cast<double>(lag) * dt;
        const double expected = std::exp(-friction * t);
        CAPTURE(lag);
        CHECK(c[lag] == doctest::Approx(expected).epsilon(0.05));
    }
}
