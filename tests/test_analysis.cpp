#include "doctest/doctest.h"
#include "moldyn/analyze/block_average.hpp"
#include "moldyn/analyze/msd.hpp"
#include "moldyn/analyze/rdf.hpp"
#include "moldyn/analyze/vacf.hpp"
#include "moldyn/core/box.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/util/random.hpp"

#include <cmath>
#include <cstddef>
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
