#include "doctest/doctest.h"
#include "moldyn/integrate/langevin.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/util/lattice.hpp"
#include "moldyn/core/system.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

using namespace moldyn;

namespace {

// Kolmogorov-Smirnov statistic for `sample` (assumed to be independent
// draws) against Normal(0, sigma): sort the sample, compute
// max(i/n - Phi(x_i)) and max(Phi(x_i) - (i-1)/n), and take the larger.
// Phi(x) = 0.5*erfc(-x/sqrt(2)).
double ksStatisticNormal(std::vector<double> sample, double sigma) {
    std::sort(sample.begin(), sample.end());
    const std::size_t n = sample.size();
    double d = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        const double phi = 0.5 * std::erfc(-sample[k] / (sigma * std::sqrt(2.0)));
        const double iOverN = static_cast<double>(k + 1) / static_cast<double>(n);
        const double iMinus1OverN = static_cast<double>(k) / static_cast<double>(n);
        d = std::max(d, std::max(iOverN - phi, phi - iMinus1OverN));
    }
    return d;
}

}  // namespace

TEST_CASE("Langevin BAOAB: equipartition holds after equilibration") {
    const int cells = 3;                 // N = 108
    const double density = 0.80;
    const double dt = 0.005;
    const double friction = 1.0;         // velocity autocorrelation time ~ 1/friction = 1.0
    const double targetT = 1.0;

    System sys = fccLattice(cells, density, targetT, 42);
    LennardJones lj(1.0, 1.0, 2.5, Truncation::LinearForceShift);
    LangevinBAOAB thermostat(dt, friction, targetT, 777);
    auto forces = [&](System& s) { return lj.computeForces(s).energy; };
    forces(sys);

    const std::size_t dof = 3 * sys.size() - 3;

    // Equilibrate for well over a correlation time before sampling.
    const int equilSteps = 3000;
    for (int i = 0; i < equilSteps; ++i) thermostat.step(sys, forces);

    // Running mean of the instantaneous temperature over a long window --
    // this only needs many samples, not independent ones.
    const int sampleSteps = 10000;
    double sumT = 0.0;
    int nSamples = 0;
    for (int i = 0; i < sampleSteps; ++i) {
        thermostat.step(sys, forces);
        sumT += sys.temperature(dof);
        ++nSamples;
    }
    const double meanT = sumT / static_cast<double>(nSamples);

    CAPTURE(meanT);
    CAPTURE(targetT);
    CHECK(meanT > targetT * 0.98);
    CHECK(meanT < targetT * 1.02);
}

TEST_CASE("Langevin BAOAB: velocity distribution matches Maxwell-Boltzmann") {
    const int cells = 3;                 // N = 108
    const double density = 0.80;
    const double dt = 0.005;
    const double friction = 1.0;         // velocity autocorrelation time ~ 1/friction = 1.0
    const double targetT = 1.0;
    const double mass = 1.0;             // fccLattice gives every atom mass 1

    System sys = fccLattice(cells, density, targetT, 43);
    LennardJones lj(1.0, 1.0, 2.5, Truncation::LinearForceShift);
    LangevinBAOAB thermostat(dt, friction, targetT, 778);
    auto forces = [&](System& s) { return lj.computeForces(s).energy; };
    forces(sys);

    const int equilSteps = 3000;
    for (int i = 0; i < equilSteps; ++i) thermostat.step(sys, forces);

    // Sample frames separated by 5/friction (5 correlation times) so
    // consecutive frames of the same atom are statistically independent --
    // pooling every frame would feed the KS test correlated data.
    const int stepsPerCorrelationTime = static_cast<int>(std::lround(1.0 / friction / dt));
    const int frameSpacing = 5 * stepsPerCorrelationTime;
    const int numFrames = 60;

    std::vector<double> velocityComponents;
    velocityComponents.reserve(static_cast<std::size_t>(numFrames) * sys.size() * 3);

    for (int f = 0; f < numFrames; ++f) {
        for (int i = 0; i < frameSpacing; ++i) thermostat.step(sys, forces);
        for (std::size_t a = 0; a < sys.size(); ++a) {
            const Vec3 v = sys.velocity(a);
            velocityComponents.push_back(v.x);
            velocityComponents.push_back(v.y);
            velocityComponents.push_back(v.z);
        }
    }

    const std::size_t n = velocityComponents.size();
    const double sigma = std::sqrt(targetT / mass);
    const double ks = ksStatisticNormal(velocityComponents, sigma);
    const double criticalValue = 1.63 / std::sqrt(static_cast<double>(n));

    CAPTURE(frameSpacing);
    CAPTURE(numFrames);
    CAPTURE(n);
    CAPTURE(ks);
    CAPTURE(criticalValue);
    CHECK(ks < criticalValue);
}

TEST_CASE("Langevin BAOAB relaxes a hot start to the target temperature") {
    const int cells = 3;                 // N = 108
    const double density = 0.80;
    const double dt = 0.005;
    const double friction = 1.0;
    const double hotT = 3.0;
    const double targetT = 1.0;

    System sys = fccLattice(cells, density, hotT, 99);
    LennardJones lj(1.0, 1.0, 2.5, Truncation::LinearForceShift);
    LangevinBAOAB thermostat(dt, friction, targetT, 779);
    auto forces = [&](System& s) { return lj.computeForces(s).energy; };
    forces(sys);

    const std::size_t dof = 3 * sys.size() - 3;

    // Let the hot start relax, then confirm it settled near the target
    // rather than staying near the hot starting temperature.
    const int relaxSteps = 3000;
    for (int i = 0; i < relaxSteps; ++i) thermostat.step(sys, forces);

    const int sampleSteps = 4000;
    double sumT = 0.0;
    for (int i = 0; i < sampleSteps; ++i) {
        thermostat.step(sys, forces);
        sumT += sys.temperature(dof);
    }
    const double finalT = sumT / static_cast<double>(sampleSteps);

    CAPTURE(finalT);
    CAPTURE(hotT);
    CAPTURE(targetT);
    CHECK(finalT > targetT * 0.9);
    CHECK(finalT < targetT * 1.1);
    CHECK(finalT < hotT * 0.5);
}
