#include "doctest/doctest.h"
#include "moldyn/integrate/velocity_verlet.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/util/lattice.hpp"
#include <cmath>
#include <vector>

using namespace moldyn;

TEST_CASE("velocity-Verlet energy error scales as dt^2") {
    const std::vector<double> dts = {0.001, 0.002, 0.004, 0.008};
    const double totalTime = 1.0;          // same physical time for every dt
    std::vector<double> logDt, logErr;

    for (double dt : dts) {
        System sys = fccLattice(3, 0.80, 1.0, 12345);   // N = 108
        LennardJones lj(1.0, 1.0, 2.5, Truncation::LinearForceShift);
        VelocityVerlet integrator(dt);

        auto forces = [&](System& s) { return lj.computeForces(s).energy; };
        double pe = forces(sys);
        const double e0 = pe + sys.kineticEnergy();

        double maxDev = 0.0;
        const int steps = static_cast<int>(std::lround(totalTime / dt));
        for (int n = 0; n < steps; ++n) {
            integrator.step(sys, forces);
            pe = lj.computeForces(sys).energy;
            const double e = pe + sys.kineticEnergy();
            maxDev = std::max(maxDev, std::abs(e - e0));
        }
        logDt.push_back(std::log(dt));
        logErr.push_back(std::log(maxDev / std::abs(e0)));
    }

    // least-squares slope of log(error) against log(dt)
    double mx = 0, my = 0;
    for (std::size_t i = 0; i < dts.size(); ++i) { mx += logDt[i]; my += logErr[i]; }
    mx /= static_cast<double>(dts.size());
    my /= static_cast<double>(dts.size());
    double num = 0, den = 0;
    for (std::size_t i = 0; i < dts.size(); ++i) {
        num += (logDt[i] - mx) * (logErr[i] - my);
        den += (logDt[i] - mx) * (logDt[i] - mx);
    }
    const double slope = num / den;
    CAPTURE(slope);
    CHECK(slope > 1.85);
    CHECK(slope < 2.15);
}

TEST_CASE("velocity-Verlet is time reversible") {
    System sys = fccLattice(3, 0.80, 1.0, 999);
    LennardJones lj(1.0, 1.0, 2.5, Truncation::LinearForceShift);
    VelocityVerlet integrator(0.002);
    auto forces = [&](System& s) { return lj.computeForces(s).energy; };

    std::vector<Vec3> start;
    for (std::size_t i = 0; i < sys.size(); ++i) start.push_back(sys.position(i));

    forces(sys);
    for (int n = 0; n < 200; ++n) integrator.step(sys, forces);
    for (std::size_t i = 0; i < sys.size(); ++i)
        sys.setVelocity(i, sys.velocity(i) * -1.0);
    for (int n = 0; n < 200; ++n) integrator.step(sys, forces);

    // Under periodic boundary conditions a position is only defined modulo
    // L, so raw per-component coordinate equality is the wrong check: an
    // atom that starts exactly on Box::wrap's branch cut (fccLattice puts
    // several there, at -L/2) can come back with its stored coordinate
    // flipped to the opposite face -- the same physical point, differing
    // by exactly one box length -- from a single ULP of rounding. Compare
    // via minimum-image displacement instead, which is well-defined on the
    // physical (mod L) positions. This is not a loosened tolerance: a
    // genuinely broken integrator produces an arbitrary, non-zero
    // displacement, and minimumImage() does not make that vanish. The only
    // thing it stops flagging is a change of periodic representative.
    for (std::size_t i = 0; i < sys.size(); ++i) {
        const Vec3 delta = sys.box().minimumImage(sys.position(i) - start[i]);
        CAPTURE(i);
        CHECK(norm(delta) < 1e-8);
    }
}
