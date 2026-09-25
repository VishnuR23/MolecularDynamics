#include "doctest/doctest.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "moldyn/constrain/rigid_water.hpp"
#include "moldyn/core/box.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/units.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/forcefield/spce_forcefield.hpp"
#include "moldyn/forcefield/spce.hpp"
#include "moldyn/io/nist_config.hpp"
#include "moldyn/util/random.hpp"

using namespace moldyn;

namespace {

SpceSystem loadWater() {
    const std::string path =
        std::string(MOLDYN_DATA_DIR) + "/nist/spce/spce_sample_config_periodic1.txt";
    return buildSpceSystem(readNistConfig(path));
}

RigidWater makeConstraint() {
    const SpceParameters p;
    return RigidWater(p.rOH, p.angleHOH);
}

// Worst relative deviation of any constrained distance from its target.
double worstGeometryError(const System& sys, const RigidWater& rw) {
    const Box& box = sys.box();
    const double targets[3] = {rw.ohDistance(), rw.ohDistance(), rw.hhDistance()};
    const std::size_t pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    double worst = 0.0;
    for (std::size_t m = 0; m < sys.size(); m += 3) {
        for (std::size_t c = 0; c < 3; ++c) {
            const Vec3 d = box.minimumImage(sys.position(m + pairs[c][0]) -
                                            sys.position(m + pairs[c][1]));
            worst = std::max(worst, std::abs(norm(d) - targets[c]) / targets[c]);
        }
    }
    return worst;
}

// One constrained velocity-Verlet step (RATTLE form). The position reset also
// corrects the half-step velocity -- omitting that correction is a classic
// error that leaves the geometry rigid while quietly breaking conservation.
template <class ForceFn>
void constrainedStep(System& sys, double dt, const RigidWater& rw, ForceFn&& forces) {
    const std::size_t n = sys.size();
    std::vector<Vec3> reference(n);
    for (std::size_t i = 0; i < n; ++i) {
        reference[i] = sys.position(i);
    }
    for (std::size_t i = 0; i < n; ++i) {
        sys.setVelocity(i, sys.velocity(i) + sys.force(i) * (0.5 * dt / sys.mass(i)));
    }
    std::vector<Vec3> unconstrained(n);
    for (std::size_t i = 0; i < n; ++i) {
        sys.setPosition(i, sys.position(i) + sys.velocity(i) * dt);
        unconstrained[i] = sys.position(i);
    }
    rw.constrainPositions(sys, reference);
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 shift = sys.box().minimumImage(sys.position(i) - unconstrained[i]);
        sys.setVelocity(i, sys.velocity(i) + shift * (1.0 / dt));
    }
    forces(sys);
    for (std::size_t i = 0; i < n; ++i) {
        sys.setVelocity(i, sys.velocity(i) + sys.force(i) * (0.5 * dt / sys.mass(i)));
    }
    rw.constrainVelocities(sys);
}

}  // namespace

TEST_CASE("degrees of freedom count constraints, not coordinates") {
    // 9 Cartesian coordinates per molecule minus 3 constraints = 6.
    CHECK(RigidWater::degreesOfFreedom(100, false) == 600);
    CHECK(RigidWater::degreesOfFreedom(100, true) == 597);
    CHECK(RigidWater::degreesOfFreedom(1, false) == 6);
}

TEST_CASE("the internal time unit is Angstrom * sqrt(amu / kB)") {
    // 1.0966875353 ps, derived from CODATA kB and the atomic mass unit.
    CHECK(units::kTauInPicoseconds == doctest::Approx(1.0966875353).epsilon(1e-9));
    CHECK(units::femtoseconds(2.0) == doctest::Approx(0.00182367).epsilon(1e-5));
    CHECK(units::toPicoseconds(units::picoseconds(0.5)) == doctest::Approx(0.5));
}

TEST_CASE("the NIST configuration is already rigid to begin with") {
    // If it were not, every later test would be measuring the wrong thing.
    auto w = loadWater();
    CHECK(worstGeometryError(w.full, makeConstraint()) < 1e-6);
}

TEST_CASE("constraining an already-rigid molecule changes nothing") {
    // Idempotence. A projection that "fixes" a configuration already on the
    // constraint surface is doing something it should not.
    auto w = loadWater();
    const RigidWater rw = makeConstraint();
    std::vector<Vec3> before(w.full.size());
    for (std::size_t i = 0; i < w.full.size(); ++i) {
        before[i] = w.full.position(i);
    }
    rw.constrainPositions(w.full, before);
    double worst = 0.0;
    for (std::size_t i = 0; i < w.full.size(); ++i) {
        worst = std::max(worst, norm(w.full.position(i) - before[i]));
    }
    CAPTURE(worst);
    CHECK(worst < 1e-9);
}

TEST_CASE("the constraint conserves each molecule's linear momentum") {
    // Constraint forces are internal, so they cannot move the centre of mass.
    // A correction applied with the wrong mass weighting still produces rigid
    // geometry but pushes the molecule, which this catches.
    auto w = loadWater();
    const RigidWater rw = makeConstraint();
    const std::size_t n = w.full.size();
    std::vector<Vec3> reference(n);
    for (std::size_t i = 0; i < n; ++i) {
        reference[i] = w.full.position(i);
    }
    // Displace every site so the constraints are genuinely violated.
    for (std::size_t i = 0; i < n; ++i) {
        const double s = 0.01 * static_cast<double>((i % 7)) - 0.03;
        w.full.setPosition(i, w.full.position(i) + Vec3{s, -0.5 * s, 0.25 * s});
    }
    std::vector<Vec3> comBefore;
    for (std::size_t m = 0; m < n; m += 3) {
        Vec3 c{0, 0, 0};
        double mass = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            c += w.full.position(m + k) * w.full.mass(m + k);
            mass += w.full.mass(m + k);
        }
        comBefore.push_back(c * (1.0 / mass));
    }
    rw.constrainPositions(w.full, reference);
    double worstShift = 0.0;
    std::size_t idx = 0;
    for (std::size_t m = 0; m < n; m += 3, ++idx) {
        Vec3 c{0, 0, 0};
        double mass = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            c += w.full.position(m + k) * w.full.mass(m + k);
            mass += w.full.mass(m + k);
        }
        worstShift = std::max(worstShift, norm(c * (1.0 / mass) - comBefore[idx]));
    }
    CAPTURE(worstShift);
    CHECK(worstShift < 1e-10);
    CHECK(worstGeometryError(w.full, rw) < 1e-10);
}

TEST_CASE("geometry is held over a trajectory under real forces") {
    auto w = loadWater();
    const RigidWater rw = makeConstraint();
    const SpceForceField ff;
    auto forces = [&](System& s) { ff.computeForces(s); };
    forces(w.full);

    const double dt = units::femtoseconds(2.0);
    double worst = 0.0;
    for (int step = 0; step < 50; ++step) {
        constrainedStep(w.full, dt, rw, forces);
        worst = std::max(worst, worstGeometryError(w.full, rw));
    }
    CAPTURE(worst);
    CHECK(worst < 1e-10);
}

TEST_CASE("energy error still scales as dt^2 with constraints active") {
    // The discriminating test for a constraint: a wrong one still holds bonds
    // rigid and still looks stable, but cannot preserve the integrator's
    // second-order convergence. The equivalent test on the Nose-Hoover chain
    // caught a 56% drift in Phase 1 while its temperature looked perfect.
    //
    // The force field here is deliberately NOT the full SPC/E interaction. It
    // is oxygen-oxygen Lennard-Jones with the force shifted to zero at the
    // cutoff, and no electrostatics. The reason is that NIST's SPC/E
    // parameters are chosen to reproduce static energies, not to conserve
    // energy in dynamics: the dispersion term is truncated unshifted, and the
    // Ewald screening leaves erfc(alpha*rc) = erfc(2.8) = 7.5e-5 of the
    // real-space term still present at the cutoff. Each pair crossing the
    // cutoff therefore steps the energy -- by 0.31 K for dispersion and 0.45 K
    // for a charged pair -- and the number of crossings depends on elapsed
    // physical time, not on dt. That puts a dt-independent floor under the
    // error and flattens the fitted slope no matter how good the constraint
    // is. Measured with the full interaction the slope comes out near 1.2.
    //
    // A continuous force field removes the floor and lets this test measure
    // what it is supposed to measure: the constraint. See
    // docs/findings/ for the truncation note.
    const std::vector<double> stepsFs = {0.5, 1.0, 2.0, 4.0};
    const double totalTimeFs = 100.0;

    std::vector<double> logDt;
    std::vector<double> logErr;

    for (double fs : stepsFs) {
        auto w = loadWater();
        const RigidWater rw = makeConstraint();
        const SpceParameters params;
        const LennardJones lj(params.sigmaO, params.epsilonO, 9.0,
                              Truncation::LinearForceShift);

        // Oxygen-only dispersion, gathered and scattered the way
        // SpceForceField does, but continuous at the cutoff.
        auto forces = [&](System& sys) {
            System ox(sys.box());
            ox.reserve(sys.size() / 3);
            for (std::size_t i = 0; i < sys.size(); i += 3) {
                ox.addAtom(sys.position(i), Vec3{0, 0, 0}, params.massOxygen);
            }
            const double e = lj.computeForces(ox).energy;
            sys.zeroForces();
            for (std::size_t m = 0, o = 0; m < sys.size(); m += 3, ++o) {
                sys.addForce(m, ox.force(o));
            }
            return e;
        };

        Rng rng(2026);
        const double targetT = 298.0;
        for (std::size_t i = 0; i < w.full.size(); ++i) {
            const double sp = std::sqrt(targetT / w.full.mass(i));
            w.full.setVelocity(i, Vec3{sp * rng.normal(), sp * rng.normal(), sp * rng.normal()});
        }
        w.full.removeCenterOfMassMotion();
        rw.constrainVelocities(w.full);

        const double dt = units::femtoseconds(fs);
        double pe = forces(w.full);
        const double e0 = pe + w.full.kineticEnergy();

        double maxDev = 0.0;
        const int steps = static_cast<int>(std::lround(totalTimeFs / fs));
        for (int n = 0; n < steps; ++n) {
            constrainedStep(w.full, dt, rw, forces);
            pe = forces(w.full);
            maxDev = std::max(maxDev, std::abs(pe + w.full.kineticEnergy() - e0));
        }
        logDt.push_back(std::log(fs));
        logErr.push_back(std::log(maxDev / std::abs(e0)));
    }

    double mx = 0.0;
    double my = 0.0;
    for (std::size_t i = 0; i < stepsFs.size(); ++i) {
        mx += logDt[i];
        my += logErr[i];
    }
    mx /= static_cast<double>(stepsFs.size());
    my /= static_cast<double>(stepsFs.size());
    double num = 0.0;
    double den = 0.0;
    for (std::size_t i = 0; i < stepsFs.size(); ++i) {
        num += (logDt[i] - mx) * (logErr[i] - my);
        den += (logDt[i] - mx) * (logDt[i] - mx);
    }
    const double slope = num / den;

    MESSAGE("constrained conservation slope = " << slope);
    for (std::size_t i = 0; i < stepsFs.size(); ++i) {
        MESSAGE("  dt = " << stepsFs[i] << " fs -> relative error " << std::exp(logErr[i]));
    }
    CHECK(slope > 1.85);
    CHECK(slope < 2.15);
}
