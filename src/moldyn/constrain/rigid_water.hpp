#pragma once

#include <cstddef>
#include <vector>

#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"

namespace moldyn {

// Rigid three-site water, held by three distance constraints per molecule:
// O-H1, O-H2 and H1-H2. Positions are reset with SHAKE, velocities with
// RATTLE.
//
// SPC/E is a rigid model -- every published result for it assumes rigidity --
// and the constraint is also what makes a sane timestep possible. An
// unconstrained O-H stretch has a period around 9 fs, forcing dt near 0.1 fs;
// with the bonds frozen, 2 fs is routine.
//
// SETTLE (Miyamoto & Kollman 1992) solves these same three constraints in
// closed form by exploiting the specific geometry, which is faster. It is a
// performance specialisation, not different physics, and it belongs with the
// rest of the performance work. When it lands it must reproduce this class to
// machine precision -- that agreement is how it gets validated.
//
// Molecules are triples of consecutive atoms in O, H, H order, which is what
// buildSpceSystem produces and what the NIST configuration files store.
class RigidWater {
public:
    // `hohAngleDegrees` is the H-O-H angle; the H-H target follows from it.
    // `tolerance` is on the relative constraint violation, |r^2 - d^2| / d^2.
    RigidWater(double ohDistance,
               double hohAngleDegrees,
               double tolerance = 1e-12,
               int maxIterations = 500);

    // Reset positions onto the constraint surface. `reference` holds the
    // positions at the start of the step, which must already satisfy the
    // constraints -- SHAKE's corrections act along those reference bond
    // directions, which is what makes this equivalent to constraint forces
    // rather than an arbitrary geometric fit.
    //
    // Returns the iteration count. Throws std::runtime_error if it does not
    // converge, rather than silently returning an unconstrained system.
    int constrainPositions(System& sys, const std::vector<Vec3>& reference) const;

    // Remove the velocity component along each constraint direction. Skipping
    // this leaves the geometry rigid while velocities drift off the constraint
    // surface, which quietly corrupts the temperature and anything computed
    // from velocities.
    int constrainVelocities(System& sys) const;

    double ohDistance() const { return oh_; }
    double hhDistance() const { return hh_; }

    // A rigid three-site molecule has 9 Cartesian coordinates and 3
    // constraints, so 6 degrees of freedom -- not 9. Removing centre-of-mass
    // motion freezes 3 more.
    static std::size_t degreesOfFreedom(std::size_t molecules, bool comRemoved);

private:
    double oh_;
    double hh_;
    double tolerance_;
    int maxIterations_;
};

}  // namespace moldyn
