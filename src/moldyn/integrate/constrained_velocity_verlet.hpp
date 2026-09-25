#pragma once

#include <cstddef>
#include <vector>

#include "moldyn/constrain/rigid_water.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"

namespace moldyn {

// Velocity-Verlet with rigid-body constraints -- the RATTLE form.
//
// The ordering matters and one part of it is easy to leave out: after the
// position reset, the half-step velocity has to be corrected by the reset
// displacement over dt. Without that correction the geometry still comes out
// rigid and the trajectory still looks plausible, but energy is no longer
// conserved to the integrator's order, because the velocities no longer
// correspond to the positions the constraint actually produced.
//
// Forces must be current on entry, as for the unconstrained integrator.
class ConstrainedVelocityVerlet {
public:
    ConstrainedVelocityVerlet(double dt, RigidWater constraint)
        : dt_(dt), constraint_(constraint) {}

    template <class ForceFn>
    void step(System& sys, ForceFn&& computeForces) {
        const std::size_t n = sys.size();
        reference_.resize(n);
        unconstrained_.resize(n);

        for (std::size_t i = 0; i < n; ++i) {
            reference_[i] = sys.position(i);
        }
        for (std::size_t i = 0; i < n; ++i) {
            sys.setVelocity(i, sys.velocity(i) + sys.force(i) * (0.5 * dt_ / sys.mass(i)));
        }
        for (std::size_t i = 0; i < n; ++i) {
            sys.setPosition(i, sys.position(i) + sys.velocity(i) * dt_);
            unconstrained_[i] = sys.position(i);
        }

        constraint_.constrainPositions(sys, reference_);

        for (std::size_t i = 0; i < n; ++i) {
            const Vec3 shift = sys.box().minimumImage(sys.position(i) - unconstrained_[i]);
            sys.setVelocity(i, sys.velocity(i) + shift * (1.0 / dt_));
        }

        computeForces(sys);

        for (std::size_t i = 0; i < n; ++i) {
            sys.setVelocity(i, sys.velocity(i) + sys.force(i) * (0.5 * dt_ / sys.mass(i)));
        }
        constraint_.constrainVelocities(sys);
    }

    double dt() const { return dt_; }
    const RigidWater& constraint() const { return constraint_; }

private:
    double dt_;
    RigidWater constraint_;
    // Reused across steps so a long run does not reallocate every step.
    mutable std::vector<Vec3> reference_;
    mutable std::vector<Vec3> unconstrained_;
};

}  // namespace moldyn
