#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/util/random.hpp"

namespace moldyn {

// Langevin thermostat integrated with the BAOAB splitting (Leimkuhler &
// Matthews 2013), which has the smallest configurational sampling bias of
// the common Langevin splittings:
//
//   B: v += (dt/2) * f/m
//   A: r += (dt/2) * v ; wrap into the box
//   O: v = c1*v + c2*sqrt(kT/m)*xi     c1 = exp(-friction*dt), c2 = sqrt(1-c1^2)
//   A: r += (dt/2) * v ; wrap into the box
//      f = computeForces(sys)
//   B: v += (dt/2) * f/m
//
// step() assumes forces already in sys are current for the position sys
// holds on entry (as VelocityVerlet::step does) -- callers must call
// computeForces(sys) once before the first step(). Every subsequent call
// leaves sys with forces current for its new position.
class LangevinBAOAB {
public:
    LangevinBAOAB(double dt, double friction, double temperature, uint64_t seed);

    // One integration step of length dt(). `computeForces` must fill sys's
    // per-atom forces and return the potential energy; it is called exactly
    // once per step, between the second A and the final B.
    template <class ForceFn>
    void step(System& sys, ForceFn&& computeForces) {
        const std::size_t n = sys.size();
        const double halfDt = 0.5 * dt_;
        const double c1 = std::exp(-friction_ * dt_);
        const double c2 = std::sqrt(1.0 - c1 * c1);

        // B: half-kick using the forces already current in sys on entry.
        for (std::size_t i = 0; i < n; ++i) {
            const double invMass = 1.0 / sys.mass(i);
            sys.setVelocity(i, sys.velocity(i) + (halfDt * invMass) * sys.force(i));
        }

        // A: half-drift, wrapped into the box.
        for (std::size_t i = 0; i < n; ++i) {
            sys.setPosition(i, sys.box().wrap(sys.position(i) + halfDt * sys.velocity(i)));
        }

        // O: Ornstein-Uhlenbeck velocity randomization, per component.
        for (std::size_t i = 0; i < n; ++i) {
            const double sigma = std::sqrt(temperature_ / sys.mass(i));
            const Vec3 v = sys.velocity(i);
            const Vec3 xi{rng_.normal(), rng_.normal(), rng_.normal()};
            sys.setVelocity(i, Vec3{c1 * v.x, c1 * v.y, c1 * v.z} + (c2 * sigma) * xi);
        }

        // A: half-drift, wrapped into the box.
        for (std::size_t i = 0; i < n; ++i) {
            sys.setPosition(i, sys.box().wrap(sys.position(i) + halfDt * sys.velocity(i)));
        }

        // Forces at the new positions.
        computeForces(sys);

        // B: half-kick with the freshly recomputed forces.
        for (std::size_t i = 0; i < n; ++i) {
            const double invMass = 1.0 / sys.mass(i);
            sys.setVelocity(i, sys.velocity(i) + (halfDt * invMass) * sys.force(i));
        }
    }

    double temperature() const;

private:
    double dt_;
    double friction_;
    double temperature_;
    Rng rng_;
};

}  // namespace moldyn
