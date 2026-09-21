#pragma once

#include <cmath>
#include <cstddef>

#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"

namespace moldyn {

namespace detail {
// a*b + c with a single rounding (std::fma), componentwise. Used for the
// kick and drift updates below so each step introduces one rounding per
// component instead of two (a separate multiply then add), which keeps
// velocity-Verlet's algebraic time-reversibility as close to exact as
// double precision allows.
inline Vec3 fmaVec3(double a, Vec3 b, Vec3 c) {
    return Vec3{
        std::fma(a, b.x, c.x),
        std::fma(a, b.y, c.y),
        std::fma(a, b.z, c.z),
    };
}
}  // namespace detail

// Kick-drift-kick velocity-Verlet, the standard second-order symplectic
// integrator for Newtonian dynamics:
//
//   v += (dt/2) * f/m
//   r += dt * v ; wrap into the box
//   f  = computeForces(sys)     // recomputed HERE, between the half-kicks
//   v += (dt/2) * f/m
//
// step() assumes forces already in sys are current for the position sys
// holds on entry -- callers must call computeForces(sys) once before the
// first step() (see tests/test_integrator.cpp). Every subsequent call
// leaves sys with forces current for its new position, so steps chain
// without recomputing forces outside the loop.
class VelocityVerlet {
public:
    explicit VelocityVerlet(double dt);

    // One integration step of length dt(). `computeForces` must fill sys's
    // per-atom forces (sys.setForce/addForce) and return the potential
    // energy; it is called exactly once per step, between the two
    // half-kicks, so the force it computes is the current one used for the
    // second half-kick and left current on exit.
    template <class ForceFn>
    void step(System& sys, ForceFn&& computeForces) const {
        const std::size_t n = sys.size();
        const double halfDt = 0.5 * dt_;

        // First half-kick: v += (dt/2) * f/m, using the forces already
        // current in sys on entry.
        for (std::size_t i = 0; i < n; ++i) {
            const double invMass = 1.0 / sys.mass(i);
            sys.setVelocity(i, detail::fmaVec3(halfDt * invMass, sys.force(i), sys.velocity(i)));
        }

        // Drift: r += dt * v, then wrap into the box.
        for (std::size_t i = 0; i < n; ++i) {
            Vec3 r = detail::fmaVec3(dt_, sys.velocity(i), sys.position(i));
            sys.setPosition(i, sys.box().wrap(r));
        }

        // Forces at the new positions.
        computeForces(sys);

        // Second half-kick, with the freshly recomputed forces.
        for (std::size_t i = 0; i < n; ++i) {
            const double invMass = 1.0 / sys.mass(i);
            sys.setVelocity(i, detail::fmaVec3(halfDt * invMass, sys.force(i), sys.velocity(i)));
        }
    }

    double dt() const;

private:
    double dt_;
};

}  // namespace moldyn
