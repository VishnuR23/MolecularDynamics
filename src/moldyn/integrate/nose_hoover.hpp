#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "moldyn/core/system.hpp"
#include "moldyn/integrate/velocity_verlet.hpp"

namespace moldyn {

// Nose-Hoover chain thermostat (Martyna, Tuckerman, Tobias & Klein, 1996),
// coupled to a velocity-Verlet core. The chain is propagated with a
// reversible, Suzuki-Yoshida decomposed update (Yoshida 1990) applied as a
// half-step of length dt/2 before and after the velocity-Verlet core step,
// following the standard second-order Trotter factorisation of the
// extended-system Liouville operator:
//
//   step(dt) = chain(dt/2) . velocityVerlet(dt) . chain(dt/2)
//
// Each chain half-step is itself split into the (nested, reversible)
// updates of the thermostat momenta p_xi, a drift of the thermostat
// positions xi, and a rescaling of the physical velocities by
// exp(-delta * p_xi_1 / Q_1); see nose_hoover.cpp for the derivation this
// follows.
class NoseHooverChain {
public:
    // chainLength >= 1; 3 is the usual choice. tau is the thermostat
    // period; Q_1 = degreesOfFreedom*kT*tau^2, Q_i = kT*tau^2 for i > 1.
    NoseHooverChain(double dt, double temperature, double tau, int chainLength,
                     std::size_t degreesOfFreedom);

    // One integration step of length dt(). `computeForces` must fill sys's
    // per-atom forces and return the potential energy; it is called exactly
    // once per step, inside the velocity-Verlet core.
    template <class ForceFn>
    void step(System& sys, ForceFn&& computeForces) {
        integrateChain(sys, 0.5 * dt_);
        core_.step(sys, std::forward<ForceFn>(computeForces));
        integrateChain(sys, 0.5 * dt_);
    }

    // The extended-system conserved quantity contribution. Add this to the
    // physical total energy (PE + KE) and the sum must be conserved:
    //   sum_i p_xi_i^2 / (2 Q_i) + Nf*kT*xi_1 + kT * sum_{i>1} xi_i
    double conservedQuantityContribution() const;

private:
    void integrateChain(System& sys, double delta);
    void integrateChainSubstep(System& sys, double delta);
    void integratePEtaLink(System& sys, int j, double delta2, double delta4);
    void integrateEta(double delta);
    void scaleVelocities(System& sys, double delta) const;

    double dt_;
    double temperature_;
    double tau_;
    int chainLength_;
    std::size_t dof_;

    VelocityVerlet core_;

    std::vector<double> Q_;    // thermostat "masses"
    std::vector<double> xi_;   // thermostat positions
    std::vector<double> pxi_;  // thermostat momenta
};

}  // namespace moldyn
