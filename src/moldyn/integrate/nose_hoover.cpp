#include "moldyn/integrate/nose_hoover.hpp"

#include <array>
#include <cmath>

namespace moldyn {

namespace {

// Fourth-order Suzuki-Yoshida weights (Yoshida 1990), summing to 1. Used to
// decompose one chain half-step into three reversible sub-steps, each of
// which is itself the nested Trotter split implemented in
// integrateChainSubstep/integratePEtaLink below.
std::array<double, 3> yoshidaWeights() {
    const double w1 = 1.0 / (2.0 - std::cbrt(2.0));
    return {w1, 1.0 - 2.0 * w1, w1};
}

}  // namespace

NoseHooverChain::NoseHooverChain(double dt, double temperature, double tau, int chainLength,
                                  std::size_t degreesOfFreedom)
    : dt_(dt),
      temperature_(temperature),
      tau_(tau),
      chainLength_(chainLength),
      dof_(degreesOfFreedom),
      core_(dt),
      Q_(static_cast<std::size_t>(chainLength)),
      xi_(static_cast<std::size_t>(chainLength), 0.0),
      pxi_(static_cast<std::size_t>(chainLength), 0.0) {
    Q_[0] = static_cast<double>(dof_) * temperature_ * tau_ * tau_;
    for (std::size_t k = 1; k < Q_.size(); ++k) {
        Q_[k] = temperature_ * tau_ * tau_;
    }
}

double NoseHooverChain::conservedQuantityContribution() const {
    double sum = static_cast<double>(dof_) * temperature_ * xi_[0];
    for (std::size_t k = 1; k < xi_.size(); ++k) {
        sum += temperature_ * xi_[k];
    }
    for (std::size_t k = 0; k < pxi_.size(); ++k) {
        sum += 0.5 * pxi_[k] * pxi_[k] / Q_[k];
    }
    return sum;
}

// Propagates the chain Liouville operator for total time `delta` (dt/2, as
// called from step()), via the Suzuki-Yoshida decomposition: apply the same
// reversible sub-step three times, at fractions w1, w2, w1 of delta.
void NoseHooverChain::integrateChain(System& sys, double delta) {
    for (double w : yoshidaWeights()) {
        integrateChainSubstep(sys, w * delta);
    }
}

// One reversible sub-step of length `delta`, following the standard MTK96
// NHC integrator: a backward pass over the chain links (M-1 down to 0)
// updating the thermostat momenta from the coupling/force terms, a drift of
// the thermostat positions, a rescaling of the physical velocities, and a
// mirrored forward pass (0 up to M-1) that picks up the change in kinetic
// energy the rescaling just produced.
void NoseHooverChain::integrateChainSubstep(System& sys, double delta) {
    const double delta2 = 0.5 * delta;
    const double delta4 = 0.25 * delta;

    for (int j = chainLength_ - 1; j >= 0; --j) integratePEtaLink(sys, j, delta2, delta4);
    integrateEta(delta);
    scaleVelocities(sys, delta);
    for (int j = 0; j < chainLength_; ++j) integratePEtaLink(sys, j, delta2, delta4);
}

// Updates p_xi_j: multiplicatively damped by the next link in the chain
// (exp(-delta4 * p_xi_{j+1}/Q_{j+1}) on both sides, held fixed across this
// update -- the standard trick that gives this sub-step an exact solution),
// with the force term G_j added between the two damping factors. Per the
// MTK96 equations of motion, G_j is NOT divided by Q_j -- only the
// coupling term (the exp() factors above) carries a 1/Q_{j+1}:
//   G_0 = 2*KE - Nf*kT
//   G_j = p_xi_{j-1}^2/Q_{j-1} - kT,  j > 0
void NoseHooverChain::integratePEtaLink(System& sys, int j, double delta2, double delta4) {
    const std::size_t jj = static_cast<std::size_t>(j);
    const bool hasNext = j < chainLength_ - 1;

    if (hasNext) {
        pxi_[jj] *= std::exp(-delta4 * pxi_[jj + 1] / Q_[jj + 1]);
    }

    double g = 0.0;
    if (j == 0) {
        g = 2.0 * sys.kineticEnergy() - static_cast<double>(dof_) * temperature_;
    } else {
        g = pxi_[jj - 1] * pxi_[jj - 1] / Q_[jj - 1] - temperature_;
    }
    pxi_[jj] += delta2 * g;

    if (hasNext) {
        pxi_[jj] *= std::exp(-delta4 * pxi_[jj + 1] / Q_[jj + 1]);
    }
}

// Drifts every thermostat position by its own velocity p_xi_k/Q_k.
void NoseHooverChain::integrateEta(double delta) {
    for (std::size_t k = 0; k < xi_.size(); ++k) {
        xi_[k] += delta * pxi_[k] / Q_[k];
    }
}

// Rescales every physical velocity by exp(-delta * p_xi_1/Q_1); this is the
// operator that actually couples the chain to the physical system.
void NoseHooverChain::scaleVelocities(System& sys, double delta) const {
    const double factor = std::exp(-delta * pxi_[0] / Q_[0]);
    for (std::size_t i = 0; i < sys.size(); ++i) {
        sys.setVelocity(i, sys.velocity(i) * factor);
    }
}

}  // namespace moldyn
