#include "moldyn/forcefield/spce_forcefield.hpp"

#include <cstddef>
#include <stdexcept>

namespace moldyn {

SpceForceField::SpceForceField(const SpceParameters& params, double cutoff, int kmax)
    : params_(params),
      lj_(params.sigmaO, params.epsilonO, cutoff, Truncation::Truncated),
      ewald_(cutoff, kmax) {}

System SpceForceField::gatherOxygens(const System& sys) const {
    if (sys.size() % 3 != 0) {
        throw std::invalid_argument(
            "SpceForceField: atom count is not a multiple of 3; expected O,H,H triples");
    }
    System oxygens(sys.box());
    oxygens.reserve(sys.size() / 3);
    for (std::size_t i = 0; i < sys.size(); i += 3) {
        oxygens.addAtom(sys.position(i), Vec3{0.0, 0.0, 0.0}, params_.massOxygen);
    }
    return oxygens;
}

double SpceForceField::computeForces(System& sys) const {
    System oxygens = gatherOxygens(sys);
    const EnergyVirial disp = lj_.computeForces(oxygens);
    const double elrc = lj_.longRangeCorrection(oxygens);

    sys.zeroForces();
    // Scatter the dispersion force back onto the oxygen of each molecule.
    for (std::size_t m = 0, o = 0; m < sys.size(); m += 3, ++o) {
        sys.addForce(m, oxygens.force(o));
    }

    // Ewald accumulates onto whatever is already there, so this must follow
    // the scatter rather than zeroing over it.
    const EwaldTerms elec = ewald_.accumulateForces(sys);
    return disp.energy + elrc + elec.total();
}

double SpceForceField::computeEnergy(const System& sys) const {
    const System oxygens = gatherOxygens(sys);
    const double disp = lj_.computeEnergyVirial(oxygens).energy;
    const double elrc = lj_.longRangeCorrection(oxygens);
    return disp + elrc + ewald_.computeEnergies(sys).total();
}

}  // namespace moldyn
