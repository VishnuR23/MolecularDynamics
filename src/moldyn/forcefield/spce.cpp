#include "moldyn/forcefield/spce.hpp"

#include <cstddef>
#include <stdexcept>

namespace moldyn {

SpceSystem buildSpceSystem(const NistConfig& cfg, const SpceParameters& params) {
    const std::size_t n = cfg.system.size();
    if (n == 0 || n % 3 != 0) {
        throw std::invalid_argument(
            "buildSpceSystem: atom count must be a positive multiple of 3 (O, H, H per molecule)");
    }
    if (cfg.types.size() != n) {
        throw std::invalid_argument("buildSpceSystem: config has no per-atom types");
    }

    const double qO = -2.0 * params.qH;

    SpceSystem spce{System(cfg.system.box()), System(cfg.system.box()), n / 3};
    spce.full.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t moleculeId = i / 3;
        const Vec3 position = cfg.system.position(i);
        const Vec3 velocity = cfg.system.velocity(i);

        if (cfg.types[i] == "O") {
            spce.full.addAtom(position, velocity, params.massOxygen, qO, moleculeId);
            spce.oxygens.addAtom(position, velocity, params.massOxygen);
        } else if (cfg.types[i] == "H") {
            spce.full.addAtom(position, velocity, params.massHydrogen, params.qH, moleculeId);
        } else {
            throw std::invalid_argument("buildSpceSystem: unrecognised atom type '" + cfg.types[i] +
                                         "' (expected \"O\" or \"H\")");
        }
    }

    return spce;
}

}  // namespace moldyn
