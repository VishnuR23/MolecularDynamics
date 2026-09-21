#include "moldyn/util/lattice.hpp"

#include <array>
#include <cmath>

#include "moldyn/core/box.hpp"
#include "moldyn/util/random.hpp"

namespace moldyn {

namespace {
constexpr double kMass = 1.0;  // reduced-unit convention (see io/nist_config.hpp)
}  // namespace

System fccLattice(int cells, double density, double temperature, uint64_t seed) {
    const auto cellsD = static_cast<double>(cells);
    const int n = 4 * cells * cells * cells;
    const double length = std::cbrt(static_cast<double>(n) / density);
    const double a = length / cellsD;

    // Four-atom fcc basis, in units of the lattice constant.
    const std::array<Vec3, 4> basis = {{
        Vec3{0.0, 0.0, 0.0},
        Vec3{0.0, 0.5, 0.5},
        Vec3{0.5, 0.0, 0.5},
        Vec3{0.5, 0.5, 0.0},
    }};

    Box box(Vec3{length, length, length});
    System sys(box);
    sys.reserve(static_cast<std::size_t>(n));

    Rng rng(seed);
    const double sigmaV = std::sqrt(temperature / kMass);

    for (int ix = 0; ix < cells; ++ix) {
        for (int iy = 0; iy < cells; ++iy) {
            for (int iz = 0; iz < cells; ++iz) {
                for (const Vec3& b : basis) {
                    // Fills [0, L) per axis, then shift to the project's
                    // [-L/2, +L/2) convention; wrap folds the -L/2 boundary
                    // case correctly (Box::wrap is now exact there).
                    Vec3 r{
                        a * (static_cast<double>(ix) + b.x) - 0.5 * length,
                        a * (static_cast<double>(iy) + b.y) - 0.5 * length,
                        a * (static_cast<double>(iz) + b.z) - 0.5 * length,
                    };
                    r = box.wrap(r);

                    Vec3 v{
                        sigmaV * rng.normal(),
                        sigmaV * rng.normal(),
                        sigmaV * rng.normal(),
                    };

                    sys.addAtom(r, v, kMass);
                }
            }
        }
    }

    sys.removeCenterOfMassMotion();
    return sys;
}

}  // namespace moldyn
