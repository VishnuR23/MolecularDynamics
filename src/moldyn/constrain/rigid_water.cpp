#include "moldyn/constrain/rigid_water.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "moldyn/core/box.hpp"

namespace moldyn {

namespace {

// The three constraints of a rigid three-site molecule, as offsets within the
// triple: O is 0, H1 is 1, H2 is 2.
struct Pair {
    std::size_t i;
    std::size_t j;
};
constexpr Pair kConstraints[3] = {{0, 1}, {0, 2}, {1, 2}};

// M_PI is a POSIX extension, not standard C++, and this project compiles with
// CMAKE_CXX_EXTENSIONS OFF -- which defines __STRICT_ANSI__, under which glibc
// need not declare it at all.
constexpr double kPi = 3.14159265358979323846;

}  // namespace

RigidWater::RigidWater(double ohDistance,
                       double hohAngleDegrees,
                       double tolerance,
                       int maxIterations)
    : oh_(ohDistance),
      hh_(0.0),
      tolerance_(tolerance),
      maxIterations_(maxIterations) {
    if (!(ohDistance > 0.0)) {
        throw std::invalid_argument("RigidWater: ohDistance must be positive");
    }
    if (!(hohAngleDegrees > 0.0) || !(hohAngleDegrees < 180.0)) {
        throw std::invalid_argument("RigidWater: hohAngleDegrees must lie in (0, 180)");
    }
    if (!(tolerance > 0.0)) {
        throw std::invalid_argument("RigidWater: tolerance must be positive");
    }
    if (maxIterations < 1) {
        throw std::invalid_argument("RigidWater: maxIterations must be at least 1");
    }
    // Law of cosines across the H-O-H angle.
    const double halfAngle = 0.5 * hohAngleDegrees * kPi / 180.0;
    hh_ = 2.0 * ohDistance * std::sin(halfAngle);
}

std::size_t RigidWater::degreesOfFreedom(std::size_t molecules, bool comRemoved) {
    const std::size_t perMolecule = 6;  // 9 coordinates - 3 constraints
    const std::size_t total = perMolecule * molecules;
    return comRemoved ? (total >= 3 ? total - 3 : 0) : total;
}

int RigidWater::constrainPositions(System& sys, const std::vector<Vec3>& reference) const {
    const std::size_t n = sys.size();
    if (n % 3 != 0) {
        throw std::runtime_error("RigidWater: atom count is not a multiple of 3");
    }
    if (reference.size() != n) {
        throw std::runtime_error("RigidWater: reference size does not match the system");
    }
    const Box& box = sys.box();
    const double targets[3] = {oh_, oh_, hh_};

    int iterations = 0;
    for (; iterations < maxIterations_; ++iterations) {
        double worst = 0.0;
        for (std::size_t m = 0; m < n; m += 3) {
            for (std::size_t c = 0; c < 3; ++c) {
                const std::size_t i = m + kConstraints[c].i;
                const std::size_t j = m + kConstraints[c].j;
                const double d2 = targets[c] * targets[c];

                // Minimum image throughout: a molecule straddling a periodic
                // boundary has its sites stored on opposite faces, so raw
                // separations would be a box length rather than a bond length.
                const Vec3 rij = box.minimumImage(sys.position(i) - sys.position(j));
                const double s = norm2(rij) - d2;
                const double violation = std::abs(s) / d2;
                if (violation > worst) {
                    worst = violation;
                }
                if (violation <= tolerance_) {
                    continue;
                }

                // The correction acts along the reference bond direction, not
                // the current one. That is what makes this the action of a
                // constraint force rather than a geometric projection.
                const Vec3 rref = box.minimumImage(reference[i] - reference[j]);
                const double invMassSum = 1.0 / sys.mass(i) + 1.0 / sys.mass(j);
                const double denom = 2.0 * invMassSum * dot(rij, rref);
                if (std::abs(denom) < 1e-30) {
                    throw std::runtime_error(
                        "RigidWater: constraint became singular (bond rotated ~90 degrees "
                        "in one step); the timestep is far too large");
                }
                const double g = s / denom;
                sys.setPosition(i, sys.position(i) - rref * (g / sys.mass(i)));
                sys.setPosition(j, sys.position(j) + rref * (g / sys.mass(j)));
            }
        }
        if (worst <= tolerance_) {
            return iterations + 1;
        }
    }
    throw std::runtime_error("RigidWater: SHAKE failed to converge in " +
                             std::to_string(maxIterations_) + " iterations");
}

int RigidWater::constrainVelocities(System& sys) const {
    const std::size_t n = sys.size();
    if (n % 3 != 0) {
        throw std::runtime_error("RigidWater: atom count is not a multiple of 3");
    }
    const Box& box = sys.box();
    const double targets[3] = {oh_, oh_, hh_};

    int iterations = 0;
    for (; iterations < maxIterations_; ++iterations) {
        double worst = 0.0;
        for (std::size_t m = 0; m < n; m += 3) {
            for (std::size_t c = 0; c < 3; ++c) {
                const std::size_t i = m + kConstraints[c].i;
                const std::size_t j = m + kConstraints[c].j;
                const double d2 = targets[c] * targets[c];

                const Vec3 rij = box.minimumImage(sys.position(i) - sys.position(j));
                const Vec3 vij = sys.velocity(i) - sys.velocity(j);

                // A rigid bond cannot lengthen, so the relative velocity must
                // have no component along it: r . v == 0.
                const double rv = dot(rij, vij);
                const double scaled = std::abs(rv) / d2;
                if (scaled > worst) {
                    worst = scaled;
                }
                if (scaled <= tolerance_) {
                    continue;
                }

                const double invMassSum = 1.0 / sys.mass(i) + 1.0 / sys.mass(j);
                const double k = rv / (invMassSum * norm2(rij));
                sys.setVelocity(i, sys.velocity(i) - rij * (k / sys.mass(i)));
                sys.setVelocity(j, sys.velocity(j) + rij * (k / sys.mass(j)));
            }
        }
        if (worst <= tolerance_) {
            return iterations + 1;
        }
    }
    throw std::runtime_error("RigidWater: RATTLE failed to converge in " +
                             std::to_string(maxIterations_) + " iterations");
}

}  // namespace moldyn
