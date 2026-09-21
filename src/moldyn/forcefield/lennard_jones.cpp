#include "moldyn/forcefield/lennard_jones.hpp"

#include <cmath>
#include <cstddef>

#include "moldyn/core/vec3.hpp"
#include "moldyn/neighbor/cell_list.hpp"

namespace moldyn {

namespace {
constexpr double kPi = 3.14159265358979323846;

// Pair source for the original O(N^2) path: every i<j pair. Kept local to
// this translation unit so the public brute-force API is unchanged.
struct BruteForcePairs {
    std::size_t n;

    template <class Fn>
    void forEachPair(Fn&& fn) const {
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                fn(i, j);
            }
        }
    }
};
}  // namespace

LennardJones::LennardJones(double sigma, double epsilon, double cutoff, Truncation truncation)
    : sigma_(sigma),
      epsilon_(epsilon),
      cutoff_(cutoff),
      truncation_(truncation),
      uAtCutoff_(0.0),
      dudrAtCutoff_(0.0) {
    const double sig2 = sigma_ * sigma_;
    const double rc2 = cutoff_ * cutoff_;
    const double s2 = sig2 / rc2;
    const double s6 = s2 * s2 * s2;

    uAtCutoff_ = 4.0 * epsilon_ * (s6 * s6 - s6);
    const double wAtCutoff = 24.0 * epsilon_ * (2.0 * s6 * s6 - s6);
    // w(r) = -r * du/dr  =>  du/dr|rc = -w(rc) / rc
    dudrAtCutoff_ = -wAtCutoff / cutoff_;
}

template <bool AccumulateForces, class PairSource>
EnergyVirial LennardJones::computePairSum(const System& sys, System* forceSink,
                                           const PairSource& pairs) const {
    EnergyVirial result;

    const double sig2 = sigma_ * sigma_;
    const double rc2 = cutoff_ * cutoff_;

    pairs.forEachPair([&](std::size_t i, std::size_t j) {
        const Vec3 d = sys.box().minimumImage(sys.position(i) - sys.position(j));
        const double r2 = norm2(d);
        if (!(r2 < rc2)) {
            return;
        }

        const double s2 = sig2 / r2;
        const double s6 = s2 * s2 * s2;
        const double uLJ = 4.0 * epsilon_ * (s6 * s6 - s6);
        const double wLJ = 24.0 * epsilon_ * (2.0 * s6 * s6 - s6);

        double u = uLJ;
        double w = wLJ;
        if (truncation_ == Truncation::LinearForceShift) {
            const double r = std::sqrt(r2);
            u = uLJ - uAtCutoff_ - (r - cutoff_) * dudrAtCutoff_;
            w = wLJ + r * dudrAtCutoff_;
        }
        result.energy += u;
        result.virial += w;

        if constexpr (AccumulateForces) {
            const Vec3 f = (w / r2) * d;
            forceSink->addForce(i, f);
            forceSink->addForce(j, -1.0 * f);
        }
    });

    return result;
}

EnergyVirial LennardJones::computeEnergyVirial(const System& sys) const {
    return computePairSum<false>(sys, nullptr, BruteForcePairs{sys.size()});
}

EnergyVirial LennardJones::computeForces(System& sys) const {
    sys.zeroForces();
    return computePairSum<true>(sys, &sys, BruteForcePairs{sys.size()});
}

EnergyVirial LennardJones::computeForces(System& sys, const CellList& cells) const {
    sys.zeroForces();
    return computePairSum<true>(sys, &sys, cells);
}

double LennardJones::longRangeCorrection(const System& sys) const {
    if (truncation_ == Truncation::LinearForceShift) {
        return 0.0;
    }

    const double n = static_cast<double>(sys.size());
    const double rho = n / sys.box().volume();
    const double sigOverRc = sigma_ / cutoff_;
    const double sigOverRc3 = sigOverRc * sigOverRc * sigOverRc;
    const double sigOverRc9 = sigOverRc3 * sigOverRc3 * sigOverRc3;
    const double sig3 = sigma_ * sigma_ * sigma_;

    return (8.0 / 3.0) * kPi * n * rho * epsilon_ * sig3 *
           ((1.0 / 3.0) * sigOverRc9 - sigOverRc3);
}

double LennardJones::sigma() const {
    return sigma_;
}

double LennardJones::epsilon() const {
    return epsilon_;
}

double LennardJones::cutoff() const {
    return cutoff_;
}

Truncation LennardJones::truncation() const {
    return truncation_;
}

}  // namespace moldyn
