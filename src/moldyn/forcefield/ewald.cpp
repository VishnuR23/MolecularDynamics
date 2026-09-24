#include "moldyn/forcefield/ewald.hpp"

#include <cmath>
#include <cstddef>
#include <unordered_map>
#include <vector>

#include "moldyn/core/vec3.hpp"

namespace moldyn {

namespace {
constexpr double kPi = 3.14159265358979323846;
}  // namespace

double coulombPrefactorKelvinAngstrom() {
    return (kElementaryChargeCoulombs * kElementaryChargeCoulombs) /
           (4.0 * kPi * kVacuumPermittivityFaradPerMeter) * kAngstromsPerMeter /
           kBoltzmannJoulesPerKelvin;
}

Ewald::Ewald(double cutoff, int kmax) : cutoff_(cutoff), kmax_(kmax) {}

double Ewald::cutoff() const {
    return cutoff_;
}

int Ewald::kmax() const {
    return kmax_;
}

template <bool AccumulateForces>
EwaldTerms Ewald::compute(const System& sys, System* forceSink) const {
    EwaldTerms terms;

    const std::size_t n = sys.size();
    const Vec3 lengths = sys.box().lengths();
    const double minLength = std::min({lengths.x, lengths.y, lengths.z});
    const double alpha = 5.6 / minLength;
    const double coulomb = coulombPrefactorKelvinAngstrom();
    const double rc2 = cutoff_ * cutoff_;
    const double twoAlphaOverSqrtPi = 2.0 * alpha / std::sqrt(kPi);

    // --- Real space: intermolecular pairs only, r < cutoff, minimum image.
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (sys.moleculeId(i) == sys.moleculeId(j)) {
                continue;  // handled by the intramolecular correction
            }
            const Vec3 d = sys.box().minimumImage(sys.position(i) - sys.position(j));
            const double r2 = norm2(d);
            if (!(r2 < rc2)) {
                continue;
            }
            const double r = std::sqrt(r2);
            const double qq = sys.charge(i) * sys.charge(j);
            const double erfcAr = std::erfc(alpha * r);
            terms.real += coulomb * qq * erfcAr / r;

            if constexpr (AccumulateForces) {
                const double w =
                    coulomb * qq * (twoAlphaOverSqrtPi * std::exp(-alpha * alpha * r2) + erfcAr / r);
                const Vec3 f = (w / r2) * d;
                forceSink->addForce(i, f);
                forceSink->addForce(j, -1.0 * f);
            }
        }
    }

    // --- Self term: no position dependence, so no force contribution.
    double sumQ2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double q = sys.charge(i);
        sumQ2 += q * q;
    }
    terms.self = -coulomb * alpha / std::sqrt(kPi) * sumQ2;

    // --- Intramolecular correction: minimum-imaged distances within each
    // molecule. Required because some NIST sample configurations store
    // molecules split across the periodic boundary (oxygen on one side,
    // hydrogens wrapped to the other) -- raw separations there reach ~20 A
    // where the bond is 1 A.
    std::unordered_map<std::size_t, std::vector<std::size_t>> molecules;
    for (std::size_t i = 0; i < n; ++i) {
        molecules[sys.moleculeId(i)].push_back(i);
    }
    for (const auto& kv : molecules) {
        const std::vector<std::size_t>& atoms = kv.second;
        for (std::size_t a = 0; a < atoms.size(); ++a) {
            for (std::size_t b = a + 1; b < atoms.size(); ++b) {
                const std::size_t i = atoms[a];
                const std::size_t j = atoms[b];
                const Vec3 d = sys.box().minimumImage(sys.position(i) - sys.position(j));
                const double r2 = norm2(d);
                const double r = std::sqrt(r2);
                const double qq = sys.charge(i) * sys.charge(j);
                const double erfAr = std::erf(alpha * r);
                terms.intra += -coulomb * qq * erfAr / r;

                if constexpr (AccumulateForces) {
                    const double w = coulomb * qq *
                                      (twoAlphaOverSqrtPi * std::exp(-alpha * alpha * r2) - erfAr / r);
                    const Vec3 f = (w / r2) * d;
                    forceSink->addForce(i, f);
                    forceSink->addForce(j, -1.0 * f);
                }
            }
        }
    }

    // --- Reciprocal space: k = 2*pi*n/L, n in [-kmax, kmax]^3 \ {0},
    // included only when |n|^2 < kmax^2 + 2 (27 for kmax = 5, NIST's stated
    // convention).
    const double volume = sys.box().volume();
    const double kx0 = 2.0 * kPi / lengths.x;
    const double ky0 = 2.0 * kPi / lengths.y;
    const double kz0 = 2.0 * kPi / lengths.z;
    const int n2Max = kmax_ * kmax_ + 2;

    const double* xs = sys.xs();
    const double* ys = sys.ys();
    const double* zs = sys.zs();
    const double* qs = sys.qs();

    std::vector<double> cosPhi;
    std::vector<double> sinPhi;
    if constexpr (AccumulateForces) {
        cosPhi.resize(n);
        sinPhi.resize(n);
    }

    for (int nx = -kmax_; nx <= kmax_; ++nx) {
        for (int ny = -kmax_; ny <= kmax_; ++ny) {
            for (int nz = -kmax_; nz <= kmax_; ++nz) {
                const int n2 = nx * nx + ny * ny + nz * nz;
                if (n2 == 0 || n2 >= n2Max) {
                    continue;
                }
                const double kxv = kx0 * static_cast<double>(nx);
                const double kyv = ky0 * static_cast<double>(ny);
                const double kzv = kz0 * static_cast<double>(nz);
                const double k2 = kxv * kxv + kyv * kyv + kzv * kzv;

                double re = 0.0;
                double im = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    const double phase = kxv * xs[i] + kyv * ys[i] + kzv * zs[i];
                    const double c = std::cos(phase);
                    const double s = std::sin(phase);
                    re += qs[i] * c;
                    im += qs[i] * s;
                    if constexpr (AccumulateForces) {
                        cosPhi[i] = c;
                        sinPhi[i] = s;
                    }
                }

                const double factor =
                    coulomb * (2.0 * kPi / volume) * std::exp(-k2 / (4.0 * alpha * alpha)) / k2;
                terms.reciprocal += factor * (re * re + im * im);

                if constexpr (AccumulateForces) {
                    for (std::size_t i = 0; i < n; ++i) {
                        const double coeff =
                            factor * 2.0 * qs[i] * (re * sinPhi[i] - im * cosPhi[i]);
                        forceSink->addForce(i, Vec3{coeff * kxv, coeff * kyv, coeff * kzv});
                    }
                }
            }
        }
    }

    return terms;
}

EwaldTerms Ewald::computeEnergies(const System& sys) const {
    return compute<false>(sys, nullptr);
}

EwaldTerms Ewald::computeForces(System& sys) const {
    sys.zeroForces();
    return compute<true>(sys, &sys);
}

}  // namespace moldyn
