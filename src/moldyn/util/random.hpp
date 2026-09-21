#pragma once

#include <cmath>
#include <cstdint>
#include <random>

namespace moldyn {

// Deterministic, cross-platform random source.
//
// std::mt19937_64 is fully specified by the C++ standard, so a given seed
// produces the same bit sequence on every platform and standard library.
// std::normal_distribution is NOT specified that way -- libstdc++ (Linux)
// and libc++ (macOS) use different algorithms, so the same seed and the
// same std::mt19937_64 stream would still yield different normal deviates
// on the two platforms this project's CI runs. Since reproducible runs are
// this project's premise, normal() is implemented here directly as an
// explicit Box-Muller transform over uniform(), which only depends on
// std::mt19937_64 and is therefore identical on every platform.
//
// It uses Marsaglia's polar variant (rejection sampling in the unit disc)
// rather than the textbook trigonometric form. Both are standard Box-Muller
// variants; the polar form needs only sqrt and log, not sin/cos. That
// matters here: IEEE 754 requires sqrt to be correctly rounded on every
// conformant platform, but log, sin and cos are not -- libm implementations
// are only required to be within about 1 ULP, so their last bit can
// legitimately differ between platforms. Dropping sin/cos removes two of
// the three sources of that residual, platform-dependent last-bit risk.
// normal() still calls log, though, so this project's bit-for-bit
// reproducibility claim holds for a given seed on the same binary, but not
// necessarily across different libm implementations (e.g. Linux vs
// macOS): a last-ULP difference in log() there can still perturb the
// initial velocities by 1 ULP. The statistical checks this project cares
// about (dt^2 convergence, reversibility to machine precision) are
// insensitive to that; exact bit-for-bit cross-platform reproduction of a
// given seed's trajectory is not claimed.
class Rng {
public:
    explicit Rng(uint64_t seed) : engine_(seed) {}

    // Uniform deviate in (0, 1]. std::uniform_real_distribution draws from
    // [0, 1), which can return exactly 0; flip to 1 - u so the result is in
    // (0, 1] instead, since Box-Muller below needs log(u) with u > 0.
    double uniform() {
        return 1.0 - dist_(engine_);
    }

    // Standard normal deviate, N(0, 1), via Marsaglia's polar Box-Muller.
    // The rejection loop's trip count is itself a deterministic function of
    // the uniform() stream, so for a given seed it consumes the same
    // number of uniforms (and produces the same deviates) every run.
    double normal() {
        if (haveSpare_) {
            haveSpare_ = false;
            return spare_;
        }
        double u = 0.0;
        double v = 0.0;
        double s = 0.0;
        do {
            u = 2.0 * uniform() - 1.0;
            v = 2.0 * uniform() - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        const double multiplier = std::sqrt(-2.0 * std::log(s) / s);
        spare_ = v * multiplier;
        haveSpare_ = true;
        return u * multiplier;
    }

private:
    std::mt19937_64 engine_;
    std::uniform_real_distribution<double> dist_{0.0, 1.0};  // [0, 1)
    double spare_ = 0.0;
    bool haveSpare_ = false;
};

}  // namespace moldyn
