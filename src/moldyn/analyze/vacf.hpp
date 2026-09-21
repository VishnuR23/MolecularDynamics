#pragma once

#include <cstddef>
#include <vector>

#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"

namespace moldyn {

// Velocity autocorrelation function, accumulated over all time origins for
// each lag.
class Vacf {
public:
    Vacf(std::size_t natoms, std::size_t maxLag);

    void push(const System& sys);

    std::vector<double> vacf() const;  // normalised to 1 at lag 0, size maxLag + 1

    // Green-Kubo: D = (1/3) * integral of the UNNORMALISED VACF over time,
    // by the trapezoid rule.
    double diffusionCoefficient(double dt) const;

private:
    std::vector<double> unnormalised() const;

    std::size_t natoms_;
    std::size_t maxLag_;
    std::size_t period_;  // maxLag_ + 1: the ring buffer length

    std::vector<std::vector<Vec3>> ring_;  // ring_[slot][atom]: velocity at that push
    std::vector<double> sumDot_;           // indexed by lag
    std::vector<std::size_t> count_;       // indexed by lag
    std::size_t nPushed_ = 0;
};

}  // namespace moldyn
