#include "moldyn/analyze/vacf.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "moldyn/core/vec3.hpp"

namespace moldyn {

Vacf::Vacf(std::size_t natoms, std::size_t maxLag)
    : natoms_(natoms),
      maxLag_(maxLag),
      period_(maxLag + 1),
      ring_(period_, std::vector<Vec3>(natoms_)),
      sumDot_(period_, 0.0),
      count_(period_, 0) {}

void Vacf::push(const System& sys) {
    std::vector<Vec3> v(natoms_);
    for (std::size_t i = 0; i < natoms_; ++i) {
        v[i] = sys.velocity(i);
    }

    const std::size_t slot = nPushed_ % period_;
    ring_[slot] = v;

    const std::size_t maxAvailableLag = std::min(nPushed_, maxLag_);
    for (std::size_t lag = 0; lag <= maxAvailableLag; ++lag) {
        const std::size_t originSlot = (slot + period_ - lag) % period_;
        double sum = 0.0;
        for (std::size_t i = 0; i < natoms_; ++i) {
            sum += dot(v[i], ring_[originSlot][i]);
        }
        sumDot_[lag] += sum;
        count_[lag] += natoms_;
    }

    ++nPushed_;
}

std::vector<double> Vacf::unnormalised() const {
    std::vector<double> result(period_, 0.0);
    for (std::size_t lag = 0; lag < period_; ++lag) {
        if (count_[lag] > 0) {
            result[lag] = sumDot_[lag] / static_cast<double>(count_[lag]);
        }
    }
    return result;
}

std::vector<double> Vacf::vacf() const {
    const std::vector<double> u = unnormalised();
    std::vector<double> result(u.size(), 0.0);
    if (u.empty() || u[0] == 0.0) {
        return result;
    }
    const double c0 = u[0];
    for (std::size_t lag = 0; lag < u.size(); ++lag) {
        result[lag] = u[lag] / c0;
    }
    return result;
}

double Vacf::diffusionCoefficient(double dt) const {
    const std::vector<double> u = unnormalised();

    // Integrate only over lags that actually accumulated data.
    std::size_t validLags = 0;
    for (std::size_t lag = 0; lag < period_; ++lag) {
        if (count_[lag] == 0) {
            break;
        }
        validLags = lag + 1;
    }
    if (validLags < 2) {
        return 0.0;
    }

    double integral = 0.0;
    for (std::size_t lag = 0; lag + 1 < validLags; ++lag) {
        integral += 0.5 * (u[lag] + u[lag + 1]) * dt;
    }
    return integral / 3.0;
}

}  // namespace moldyn
