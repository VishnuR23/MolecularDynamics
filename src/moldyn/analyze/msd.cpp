#include "moldyn/analyze/msd.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "moldyn/core/box.hpp"
#include "moldyn/core/vec3.hpp"

namespace moldyn {

Msd::Msd(std::size_t natoms, std::size_t maxLag)
    : natoms_(natoms),
      maxLag_(maxLag),
      period_(maxLag + 1),
      prevWrapped_(natoms_),
      unwrapped_(natoms_),
      ring_(period_, std::vector<Vec3>(natoms_)),
      sumSq_(period_, 0.0),
      count_(period_, 0) {}

void Msd::push(const System& sys) {
    const Box& box = sys.box();

    for (std::size_t i = 0; i < natoms_; ++i) {
        const Vec3 wrapped = sys.position(i);
        if (nPushed_ == 0) {
            unwrapped_[i] = wrapped;
        } else {
            const Vec3 disp = box.minimumImage(wrapped - prevWrapped_[i]);
            unwrapped_[i] += disp;
        }
        prevWrapped_[i] = wrapped;
    }

    const std::size_t slot = nPushed_ % period_;
    ring_[slot] = unwrapped_;

    const std::size_t maxAvailableLag = std::min(nPushed_, maxLag_);
    for (std::size_t lag = 0; lag <= maxAvailableLag; ++lag) {
        const std::size_t originSlot = (slot + period_ - lag) % period_;
        double sq = 0.0;
        for (std::size_t i = 0; i < natoms_; ++i) {
            const Vec3 d = unwrapped_[i] - ring_[originSlot][i];
            sq += norm2(d);
        }
        sumSq_[lag] += sq;
        count_[lag] += natoms_;
    }

    ++nPushed_;
}

std::vector<double> Msd::msd() const {
    std::vector<double> result(period_, 0.0);
    for (std::size_t lag = 0; lag < period_; ++lag) {
        if (count_[lag] > 0) {
            result[lag] = sumSq_[lag] / static_cast<double>(count_[lag]);
        }
    }
    return result;
}

double Msd::diffusionCoefficient(double dt, std::size_t fitFrom, std::size_t fitTo) const {
    const std::vector<double> m = msd();
    const std::size_t upper = std::min(fitTo, m.size());

    double sumT = 0.0;
    double sumM = 0.0;
    double sumTT = 0.0;
    double sumTM = 0.0;
    std::size_t n = 0;
    for (std::size_t k = fitFrom; k < upper; ++k) {
        const double t = static_cast<double>(k) * dt;
        sumT += t;
        sumM += m[k];
        sumTT += t * t;
        sumTM += t * m[k];
        ++n;
    }
    if (n < 2) {
        return 0.0;
    }

    const double nD = static_cast<double>(n);
    const double denom = nD * sumTT - sumT * sumT;
    if (denom == 0.0) {
        return 0.0;
    }
    const double slope = (nD * sumTM - sumT * sumM) / denom;
    return slope / 6.0;
}

}  // namespace moldyn
