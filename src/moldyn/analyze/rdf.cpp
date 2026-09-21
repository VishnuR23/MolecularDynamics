#include "moldyn/analyze/rdf.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "moldyn/core/box.hpp"
#include "moldyn/core/vec3.hpp"

namespace moldyn {

namespace {
constexpr double kPi = 3.14159265358979323846;
}  // namespace

Rdf::Rdf(double rMax, std::size_t bins)
    : rMax_(rMax), bins_(bins), binWidth_(rMax / static_cast<double>(bins)), hist_(bins, 0.0) {}

void Rdf::accumulate(const System& sys) {
    const Box& box = sys.box();
    const Vec3 lengths = box.lengths();
    const double smallestLength = std::min({lengths.x, lengths.y, lengths.z});
    if (rMax_ > 0.5 * smallestLength) {
        throw std::invalid_argument(
            "Rdf: rMax exceeds half the smallest box length; minimum-image "
            "distances beyond that are ambiguous");
    }

    const std::size_t n = sys.size();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const Vec3 d = box.minimumImage(sys.position(i) - sys.position(j));
            const double r = norm(d);
            if (r < rMax_) {
                const std::size_t bin = static_cast<std::size_t>(r / binWidth_);
                if (bin < bins_) {
                    hist_[bin] += 1.0;
                }
            }
        }
    }

    natoms_ = n;
    sumRho_ += static_cast<double>(n) / box.volume();
    ++frames_;
}

std::vector<double> Rdf::binCentres() const {
    std::vector<double> centres(bins_);
    for (std::size_t k = 0; k < bins_; ++k) {
        centres[k] = (static_cast<double>(k) + 0.5) * binWidth_;
    }
    return centres;
}

std::vector<double> Rdf::g() const {
    std::vector<double> result(bins_, 0.0);
    if (frames_ == 0) {
        return result;
    }
    const double rho = sumRho_ / static_cast<double>(frames_);
    for (std::size_t k = 0; k < bins_; ++k) {
        const double rLo = static_cast<double>(k) * binWidth_;
        const double rHi = static_cast<double>(k + 1) * binWidth_;
        const double shellVolume = (4.0 / 3.0) * kPi * (rHi * rHi * rHi - rLo * rLo * rLo);
        const double denom = static_cast<double>(frames_) * 0.5 *
                              static_cast<double>(natoms_) * rho * shellVolume;
        result[k] = denom > 0.0 ? hist_[k] / denom : 0.0;
    }
    return result;
}

std::size_t Rdf::frames() const {
    return frames_;
}

}  // namespace moldyn
