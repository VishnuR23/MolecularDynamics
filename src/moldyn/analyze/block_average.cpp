#include "moldyn/analyze/block_average.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace moldyn {

namespace {

double mean(const std::vector<double>& data) {
    double sum = 0.0;
    for (double x : data) {
        sum += x;
    }
    return data.empty() ? 0.0 : sum / static_cast<double>(data.size());
}

// One level of the Flyvbjerg-Petersen blocking transformation: the block
// count, the standard-error estimate at this block size, and that
// estimate's own statistical uncertainty (se / sqrt(2*(n-1))), which is
// what "plateaued" is measured against.
struct Level {
    std::size_t n;
    double se;
    double errSe;
};

}  // namespace

BlockResult blockAverage(const std::vector<double>& samples) {
    if (samples.empty()) {
        return BlockResult{0.0, 0.0, 0};
    }
    if (samples.size() == 1) {
        return BlockResult{samples[0], 0.0, 1};
    }

    const double overallMean = mean(samples);

    std::vector<double> data = samples;
    std::vector<Level> levels;

    while (data.size() >= 2) {
        const std::size_t n = data.size();
        const double m = mean(data);
        double sumSq = 0.0;
        for (double x : data) {
            sumSq += (x - m) * (x - m);
        }
        const double variance = sumSq / static_cast<double>(n - 1);
        const double se = std::sqrt(variance / static_cast<double>(n));
        const double errSe = se / std::sqrt(2.0 * static_cast<double>(n - 1));
        levels.push_back(Level{n, se, errSe});

        if (n < 4) {
            break;  // too few blocks left for the estimate to mean anything
        }
        const std::size_t newN = n / 2;
        std::vector<double> blocked(newN);
        for (std::size_t k = 0; k < newN; ++k) {
            blocked[k] = 0.5 * (data[2 * k] + data[2 * k + 1]);
        }
        data = std::move(blocked);
    }

    // Plateau: the first blocking level whose standard-error estimate no
    // longer moves by more than the previous level's own noise. Until
    // decorrelation is reached, se increases level over level (correlated
    // samples make the naive, unblocked se an underestimate); once blocks
    // are effectively independent, further blocking should leave se
    // unchanged within statistical noise.
    std::size_t plateauIdx = levels.size() - 1;
    for (std::size_t i = 1; i < levels.size(); ++i) {
        if (std::fabs(levels[i].se - levels[i - 1].se) <= levels[i - 1].errSe) {
            plateauIdx = i;
            break;
        }
    }

    return BlockResult{overallMean, levels[plateauIdx].se, levels[plateauIdx].n};
}

}  // namespace moldyn
