#pragma once

#include <cstddef>
#include <vector>

namespace moldyn {

// Flyvbjerg-Petersen blocking: returns the mean and the standard error,
// estimated by doubling block sizes until the standard-error estimate
// plateaus (i.e. stops changing by more than its own statistical noise).
struct BlockResult {
    double mean;
    double standardError;
    std::size_t blocks;
};

BlockResult blockAverage(const std::vector<double>& samples);

}  // namespace moldyn
