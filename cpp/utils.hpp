#pragma once
#include <vector>
#include <cmath>

namespace mdn {

inline double euclidean(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double d = a[i] - b[i];
        s += d*d;
    }
    return std::sqrt(s);
}

}
