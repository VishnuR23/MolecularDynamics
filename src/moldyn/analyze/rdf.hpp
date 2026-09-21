#pragma once

#include <cstddef>
#include <vector>

#include "moldyn/core/system.hpp"

namespace moldyn {

// Radial distribution function g(r), accumulated over frames.
//
// Normalisation, with pairs counted once (i < j):
//   g(r_k) = hist[k] / ( frames * 0.5 * N * rho * shellVolume(k) )
//   shellVolume(k) = (4/3)*pi*( r_{k+1}^3 - r_k^3 )
//   rho = N / V
class Rdf {
public:
    Rdf(double rMax, std::size_t bins);

    void accumulate(const System& sys);  // one frame

    std::vector<double> binCentres() const;
    std::vector<double> g() const;  // normalised g(r)
    std::size_t frames() const;

private:
    double rMax_;
    std::size_t bins_;
    double binWidth_;
    std::vector<double> hist_;
    std::size_t frames_ = 0;
    std::size_t natoms_ = 0;      // atoms in the most recently accumulated frame
    double sumRho_ = 0.0;         // running sum of number density, for the average
};

}  // namespace moldyn
