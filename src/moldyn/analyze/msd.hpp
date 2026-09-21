#pragma once

#include <cstddef>
#include <vector>

#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"

namespace moldyn {

// Mean-squared displacement, accumulated over all time origins for each lag.
//
// Positions are stored wrapped into the box, so a naive difference of two
// stored positions is wrong whenever the atom crossed a periodic boundary
// between frames. Msd unwraps incrementally: each push() folds the minimum-
// image displacement from the previous frame's wrapped position into a
// running unwrapped position, so it never reconstructs from frame zero in
// one step.
class Msd {
public:
    Msd(std::size_t natoms, std::size_t maxLag);

    void push(const System& sys);  // unwraps internally, using sys.box()

    std::vector<double> msd() const;  // indexed by lag, size maxLag + 1

    // Einstein relation: D = slope of msd vs 6*t, fit by least squares over
    // the lag window [fitFrom, fitTo).
    double diffusionCoefficient(double dt, std::size_t fitFrom, std::size_t fitTo) const;

private:
    std::size_t natoms_;
    std::size_t maxLag_;
    std::size_t period_;  // maxLag_ + 1: the ring buffer length

    std::vector<Vec3> prevWrapped_;        // previous frame's wrapped position, per atom
    std::vector<Vec3> unwrapped_;          // running unwrapped position, per atom
    std::vector<std::vector<Vec3>> ring_;  // ring_[slot][atom]: unwrapped position at that push
    std::vector<double> sumSq_;            // indexed by lag
    std::vector<std::size_t> count_;       // indexed by lag
    std::size_t nPushed_ = 0;
};

}  // namespace moldyn
