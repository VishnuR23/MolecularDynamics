#include "moldyn/neighbor/verlet_list.hpp"

namespace moldyn {

VerletList::VerletList(const Box& box, double cutoff, double skin)
    : box_(box), cutoff_(cutoff), skin_(skin), cells_(box, cutoff + skin), rebuildCount_(0) {}

void VerletList::build(const System& sys) {
    pairs_.clear();

    const double listCutoff2 = (cutoff_ + skin_) * (cutoff_ + skin_);
    const std::size_t n = sys.size();

    if (cells_.usable()) {
        cells_.build(sys);
        cells_.forEachPair([&](std::size_t i, std::size_t j) {
            const Vec3 d = box_.minimumImage(sys.position(i) - sys.position(j));
            if (norm2(d) < listCutoff2) {
                pairs_.emplace_back(i, j);
            }
        });
    } else {
        // Box too small to cell-decompose at cutoff+skin: fall back to a
        // brute-force scan to build the candidate list. Still O(N) per step
        // afterwards via the cached pairs_, just an O(N^2) one-time cost.
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                const Vec3 d = box_.minimumImage(sys.position(i) - sys.position(j));
                if (norm2(d) < listCutoff2) {
                    pairs_.emplace_back(i, j);
                }
            }
        }
    }

    referencePositions_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        referencePositions_[i] = sys.position(i);
    }

    ++rebuildCount_;
}

bool VerletList::needsRebuild(const System& sys) const {
    double largest = 0.0;
    double secondLargest = 0.0;

    const std::size_t n = sys.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 raw = sys.position(i) - referencePositions_[i];
        const Vec3 disp = box_.minimumImage(raw);
        const double d = norm(disp);

        if (d > largest) {
            secondLargest = largest;
            largest = d;
        } else if (d > secondLargest) {
            secondLargest = d;
        }
    }

    return (largest + secondLargest) > skin_;
}

void VerletList::rebuildIfNeeded(const System& sys) {
    if (needsRebuild(sys)) {
        build(sys);
    }
}

std::size_t VerletList::rebuildCount() const {
    return rebuildCount_;
}

}  // namespace moldyn
