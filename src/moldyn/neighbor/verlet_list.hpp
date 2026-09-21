#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "moldyn/core/box.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/neighbor/cell_list.hpp"

namespace moldyn {

// A Verlet (skin) neighbour list built on top of CellList. The candidate
// pair list is built once at a radius of cutoff+skin and reused across MD
// steps until atoms have moved far enough that a pair could have crossed the
// true cutoff without being tracked -- at which point it is rebuilt.
//
// Rebuild criterion: rebuild once the sum of the two largest atomic
// displacements since the last build exceeds the skin. Using only the
// single largest displacement is a well-known off-by-one that lets a pair
// slip inside the cutoff unnoticed (two atoms, each moving just under the
// skin but towards each other, can together close a gap of up to 2*skin).
class VerletList {
public:
    VerletList(const Box& box, double cutoff, double skin);

    void build(const System& sys);               // also records reference positions
    bool needsRebuild(const System& sys) const;   // sum of two largest displacements > skin
    void rebuildIfNeeded(const System& sys);

    std::size_t rebuildCount() const;

    // Calls fn(i, j) exactly once for each candidate pair with i < j held in
    // the list. Callers still apply the true cutoff test themselves -- the
    // list is built at cutoff+skin, so it may offer a few pairs that are
    // farther apart than the true cutoff.
    template <class Fn>
    void forEachPair(Fn&& fn) const;

private:
    Box box_;
    double cutoff_;
    double skin_;
    CellList cells_;

    std::vector<std::pair<std::size_t, std::size_t>> pairs_;
    std::vector<Vec3> referencePositions_;
    std::size_t rebuildCount_;
};

template <class Fn>
void VerletList::forEachPair(Fn&& fn) const {
    for (const auto& p : pairs_) {
        fn(p.first, p.second);
    }
}

}  // namespace moldyn
