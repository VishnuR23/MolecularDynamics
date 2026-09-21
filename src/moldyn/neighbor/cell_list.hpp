#pragma once

#include <cstddef>
#include <vector>

#include "moldyn/core/box.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"

namespace moldyn {

// Linked-cell decomposition of a periodic box for O(N) neighbour finding.
//
// The box is divided into a grid of cells each at least `cutoff` wide. Every
// pair of atoms closer than the cutoff lies either in the same cell or in one
// of its 26 neighbouring cells, so forEachPair only has to visit a cell
// against itself plus a "forward half" of 13 neighbours to enumerate every
// candidate pair exactly once.
//
// This scheme only works under the minimum image convention when there are
// at least 3 cells per side: with fewer, a cell's neighbour set would wrap
// onto itself and pairs would be double-counted or missed. usable() reports
// that condition; callers must fall back to brute force when it is false.
class CellList {
public:
    CellList(const Box& box, double cutoff);

    void build(const System& sys);

    // Calls fn(i, j) exactly once for each candidate pair with i < j whose
    // cells are adjacent (or the same cell). Callers still apply the cutoff
    // test themselves. Only valid to call when usable() is true.
    template <class Fn>
    void forEachPair(Fn&& fn) const;

    int cellsPerSide(int dim) const;  // dim in {0, 1, 2} for x, y, z
    bool usable() const;              // false when the box is too small to divide

private:
    static int wrapIndex(int v, int n);
    static int clampIndex(int v, int n);
    int cellIndex(int ix, int iy, int iz) const;

    Box box_;
    double cutoff_;
    int nx_, ny_, nz_;

    // Head-of-chain representation: head_[c] is the index of the first atom
    // in cell c (or -1 if empty); next_[i] is the next atom in i's cell's
    // chain (or -1 at the end).
    std::vector<int> head_;
    std::vector<int> next_;
};

template <class Fn>
void CellList::forEachPair(Fn&& fn) const {
    // The "forward half" of the 26 neighbouring cells: combined with a cell
    // against itself, this set visits every unordered pair of neighbouring
    // cells exactly once.
    static constexpr int kOffsets[13][3] = {
        {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {-1, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {-1, 0, 1}, {0, 1, 1}, {0, -1, 1},
        {1, 1, 1}, {1, -1, 1}, {-1, 1, 1}, {-1, -1, 1},
    };

    for (int ix = 0; ix < nx_; ++ix) {
        for (int iy = 0; iy < ny_; ++iy) {
            for (int iz = 0; iz < nz_; ++iz) {
                const int cell = cellIndex(ix, iy, iz);

                // Pairs within the same cell: each unordered pair in the
                // chain is visited exactly once by pairing every atom with
                // the ones that follow it in the chain.
                for (int i = head_[static_cast<std::size_t>(cell)]; i != -1;
                     i = next_[static_cast<std::size_t>(i)]) {
                    for (int j = next_[static_cast<std::size_t>(i)]; j != -1;
                         j = next_[static_cast<std::size_t>(j)]) {
                        if (i < j) {
                            fn(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
                        } else {
                            fn(static_cast<std::size_t>(j), static_cast<std::size_t>(i));
                        }
                    }
                }

                // Pairs against the 13 forward neighbour cells.
                for (const auto& off : kOffsets) {
                    const int jx = wrapIndex(ix + off[0], nx_);
                    const int jy = wrapIndex(iy + off[1], ny_);
                    const int jz = wrapIndex(iz + off[2], nz_);
                    const int neighbor = cellIndex(jx, jy, jz);

                    for (int i = head_[static_cast<std::size_t>(cell)]; i != -1;
                         i = next_[static_cast<std::size_t>(i)]) {
                        for (int j = head_[static_cast<std::size_t>(neighbor)]; j != -1;
                             j = next_[static_cast<std::size_t>(j)]) {
                            if (i < j) {
                                fn(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
                            } else {
                                fn(static_cast<std::size_t>(j), static_cast<std::size_t>(i));
                            }
                        }
                    }
                }
            }
        }
    }
}

}  // namespace moldyn
