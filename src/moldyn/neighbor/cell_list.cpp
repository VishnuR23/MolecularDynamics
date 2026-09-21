#include "moldyn/neighbor/cell_list.hpp"

#include <algorithm>
#include <cmath>

namespace moldyn {

CellList::CellList(const Box& box, double cutoff)
    : box_(box), cutoff_(cutoff), nx_(1), ny_(1), nz_(1) {
    const Vec3 lengths = box_.lengths();
    nx_ = std::max(1, static_cast<int>(std::floor(lengths.x / cutoff_)));
    ny_ = std::max(1, static_cast<int>(std::floor(lengths.y / cutoff_)));
    nz_ = std::max(1, static_cast<int>(std::floor(lengths.z / cutoff_)));
}

void CellList::build(const System& sys) {
    const std::size_t n = sys.size();
    const std::size_t ncells =
        static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_) * static_cast<std::size_t>(nz_);

    head_.assign(ncells, -1);
    next_.assign(n, -1);

    const Vec3 lengths = box_.lengths();
    const double cellWidthX = lengths.x / static_cast<double>(nx_);
    const double cellWidthY = lengths.y / static_cast<double>(ny_);
    const double cellWidthZ = lengths.z / static_cast<double>(nz_);

    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 r = box_.wrap(sys.position(i));

        int ix = static_cast<int>(std::floor((r.x + 0.5 * lengths.x) / cellWidthX));
        int iy = static_cast<int>(std::floor((r.y + 0.5 * lengths.y) / cellWidthY));
        int iz = static_cast<int>(std::floor((r.z + 0.5 * lengths.z) / cellWidthZ));
        ix = clampIndex(ix, nx_);
        iy = clampIndex(iy, ny_);
        iz = clampIndex(iz, nz_);

        const int cell = cellIndex(ix, iy, iz);
        next_[i] = head_[static_cast<std::size_t>(cell)];
        head_[static_cast<std::size_t>(cell)] = static_cast<int>(i);
    }
}

int CellList::cellsPerSide(int dim) const {
    if (dim == 0) return nx_;
    if (dim == 1) return ny_;
    return nz_;
}

bool CellList::usable() const {
    return nx_ >= 3 && ny_ >= 3 && nz_ >= 3;
}

int CellList::wrapIndex(int v, int n) {
    const int m = v % n;
    return m < 0 ? m + n : m;
}

int CellList::clampIndex(int v, int n) {
    if (v < 0) return 0;
    if (v >= n) return n - 1;
    return v;
}

int CellList::cellIndex(int ix, int iy, int iz) const {
    return ix + nx_ * (iy + ny_ * iz);
}

}  // namespace moldyn
