#include "doctest/doctest.h"
#include "moldyn/analyze/rdf.hpp"
#include "moldyn/core/box.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/util/random.hpp"

#include <cstddef>
#include <vector>

using namespace moldyn;

TEST_CASE("g(r) of an ideal gas is 1 everywhere") {
    Box box(Vec3{20, 20, 20});
    Rng rng(4242);
    Rdf rdf(8.0, 80);
    for (int frame = 0; frame < 40; ++frame) {
        System sys(box);
        for (int i = 0; i < 1500; ++i)
            sys.addAtom(Vec3{(rng.uniform() - 0.5) * 20.0,
                             (rng.uniform() - 0.5) * 20.0,
                             (rng.uniform() - 0.5) * 20.0}, Vec3{0,0,0}, 1.0);
        rdf.accumulate(sys);
    }
    auto g = rdf.g();
    auto r = rdf.binCentres();
    // skip the first few bins, where the shell volume is tiny and noise is large
    for (std::size_t k = 5; k < g.size(); ++k) {
        CAPTURE(r[k]);
        CHECK(g[k] == doctest::Approx(1.0).epsilon(0.06));
    }
}

TEST_CASE("g(r) of two atoms puts all weight in one bin") {
    Box box(Vec3{20, 20, 20});
    Rdf rdf(8.0, 80);
    System sys(box);
    sys.addAtom(Vec3{0.0, 0.0, 0.0}, Vec3{0, 0, 0}, 1.0);
    sys.addAtom(Vec3{3.05, 0.0, 0.0}, Vec3{0, 0, 0}, 1.0);
    rdf.accumulate(sys);

    auto g = rdf.g();
    auto r = rdf.binCentres();
    const double binWidth = 8.0 / 80.0;
    std::size_t expectedBin = static_cast<std::size_t>(3.05 / binWidth);

    for (std::size_t k = 0; k < g.size(); ++k) {
        CAPTURE(r[k]);
        if (k == expectedBin) {
            CHECK(g[k] > 0.0);
        } else {
            CHECK(g[k] == doctest::Approx(0.0));
        }
    }
}
