#include "doctest/doctest.h"
#include "moldyn/core/box.hpp"
#include <cmath>

using moldyn::Box;
using moldyn::Vec3;

TEST_CASE("minimum image returns the short displacement") {
    Box box(Vec3{10.0, 10.0, 10.0});
    // 9.0 across a 10-wide box is really -1.0
    Vec3 d = box.minimumImage(Vec3{9.0, -9.0, 0.0});
    CHECK(d.x == doctest::Approx(-1.0));
    CHECK(d.y == doctest::Approx(1.0));
    CHECK(d.z == doctest::Approx(0.0));
}

TEST_CASE("minimum image is correct across all 27 images") {
    Box box(Vec3{7.0, 11.0, 13.0});
    Vec3 base{1.3, -2.7, 4.1};
    Vec3 want = box.minimumImage(base);
    for (int ix = -1; ix <= 1; ++ix)
    for (int iy = -1; iy <= 1; ++iy)
    for (int iz = -1; iz <= 1; ++iz) {
        Vec3 shifted{base.x + ix * 7.0, base.y + iy * 11.0, base.z + iz * 13.0};
        Vec3 got = box.minimumImage(shifted);
        CHECK(got.x == doctest::Approx(want.x));
        CHECK(got.y == doctest::Approx(want.y));
        CHECK(got.z == doctest::Approx(want.z));
    }
}

TEST_CASE("minimum image never exceeds half the box") {
    Box box(Vec3{4.0, 4.0, 4.0});
    for (double t = -20.0; t <= 20.0; t += 0.013) {
        Vec3 d = box.minimumImage(Vec3{t, 0.0, 0.0});
        CHECK(std::abs(d.x) <= 2.0 + 1e-12);
    }
}

TEST_CASE("wrap folds coordinates into the primary cell") {
    Box box(Vec3{10.0, 10.0, 10.0});
    Vec3 r = box.wrap(Vec3{14.0, -14.0, 5.0});
    CHECK(r.x == doctest::Approx(4.0));
    CHECK(r.y == doctest::Approx(-4.0));
    CHECK(std::abs(r.z) <= 5.0 + 1e-12);
}

TEST_CASE("volume") {
    CHECK(Box(Vec3{2.0, 3.0, 4.0}).volume() == doctest::Approx(24.0));
}

TEST_CASE("wrap and minimumImage keep exactly +-L/2 inside the box") {
    // Asymmetric box lengths so a mistake on one axis cannot hide behind
    // the others matching by coincidence.
    Box box(Vec3{10.0, 7.0, 4.0});
    const double halves[3] = {5.0, 3.5, 2.0};

    auto component = [](Vec3 v, int axis) {
        return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
    };

    for (int axis = 0; axis < 3; ++axis) {
        const double half = halves[axis];
        const double L = 2.0 * half;

        // Build a vector with -L/2 on this axis, 0 elsewhere.
        Vec3 negHalf{0.0, 0.0, 0.0};
        Vec3 posHalf{0.0, 0.0, 0.0};
        if (axis == 0) { negHalf.x = -half; posHalf.x = half; }
        if (axis == 1) { negHalf.y = -half; posHalf.y = half; }
        if (axis == 2) { negHalf.z = -half; posHalf.z = half; }

        // wrap(-L/2) must return -L/2 on that axis, and stay in range.
        {
            Vec3 r = box.wrap(negHalf);
            double c = component(r, axis);
            CHECK(c == doctest::Approx(-half));
            CHECK((c >= -half && c < half));
        }

        // wrap(+L/2) must fold to -L/2 on that axis.
        {
            Vec3 r = box.wrap(posHalf);
            double c = component(r, axis);
            CHECK(c == doctest::Approx(-half));
            CHECK((c >= -half && c < half));
        }

        // minimumImage(-L/2) must return -L/2 on that axis, and stay in range.
        {
            Vec3 d = box.minimumImage(negHalf);
            double c = component(d, axis);
            CHECK(c == doctest::Approx(-half));
            CHECK((c >= -half && c < half));
        }

        // minimumImage(+L/2) must fold to -L/2 on that axis.
        {
            Vec3 d = box.minimumImage(posHalf);
            double c = component(d, axis);
            CHECK(c == doctest::Approx(-half));
            CHECK((c >= -half && c < half));
        }

        // A full box away in either direction must also land in range
        // after wrapping: the boundary must behave correctly after an
        // arbitrary number of wraps, not just on the first fold.
        {
            Vec3 farNeg{0.0, 0.0, 0.0};
            Vec3 farPos{0.0, 0.0, 0.0};
            if (axis == 0) { farNeg.x = -half - L; farPos.x = half + L; }
            if (axis == 1) { farNeg.y = -half - L; farPos.y = half + L; }
            if (axis == 2) { farNeg.z = -half - L; farPos.z = half + L; }

            Vec3 rNeg = box.wrap(farNeg);
            double cNeg = component(rNeg, axis);
            CHECK((cNeg >= -half && cNeg < half));

            Vec3 rPos = box.wrap(farPos);
            double cPos = component(rPos, axis);
            CHECK((cPos >= -half && cPos < half));
        }
    }
}
