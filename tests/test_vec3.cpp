#include "doctest/doctest.h"
#include "moldyn/core/vec3.hpp"

using moldyn::Vec3;

TEST_CASE("Vec3 arithmetic") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{0.5, -1.0, 2.0};

    CHECK((a + b).x == doctest::Approx(1.5));
    CHECK((a - b).y == doctest::Approx(3.0));
    CHECK((a * 2.0).z == doctest::Approx(6.0));
    CHECK((a / 2.0).x == doctest::Approx(0.5));

    Vec3 c = a;
    c += b;
    CHECK(c.z == doctest::Approx(5.0));
}

TEST_CASE("Vec3 dot and norm") {
    Vec3 a{1.0, 2.0, 2.0};
    CHECK(moldyn::dot(a, a) == doctest::Approx(9.0));
    CHECK(moldyn::norm2(a) == doctest::Approx(9.0));
    CHECK(moldyn::norm(a) == doctest::Approx(3.0));
    CHECK(moldyn::dot(Vec3{1,0,0}, Vec3{0,1,0}) == doctest::Approx(0.0));
}
