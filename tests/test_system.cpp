#include "doctest/doctest.h"
#include "moldyn/core/system.hpp"

using moldyn::Box; using moldyn::System; using moldyn::Vec3;

TEST_CASE("System stores and returns atoms") {
    System sys(Box(Vec3{10, 10, 10}));
    sys.addAtom(Vec3{1, 2, 3}, Vec3{0.1, 0.2, 0.3}, 1.0);
    sys.addAtom(Vec3{-1, -2, -3}, Vec3{-0.1, -0.2, -0.3}, 2.0);

    CHECK(sys.size() == 2);
    CHECK(sys.position(1).y == doctest::Approx(-2.0));
    CHECK(sys.velocity(0).z == doctest::Approx(0.3));
    CHECK(sys.mass(1) == doctest::Approx(2.0));
}

TEST_CASE("kinetic energy and temperature") {
    System sys(Box(Vec3{10, 10, 10}));
    // one atom, m=2, v=(1,0,0)  ->  KE = 0.5*2*1 = 1.0
    sys.addAtom(Vec3{0, 0, 0}, Vec3{1, 0, 0}, 2.0);
    CHECK(sys.kineticEnergy() == doctest::Approx(1.0));
    // kB = 1, dof = 3  ->  T = 2*KE/(3*1) = 2/3
    CHECK(sys.temperature(3) == doctest::Approx(2.0 / 3.0));
}

TEST_CASE("forces accumulate and zero") {
    System sys(Box(Vec3{10, 10, 10}));
    sys.addAtom(Vec3{0, 0, 0}, Vec3{0, 0, 0}, 1.0);
    sys.addForce(0, Vec3{1, 1, 1});
    sys.addForce(0, Vec3{2, 0, -1});
    CHECK(sys.force(0).x == doctest::Approx(3.0));
    CHECK(sys.force(0).z == doctest::Approx(0.0));
    sys.zeroForces();
    CHECK(sys.force(0).x == doctest::Approx(0.0));
}

TEST_CASE("removing centre of mass motion zeroes total momentum") {
    System sys(Box(Vec3{10, 10, 10}));
    sys.addAtom(Vec3{0, 0, 0}, Vec3{1.0, 0.5, -2.0}, 1.0);
    sys.addAtom(Vec3{1, 1, 1}, Vec3{3.0, -1.5, 0.0}, 4.0);
    sys.removeCenterOfMassMotion();
    CHECK(sys.totalMomentum().x == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(sys.totalMomentum().y == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(sys.totalMomentum().z == doctest::Approx(0.0).epsilon(1e-12));
}
