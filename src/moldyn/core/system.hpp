#pragma once

#include <cstddef>
#include <vector>

#include "moldyn/core/box.hpp"
#include "moldyn/core/vec3.hpp"

namespace moldyn {

class System {
public:
    explicit System(Box box);

    std::size_t size() const;
    const Box& box() const;

    // Charge defaults to 0 and molecule id defaults to the atom's own index
    // (i.e. every atom is its own molecule), so every existing
    // Lennard-Jones call site is unaffected.
    void addAtom(Vec3 position, Vec3 velocity, double mass);
    // Overload for charged/molecular systems (electrostatics, Task A).
    void addAtom(Vec3 position, Vec3 velocity, double mass, double charge, std::size_t moleculeId);
    void reserve(std::size_t n);

    Vec3 position(std::size_t i) const;
    Vec3 velocity(std::size_t i) const;
    Vec3 force(std::size_t i) const;
    double mass(std::size_t i) const;
    double charge(std::size_t i) const;
    std::size_t moleculeId(std::size_t i) const;

    void setPosition(std::size_t i, Vec3 r);
    void setVelocity(std::size_t i, Vec3 v);
    void setForce(std::size_t i, Vec3 f);
    void addForce(std::size_t i, Vec3 f);
    void setCharge(std::size_t i, double q);
    void setMoleculeId(std::size_t i, std::size_t moleculeId);

    void zeroForces();
    double kineticEnergy() const;                 // sum 0.5*m*v^2
    double temperature(std::size_t dof) const;    // 2*KE / (dof * kB), kB = 1
    Vec3 totalMomentum() const;
    void removeCenterOfMassMotion();

    // Raw contiguous access for the force kernels (structure of arrays).
    const double* xs() const; const double* ys() const; const double* zs() const;
    const double* qs() const;
    double* fxs(); double* fys(); double* fzs();

private:
    Box box_;

    std::vector<double> x_, y_, z_;
    std::vector<double> vx_, vy_, vz_;
    std::vector<double> fx_, fy_, fz_;
    std::vector<double> mass_;
    std::vector<double> charge_;
    std::vector<std::size_t> moleculeId_;
};

}  // namespace moldyn
