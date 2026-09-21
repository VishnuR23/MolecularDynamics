#include "moldyn/core/system.hpp"

namespace moldyn {

System::System(Box box) : box_(box) {}

std::size_t System::size() const {
    return x_.size();
}

const Box& System::box() const {
    return box_;
}

void System::addAtom(Vec3 position, Vec3 velocity, double mass) {
    x_.push_back(position.x);
    y_.push_back(position.y);
    z_.push_back(position.z);

    vx_.push_back(velocity.x);
    vy_.push_back(velocity.y);
    vz_.push_back(velocity.z);

    fx_.push_back(0.0);
    fy_.push_back(0.0);
    fz_.push_back(0.0);

    mass_.push_back(mass);
}

void System::reserve(std::size_t n) {
    x_.reserve(n);
    y_.reserve(n);
    z_.reserve(n);

    vx_.reserve(n);
    vy_.reserve(n);
    vz_.reserve(n);

    fx_.reserve(n);
    fy_.reserve(n);
    fz_.reserve(n);

    mass_.reserve(n);
}

Vec3 System::position(std::size_t i) const {
    return Vec3{x_[i], y_[i], z_[i]};
}

Vec3 System::velocity(std::size_t i) const {
    return Vec3{vx_[i], vy_[i], vz_[i]};
}

Vec3 System::force(std::size_t i) const {
    return Vec3{fx_[i], fy_[i], fz_[i]};
}

double System::mass(std::size_t i) const {
    return mass_[i];
}

void System::setPosition(std::size_t i, Vec3 r) {
    x_[i] = r.x;
    y_[i] = r.y;
    z_[i] = r.z;
}

void System::setVelocity(std::size_t i, Vec3 v) {
    vx_[i] = v.x;
    vy_[i] = v.y;
    vz_[i] = v.z;
}

void System::setForce(std::size_t i, Vec3 f) {
    fx_[i] = f.x;
    fy_[i] = f.y;
    fz_[i] = f.z;
}

void System::addForce(std::size_t i, Vec3 f) {
    fx_[i] += f.x;
    fy_[i] += f.y;
    fz_[i] += f.z;
}

void System::zeroForces() {
    for (std::size_t i = 0; i < fx_.size(); ++i) {
        fx_[i] = 0.0;
        fy_[i] = 0.0;
        fz_[i] = 0.0;
    }
}

double System::kineticEnergy() const {
    double ke = 0.0;
    for (std::size_t i = 0; i < mass_.size(); ++i) {
        const double v2 = vx_[i] * vx_[i] + vy_[i] * vy_[i] + vz_[i] * vz_[i];
        ke += 0.5 * mass_[i] * v2;
    }
    return ke;
}

double System::temperature(std::size_t dof) const {
    return 2.0 * kineticEnergy() / static_cast<double>(dof);
}

Vec3 System::totalMomentum() const {
    Vec3 p{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < mass_.size(); ++i) {
        p.x += mass_[i] * vx_[i];
        p.y += mass_[i] * vy_[i];
        p.z += mass_[i] * vz_[i];
    }
    return p;
}

void System::removeCenterOfMassMotion() {
    double totalMass = 0.0;
    for (double m : mass_) {
        totalMass += m;
    }
    if (totalMass == 0.0) {
        return;
    }

    Vec3 p = totalMomentum();
    Vec3 vCm = p / totalMass;

    for (std::size_t i = 0; i < vx_.size(); ++i) {
        vx_[i] -= vCm.x;
        vy_[i] -= vCm.y;
        vz_[i] -= vCm.z;
    }
}

const double* System::xs() const { return x_.data(); }
const double* System::ys() const { return y_.data(); }
const double* System::zs() const { return z_.data(); }

double* System::fxs() { return fx_.data(); }
double* System::fys() { return fy_.data(); }
double* System::fzs() { return fz_.data(); }

}  // namespace moldyn
