#include "moldyn/core/box.hpp"

#include <cmath>

namespace moldyn {

namespace {
double foldComponent(double value, double length) {
    return value - length * std::round(value / length);
}
}  // namespace

Box::Box(Vec3 lengths) : lengths_(lengths) {}

Vec3 Box::minimumImage(Vec3 d) const {
    return Vec3{
        foldComponent(d.x, lengths_.x),
        foldComponent(d.y, lengths_.y),
        foldComponent(d.z, lengths_.z),
    };
}

Vec3 Box::wrap(Vec3 r) const {
    return Vec3{
        foldComponent(r.x, lengths_.x),
        foldComponent(r.y, lengths_.y),
        foldComponent(r.z, lengths_.z),
    };
}

double Box::volume() const {
    return lengths_.x * lengths_.y * lengths_.z;
}

Vec3 Box::lengths() const {
    return lengths_;
}

}  // namespace moldyn
