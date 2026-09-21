#pragma once

#include "moldyn/core/vec3.hpp"

namespace moldyn {

class Box {
public:
    explicit Box(Vec3 lengths);

    Vec3 minimumImage(Vec3 d) const;  // nearest-image displacement
    Vec3 wrap(Vec3 r) const;          // fold coordinate into [-L/2, +L/2)
    double volume() const;
    Vec3 lengths() const;

private:
    Vec3 lengths_;
};

}  // namespace moldyn
