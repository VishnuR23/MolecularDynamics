#pragma once

#include <cmath>

namespace moldyn {

struct Vec3 {
    double x, y, z;
};

inline constexpr Vec3 operator+(Vec3 a, Vec3 b) {
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

inline constexpr Vec3 operator-(Vec3 a, Vec3 b) {
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

inline constexpr Vec3 operator*(Vec3 a, double s) {
    return Vec3{a.x * s, a.y * s, a.z * s};
}

inline constexpr Vec3 operator*(double s, Vec3 a) {
    return Vec3{a.x * s, a.y * s, a.z * s};
}

inline constexpr Vec3 operator/(Vec3 a, double s) {
    return Vec3{a.x / s, a.y / s, a.z / s};
}

inline constexpr Vec3& operator+=(Vec3& a, Vec3 b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

inline constexpr Vec3& operator-=(Vec3& a, Vec3 b) {
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

inline constexpr Vec3& operator*=(Vec3& a, double s) {
    a.x *= s;
    a.y *= s;
    a.z *= s;
    return a;
}

inline constexpr Vec3& operator/=(Vec3& a, double s) {
    a.x /= s;
    a.y /= s;
    a.z /= s;
    return a;
}

inline constexpr double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline constexpr double norm2(Vec3 a) {
    return dot(a, a);
}

inline double norm(Vec3 a) {
    return std::sqrt(norm2(a));
}

}  // namespace moldyn
