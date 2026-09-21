#include "moldyn/integrate/velocity_verlet.hpp"

namespace moldyn {

VelocityVerlet::VelocityVerlet(double dt) : dt_(dt) {}

double VelocityVerlet::dt() const {
    return dt_;
}

}  // namespace moldyn
