#include "moldyn/integrate/langevin.hpp"

namespace moldyn {

LangevinBAOAB::LangevinBAOAB(double dt, double friction, double temperature, uint64_t seed)
    : dt_(dt), friction_(friction), temperature_(temperature), rng_(seed) {}

double LangevinBAOAB::temperature() const {
    return temperature_;
}

}  // namespace moldyn
