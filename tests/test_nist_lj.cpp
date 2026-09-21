#include "doctest/doctest.h"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/io/nist_config.hpp"
#include <string>

using moldyn::LennardJones;
using moldyn::Truncation;

namespace {
struct Ref {
    int config; Truncation trunc; double rc;
    double upair; double wpair; double ulrc;
};

// NIST Standard Reference Simulation Website, Lennard-Jones fluid,
// cuboid cell. Values are reported to 5 significant figures.
const Ref kRefs[] = {
    {1, Truncation::Truncated,        3.0, -4.3515e3, -5.6867e2, -1.9849e2},
    {1, Truncation::Truncated,        4.0, -4.4675e3, -1.2639e3, -8.3769e1},
    {1, Truncation::LinearForceShift, 3.0, -3.8709e3,  3.1754e2,  0.0     },
    {2, Truncation::Truncated,        3.0, -6.9000e2, -5.6846e2, -2.4230e1},
    {2, Truncation::Truncated,        4.0, -7.0460e2, -6.5599e2, -1.0226e1},
    {2, Truncation::LinearForceShift, 3.0, -6.2012e2, -4.4533e2,  0.0     },
    {3, Truncation::Truncated,        3.0, -1.1467e3, -1.1649e3, -4.9622e1},
    {3, Truncation::Truncated,        4.0, -1.1754e3, -1.3371e3, -2.0942e1},
    {3, Truncation::LinearForceShift, 3.0, -1.0210e3, -9.3578e2,  0.0     },
    {4, Truncation::Truncated,        3.0, -1.6790e1, -4.6249e1, -5.4517e-1},
    {4, Truncation::Truncated,        4.0, -1.7060e1, -4.7869e1, -2.3008e-1},
    {4, Truncation::LinearForceShift, 3.0, -1.5001e1, -4.3096e1,  0.0     },
};
}  // namespace

TEST_CASE("Lennard-Jones energy and virial match the NIST reference values") {
    for (const Ref& r : kRefs) {
        CAPTURE(r.config);
        CAPTURE(r.rc);
        const std::string path = std::string(MOLDYN_DATA_DIR) +
            "/nist/lj/lj_sample_config_periodic" + std::to_string(r.config) + ".txt";
        auto cfg = moldyn::readNistConfig(path);

        LennardJones lj(1.0, 1.0, r.rc, r.trunc);
        auto ev = lj.computeEnergyVirial(cfg.system);

        // NIST publishes 5 significant figures; compare relatively.
        CHECK(ev.energy == doctest::Approx(r.upair).epsilon(5e-5));
        CHECK(ev.virial == doctest::Approx(r.wpair).epsilon(5e-5));

        const double lrc = lj.longRangeCorrection(cfg.system);
        if (r.ulrc == 0.0) {
            CHECK(lrc == doctest::Approx(0.0).epsilon(1e-12));
        } else {
            CHECK(lrc == doctest::Approx(r.ulrc).epsilon(5e-5));
        }
    }
}
