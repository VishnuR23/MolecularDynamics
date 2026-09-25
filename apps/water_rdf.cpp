// SPC/E water structure at ambient conditions.
//
// Builds liquid water at 0.997 g/cm^3, equilibrates it at 298 K with rigid
// molecules, and accumulates the oxygen-oxygen radial distribution function.
// The published first peak sits near 2.75 Angstrom, which is what this is for.
//
// The NIST sample configurations are not at liquid density -- they are 0.374,
// 0.748, 1.122 and 0.831 g/cm^3, chosen to exercise energy routines across
// conditions rather than to be water at ambient. So the state point is built
// here: the 300-molecule configuration is expanded to the right volume by
// translating each molecule rigidly, which moves the molecules apart without
// stretching any bond.

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <sys/utsname.h>

#include "moldyn/analyze/rdf.hpp"
#include "moldyn/constrain/rigid_water.hpp"
#include "moldyn/core/box.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/units.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/forcefield/spce.hpp"
#include "moldyn/forcefield/spce_forcefield.hpp"
#include "moldyn/integrate/constrained_velocity_verlet.hpp"
#include "moldyn/io/nist_config.hpp"
#include "moldyn/util/random.hpp"

#include "moldyn_git_sha.hpp"

#ifndef MOLDYN_GIT_SHA
#define MOLDYN_GIT_SHA "unknown"
#endif

using namespace moldyn;

namespace {

constexpr double kAvogadro = 6.02214076e23;
constexpr double kWaterMolarMass = 18.01528;  // g/mol

// Box length, in Angstrom, holding `molecules` of water at `gramsPerCm3`.
double boxLengthFor(std::size_t molecules, double gramsPerCm3) {
    const double grams = static_cast<double>(molecules) * kWaterMolarMass / kAvogadro;
    const double cm3 = grams / gramsPerCm3;
    return std::cbrt(cm3 * 1.0e24);
}

// Expand or contract the box, moving each molecule as a rigid unit so no bond
// is stretched. Scaling every atom independently would deform the molecules and
// leave the constraint solver to undo it.
System rescaleRigidly(const System& in, double newLength) {
    const Vec3 old = in.box().lengths();
    const double factor = newLength / old.x;
    System out(Box(Vec3{newLength, newLength, newLength}));
    out.reserve(in.size());
    for (std::size_t m = 0; m < in.size(); m += 3) {
        Vec3 com{0.0, 0.0, 0.0};
        double mass = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            com += in.position(m + k) * in.mass(m + k);
            mass += in.mass(m + k);
        }
        com = com * (1.0 / mass);
        const Vec3 shift = com * (factor - 1.0);
        for (std::size_t k = 0; k < 3; ++k) {
            out.addAtom(in.position(m + k) + shift, Vec3{0, 0, 0}, in.mass(m + k));
            out.setCharge(out.size() - 1, in.charge(m + k));
            out.setMoleculeId(out.size() - 1, in.moleculeId(m + k));
        }
    }
    return out;
}

System gatherOxygens(const System& sys) {
    System ox(sys.box());
    ox.reserve(sys.size() / 3);
    for (std::size_t i = 0; i < sys.size(); i += 3) {
        ox.addAtom(sys.position(i), Vec3{0, 0, 0}, sys.mass(i));
    }
    return ox;
}

std::string hostDescription() {
    utsname u{};
    if (uname(&u) != 0) {
        return "unknown";
    }
    return std::string(u.sysname) + " " + u.machine;
}

}  // namespace

int main(int argc, char** argv) {
    std::string dataDir = "data";
    std::string out = "results/water_rdf.csv";
    double density = 0.997;
    double temperature = 298.0;
    double dtFs = 2.0;
    // 30 ps. The starting configuration is a NIST sample expanded from
    // 1.122 to 0.997 g/cm^3, so the system has a real density change to relax
    // through. At 10 ps it had not finished: production then ran at 308.8 K
    // against a 298 K target, the extra energy coming from potential still
    // draining into kinetic. At 30 ps it settles to 294.4 K.
    int equilibrate = 15000;
    int production = 10000;
    int sampleEvery = 10;
    std::uint64_t seed = 7001;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "water_rdf: " << a << " needs a value\n";
                std::exit(1);
            }
            return argv[++i];
        };
        if (a == "--data") dataDir = next();
        else if (a == "--out") out = next();
        else if (a == "--density") density = std::stod(next());
        else if (a == "--temperature") temperature = std::stod(next());
        else if (a == "--dt-fs") dtFs = std::stod(next());
        else if (a == "--equilibrate") equilibrate = std::stoi(next());
        else if (a == "--production") production = std::stoi(next());
        else if (a == "--sample-every") sampleEvery = std::stoi(next());
        else if (a == "--seed") seed = std::stoull(next());
        else {
            std::cerr << "water_rdf: unknown flag " << a << "\n";
            return 1;
        }
    }

    std::string commandLine;
    for (int i = 0; i < argc; ++i) {
        if (i) commandLine += " ";
        commandLine += argv[i];
    }

    const SpceParameters params;
    auto cfg = readNistConfig(dataDir + "/nist/spce/spce_sample_config_periodic3.txt");
    auto built = buildSpceSystem(cfg, params);
    const std::size_t molecules = built.moleculeCount;

    const double length = boxLengthFor(molecules, density);
    System sys = rescaleRigidly(built.full, length);

    RigidWater rigid(params.rOH, params.angleHOH);
    const double dt = units::femtoseconds(dtFs);
    ConstrainedVelocityVerlet integrator(dt, rigid);
    SpceForceField ff;

    // 6 per rigid molecule, less 3 for the removed centre-of-mass motion.
    const std::size_t dof = RigidWater::degreesOfFreedom(molecules, true);

    Rng rng(seed);
    for (std::size_t i = 0; i < sys.size(); ++i) {
        const double s = std::sqrt(temperature / sys.mass(i));
        sys.setVelocity(i, Vec3{s * rng.normal(), s * rng.normal(), s * rng.normal()});
    }
    sys.removeCenterOfMassMotion();
    rigid.constrainVelocities(sys);

    std::printf("water_rdf: %zu molecules, L = %.4f A, density %.3f g/cm3, dt = %.1f fs\n",
                molecules, length, density, dtFs);
    std::printf("water_rdf: equilibrating %d steps (%.1f ps)\n",
                equilibrate, units::toPicoseconds(dt) * equilibrate);

    auto forces = [&](System& s) { return ff.computeForces(s); };
    forces(sys);

    const auto tStart = std::chrono::steady_clock::now();

    // Equilibration: rescale velocities toward the target every 50 steps. Crude
    // by design -- it is not a thermostat that samples a correct ensemble, only
    // a way to reach the target before production, which runs at constant
    // energy and reports the temperature it actually held.
    for (int step = 0; step < equilibrate; ++step) {
        integrator.step(sys, forces);
        if ((step + 1) % 50 == 0) {
            const double t = sys.temperature(dof);
            if (t > 0.0) {
                const double scale = std::sqrt(temperature / t);
                for (std::size_t i = 0; i < sys.size(); ++i) {
                    sys.setVelocity(i, sys.velocity(i) * scale);
                }
                rigid.constrainVelocities(sys);
            }
        }
    }

    std::printf("water_rdf: production %d steps (%.1f ps), sampling every %d\n",
                production, units::toPicoseconds(dt) * production, sampleEvery);

    Rdf rdf(10.0, 200);
    double tempSum = 0.0;
    std::size_t tempCount = 0;
    for (int step = 0; step < production; ++step) {
        integrator.step(sys, forces);
        tempSum += sys.temperature(dof);
        ++tempCount;
        if ((step + 1) % sampleEvery == 0) {
            const System ox = gatherOxygens(sys);
            rdf.accumulate(ox);
        }
    }
    const auto tEnd = std::chrono::steady_clock::now();
    const double wallSeconds = std::chrono::duration<double>(tEnd - tStart).count();
    const double meanT = tempSum / static_cast<double>(tempCount);

    const std::vector<double> r = rdf.binCentres();
    const std::vector<double> g = rdf.g();

    // First peak: the tallest bin, skipping the excluded-volume region where
    // g is zero and the statistics are empty anyway.
    std::size_t peak = 0;
    double peakG = 0.0;
    for (std::size_t k = 0; k < g.size(); ++k) {
        if (r[k] < 2.0) continue;
        if (g[k] > peakG) {
            peakG = g[k];
            peak = k;
        }
    }

    // The tallest bin is only accurate to half a bin width, which at 0.05 A is
    // the same size as the difference being measured against the reference.
    // Fitting a parabola through the peak bin and its two neighbours recovers
    // the position to well inside one bin, which is what makes the comparison
    // meaningful rather than just consistent.
    const double binWidth = (r.size() > 1) ? (r[1] - r[0]) : 0.0;
    double refinedR = r[peak];
    double refinedG = peakG;
    if (peak > 0 && peak + 1 < g.size()) {
        const double gm = g[peak - 1];
        const double g0 = g[peak];
        const double gp = g[peak + 1];
        const double denom = gm - 2.0 * g0 + gp;
        if (std::abs(denom) > 1e-12) {
            const double offset = 0.5 * (gm - gp) / denom;
            refinedR = r[peak] + offset * binWidth;
            refinedG = g0 - 0.25 * (gm - gp) * offset;
        }
    }

    // g(r) must approach 1 far from any molecule. If it does not, the
    // normalisation is wrong and every peak height above is meaningless.
    double tailSum = 0.0;
    std::size_t tailCount = 0;
    for (std::size_t k = 0; k < g.size(); ++k) {
        if (r[k] > 8.5) {
            tailSum += g[k];
            ++tailCount;
        }
    }
    const double tailMean = tailCount ? tailSum / static_cast<double>(tailCount) : 0.0;

    std::ofstream f(out);
    if (!f) {
        std::cerr << "water_rdf: cannot write " << out << "\n";
        return 1;
    }
    f << "# command: " << commandLine << "\n";
    f << "# git_sha: " << MOLDYN_GIT_SHA << "\n";
    f << "# host: " << hostDescription() << "\n";
    f << "# model: SPC/E, rigid (SHAKE/RATTLE), Ewald alpha=5.6/L kmax=5 rcut=10 A\n";
    f << "# molecules: " << molecules << "\n";
    f << "# box_length_angstrom: " << length << "\n";
    f << "# density_g_per_cm3: " << density << "\n";
    f << "# target_temperature_K: " << temperature << "\n";
    f << "# mean_production_temperature_K: " << meanT << "\n";
    f << "# dt_fs: " << dtFs << "\n";
    f << "# equilibration_steps: " << equilibrate << "\n";
    f << "# production_steps: " << production << "\n";
    f << "# rdf_frames: " << rdf.frames() << "\n";
    f << "# bin_width_angstrom: " << binWidth << "\n";
    f << "# first_peak_bin_r_angstrom: " << r[peak] << "\n";
    f << "# first_peak_bin_g: " << peakG << "\n";
    f << "# first_peak_r_angstrom: " << refinedR << "\n";
    f << "# first_peak_g: " << refinedG << "\n";
    f << "# reference_first_peak_r_angstrom: 2.75\n";
    f << "# reference_first_peak_g: 3.0\n";
    f << "# tail_mean_g_beyond_8.5A: " << tailMean << "\n";
    f << "# wall_seconds: " << wallSeconds << "\n";
    f << "r_angstrom,g_oo\n";
    for (std::size_t k = 0; k < g.size(); ++k) {
        f << r[k] << "," << g[k] << "\n";
    }

    std::printf("water_rdf: mean production temperature %.1f K (target %.1f)\n", meanT, temperature);
    std::printf("water_rdf: O-O first peak at r = %.4f A, g = %.3f (published ~2.75 A, g ~3.0)\n",
                refinedR, refinedG);
    std::printf("water_rdf:   (tallest bin %.3f A, bin width %.3f A; parabolic refinement above)\n",
                r[peak], binWidth);
    std::printf("water_rdf: g(r) tail beyond 8.5 A averages %.4f (must be 1)\n", tailMean);
    std::printf("water_rdf: wrote %s in %.1f s\n", out.c_str(), wallSeconds);
    return 0;
}
