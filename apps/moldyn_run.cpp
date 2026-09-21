// moldyn_run: the Phase 1 simulation driver.
//
// Two modes, selected by whether --nist-config is present:
//
//   1. NIST validation mode (--nist-config PATH): reproduces the 12-case
//      NIST Standard Reference energy/virial table for the Lennard-Jones
//      fluid and exits non-zero if any case exceeds 5e-5 relative
//      deviation. Does not simulate.
//
//   2. Simulation mode: builds an fcc lattice, equilibrates it with a
//      Langevin thermostat, runs an NVE or NVT production phase, and
//      optionally writes g(r), MSD/VACF, per-step thermodynamics and an
//      XYZ trajectory.
//
// Every CSV this program writes opens with a provenance header: the exact
// command line (reconstructed from argc/argv), the git commit the binary
// was built from (MOLDYN_GIT_SHA, injected by CMake), and a description of
// the host (OS and architecture only -- no username or hostname, since
// these files are committed to a public repository).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/utsname.h>

#include "moldyn/analyze/msd.hpp"
#include "moldyn/analyze/rdf.hpp"
#include "moldyn/analyze/vacf.hpp"
#include "moldyn/core/system.hpp"
#include "moldyn/core/vec3.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/integrate/langevin.hpp"
#include "moldyn/integrate/nose_hoover.hpp"
#include "moldyn/integrate/velocity_verlet.hpp"
#include "moldyn/io/nist_config.hpp"
#include "moldyn/neighbor/verlet_list.hpp"
#include "moldyn/util/lattice.hpp"

#ifndef MOLDYN_GIT_SHA
#define MOLDYN_GIT_SHA "unknown"
#endif

namespace {

// ---------------------------------------------------------------------
// Usage / argument parsing
// ---------------------------------------------------------------------

const char* kUsage =
    "Usage:\n"
    "  moldyn_run --nist-config PATH\n"
    "  moldyn_run --cells N --density RHO --temperature T --cutoff RC\n"
    "             --skin S --dt DT --equilibrate STEPS --production STEPS\n"
    "             --thermostat none|langevin|nose-hoover --seed S\n"
    "             [--truncation truncated|linear-force-shift]\n"
    "             [--rdf-out PATH] [--msd-out PATH] [--thermo-out PATH]\n"
    "             [--traj-out PATH]\n"
    "\n"
    "  --nist-config PATH   directory containing the NIST LJ reference\n"
    "                       configuration files (e.g. data/nist/lj), NOT a\n"
    "                       single file. Runs the 12-case validation table\n"
    "                       and exits; does not simulate.\n"
    "  --cells N            fcc cells per side (N_atoms = 4*cells^3)\n"
    "  --density RHO        reduced number density\n"
    "  --temperature T      target temperature\n"
    "  --cutoff RC          interaction cutoff\n"
    "  --skin S             Verlet skin\n"
    "  --dt DT              timestep\n"
    "  --equilibrate STEPS  steps of Langevin equilibration\n"
    "  --production STEPS   steps of NVE or NVT production\n"
    "  --thermostat KIND    none|langevin|nose-hoover (production phase)\n"
    "  --seed S             RNG seed (required for anything stochastic)\n"
    "  --truncation KIND    truncated|linear-force-shift (default: truncated).\n"
    "                       truncated: plain cutoff, energy discontinuous at\n"
    "                       rc, long-range corrected in thermo-out. Use for\n"
    "                       thermodynamics. linear-force-shift: force shifted\n"
    "                       to zero at rc, continuous, no LRC term -- isolates\n"
    "                       the integrator's own error from the cutoff\n"
    "                       discontinuity. Use to study integrator order.\n"
    "  --rdf-out PATH       write g(r) as CSV\n"
    "  --msd-out PATH       write MSD and VACF as CSV\n"
    "  --thermo-out PATH    write per-step energy, temperature, pressure as CSV\n"
    "  --traj-out PATH      write an XYZ trajectory\n";

[[noreturn]] void printUsageAndExit(const std::string& message) {
    std::cerr << "moldyn_run: " << message << "\n\n" << kUsage;
    std::exit(1);
}

const std::vector<std::string> kKnownFlags = {
    "--nist-config", "--cells",       "--density",   "--temperature", "--cutoff",
    "--skin",        "--dt",          "--equilibrate", "--production", "--thermostat",
    "--seed",        "--truncation",  "--rdf-out",   "--msd-out",   "--thermo-out",  "--traj-out",
};

std::map<std::string, std::string> parseArgs(int argc, char** argv) {
    std::map<std::string, std::string> out;
    for (int i = 1; i < argc; ++i) {
        const std::string token = argv[i];
        if (std::find(kKnownFlags.begin(), kKnownFlags.end(), token) == kKnownFlags.end()) {
            printUsageAndExit("unknown argument '" + token + "'");
        }
        if (i + 1 >= argc) {
            printUsageAndExit("missing value for '" + token + "'");
        }
        out[token] = argv[++i];
    }
    return out;
}

std::string requireFlag(const std::map<std::string, std::string>& args, const std::string& flag) {
    const auto it = args.find(flag);
    if (it == args.end()) {
        printUsageAndExit("missing required flag " + flag);
    }
    return it->second;
}

double parseDoubleArg(const std::string& flag, const std::string& value) {
    std::size_t pos = 0;
    double result = 0.0;
    try {
        result = std::stod(value, &pos);
    } catch (const std::exception&) {
        printUsageAndExit("invalid number for " + flag + ": '" + value + "'");
    }
    if (pos != value.size()) {
        printUsageAndExit("invalid number for " + flag + ": '" + value + "'");
    }
    return result;
}

long parseIntArg(const std::string& flag, const std::string& value) {
    std::size_t pos = 0;
    long result = 0;
    try {
        result = std::stol(value, &pos);
    } catch (const std::exception&) {
        printUsageAndExit("invalid integer for " + flag + ": '" + value + "'");
    }
    if (pos != value.size()) {
        printUsageAndExit("invalid integer for " + flag + ": '" + value + "'");
    }
    return result;
}

uint64_t parseSeedArg(const std::string& value) {
    if (value.empty() || value.find('-') != std::string::npos) {
        printUsageAndExit("invalid --seed '" + value + "' (must be a non-negative integer)");
    }
    std::size_t pos = 0;
    unsigned long long result = 0;
    try {
        result = std::stoull(value, &pos);
    } catch (const std::exception&) {
        printUsageAndExit("invalid --seed '" + value + "'");
    }
    if (pos != value.size()) {
        printUsageAndExit("invalid --seed '" + value + "'");
    }
    return static_cast<uint64_t>(result);
}

enum class ThermostatKind { None, Langevin, NoseHoover };

ThermostatKind parseThermostat(const std::string& value) {
    if (value == "none") return ThermostatKind::None;
    if (value == "langevin") return ThermostatKind::Langevin;
    if (value == "nose-hoover") return ThermostatKind::NoseHoover;
    printUsageAndExit("--thermostat must be one of none|langevin|nose-hoover, got '" + value + "'");
    return ThermostatKind::None;
}

moldyn::Truncation parseTruncation(const std::string& value) {
    if (value == "truncated") return moldyn::Truncation::Truncated;
    if (value == "linear-force-shift") return moldyn::Truncation::LinearForceShift;
    printUsageAndExit("--truncation must be one of truncated|linear-force-shift, got '" + value +
                       "'");
    return moldyn::Truncation::Truncated;
}

// ---------------------------------------------------------------------
// Provenance
// ---------------------------------------------------------------------

std::string gitSha() {
    const std::string sha = MOLDYN_GIT_SHA;
    return sha.empty() ? "unknown" : sha;
}

std::string hostDescription() {
    struct utsname u{};
    if (uname(&u) == 0) {
        return std::string(u.sysname) + " " + std::string(u.machine);
    }
    return "unknown";
}

std::string reconstructCommandLine(int argc, char** argv) {
    std::ostringstream oss;
    for (int i = 0; i < argc; ++i) {
        if (i != 0) oss << ' ';
        oss << argv[i];
    }
    return oss.str();
}

void writeProvenanceHeader(std::ostream& os, const std::string& commandLine) {
    os << "# command: " << commandLine << "\n";
    os << "# git_sha: " << gitSha() << "\n";
    os << "# host: " << hostDescription() << "\n";
    os << "# units: reduced Lennard-Jones units (sigma = epsilon = mass = kB = 1)\n";
}

std::ofstream openCsvOrExit(const std::string& path, const std::string& commandLine) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "moldyn_run: cannot open '" << path << "' for writing\n";
        std::exit(1);
    }
    writeProvenanceHeader(file, commandLine);
    file << std::setprecision(10);
    return file;
}

// ---------------------------------------------------------------------
// NIST validation table mode
// ---------------------------------------------------------------------

struct NistRef {
    int config;
    moldyn::Truncation trunc;
    double rc;
    double upair;
    double wpair;
    double ulrc;
};

// Same 12 reference cases as tests/test_nist_lj.cpp (NIST Standard
// Reference Simulation Website, Lennard-Jones fluid, cuboid cell).
const NistRef kNistRefs[] = {
    {1, moldyn::Truncation::Truncated, 3.0, -4.3515e3, -5.6867e2, -1.9849e2},
    {1, moldyn::Truncation::Truncated, 4.0, -4.4675e3, -1.2639e3, -8.3769e1},
    {1, moldyn::Truncation::LinearForceShift, 3.0, -3.8709e3, 3.1754e2, 0.0},
    {2, moldyn::Truncation::Truncated, 3.0, -6.9000e2, -5.6846e2, -2.4230e1},
    {2, moldyn::Truncation::Truncated, 4.0, -7.0460e2, -6.5599e2, -1.0226e1},
    {2, moldyn::Truncation::LinearForceShift, 3.0, -6.2012e2, -4.4533e2, 0.0},
    {3, moldyn::Truncation::Truncated, 3.0, -1.1467e3, -1.1649e3, -4.9622e1},
    {3, moldyn::Truncation::Truncated, 4.0, -1.1754e3, -1.3371e3, -2.0942e1},
    {3, moldyn::Truncation::LinearForceShift, 3.0, -1.0210e3, -9.3578e2, 0.0},
    {4, moldyn::Truncation::Truncated, 3.0, -1.6790e1, -4.6249e1, -5.4517e-1},
    {4, moldyn::Truncation::Truncated, 4.0, -1.7060e1, -4.7869e1, -2.3008e-1},
    {4, moldyn::Truncation::LinearForceShift, 3.0, -1.5001e1, -4.3096e1, 0.0},
};

constexpr double kNistTolerance = 5e-5;

double relativeDeviation(double calc, double ref) {
    return std::abs(calc - ref) / std::abs(ref);
}

const char* truncationLabel(moldyn::Truncation t) {
    return t == moldyn::Truncation::Truncated ? "LRC" : "LFS";
}

int runNistTable(const std::string& dir) {
    bool allPass = true;

    std::cout << std::left << std::setw(4) << "cfg" << std::setw(5) << "trunc"
              << std::right << std::setw(5) << "rc" << std::setw(14) << "u_calc"
              << std::setw(14) << "u_ref" << std::setw(12) << "u_dev" << std::setw(14)
              << "w_calc" << std::setw(14) << "w_ref" << std::setw(12) << "w_dev"
              << std::setw(14) << "ulrc_calc" << std::setw(14) << "ulrc_ref" << std::setw(12)
              << "ulrc_dev" << "  status\n";
    std::cout << std::scientific << std::setprecision(4);

    for (const NistRef& r : kNistRefs) {
        const std::string path =
            dir + "/lj_sample_config_periodic" + std::to_string(r.config) + ".txt";

        // NistConfig has no default constructor (System has none either), so
        // construct it through the optional rather than default-then-assign.
        std::optional<moldyn::NistConfig> cfgOpt;
        try {
            cfgOpt.emplace(moldyn::readNistConfig(path));
        } catch (const std::exception& e) {
            std::cerr << "moldyn_run: failed to read '" << path << "': " << e.what() << "\n";
            return 1;
        }
        const moldyn::NistConfig& cfg = *cfgOpt;

        moldyn::LennardJones lj(1.0, 1.0, r.rc, r.trunc);
        const moldyn::EnergyVirial ev = lj.computeEnergyVirial(cfg.system);
        const double lrc = lj.longRangeCorrection(cfg.system);

        const double uDev = relativeDeviation(ev.energy, r.upair);
        const double wDev = relativeDeviation(ev.virial, r.wpair);
        const bool lrcIsZero = (r.ulrc == 0.0);
        const double lrcDev = lrcIsZero ? std::abs(lrc) : relativeDeviation(lrc, r.ulrc);
        const double lrcTolerance = lrcIsZero ? 1e-9 : kNistTolerance;

        const bool pass = (uDev <= kNistTolerance) && (wDev <= kNistTolerance) &&
                           (lrcDev <= lrcTolerance);
        allPass = allPass && pass;

        std::cout << std::left << std::setw(4) << r.config << std::setw(5)
                  << truncationLabel(r.trunc) << std::right << std::setw(5) << r.rc
                  << std::setw(14) << ev.energy << std::setw(14) << r.upair << std::setw(12)
                  << uDev << std::setw(14) << ev.virial << std::setw(14) << r.wpair
                  << std::setw(12) << wDev << std::setw(14) << lrc << std::setw(14) << r.ulrc
                  << std::setw(12) << lrcDev << "  " << (pass ? "PASS" : "FAIL") << "\n";
    }

    std::cout << (allPass ? "\nALL 12 CASES WITHIN 5e-5 RELATIVE: PASS\n"
                           : "\nONE OR MORE CASES EXCEEDED 5e-5 RELATIVE: FAIL\n");
    return allPass ? 0 : 1;
}

// ---------------------------------------------------------------------
// Simulation mode
// ---------------------------------------------------------------------

struct SimArgs {
    int cells = 0;
    double density = 0.0;
    double temperature = 0.0;
    double cutoff = 0.0;
    double skin = 0.0;
    double dt = 0.0;
    long equilibrateSteps = 0;
    long productionSteps = 0;
    ThermostatKind thermostat = ThermostatKind::None;
    uint64_t seed = 0;
    moldyn::Truncation truncation = moldyn::Truncation::Truncated;
    std::optional<std::string> rdfOut;
    std::optional<std::string> msdOut;
    std::optional<std::string> thermoOut;
    std::optional<std::string> trajOut;
};

std::optional<std::string> optionalFlag(const std::map<std::string, std::string>& args,
                                         const std::string& flag) {
    const auto it = args.find(flag);
    if (it == args.end()) return std::nullopt;
    return it->second;
}

SimArgs parseSimArgs(const std::map<std::string, std::string>& args) {
    SimArgs s;

    const long cells = parseIntArg("--cells", requireFlag(args, "--cells"));
    if (cells < 1) printUsageAndExit("--cells must be a positive integer");
    s.cells = static_cast<int>(cells);

    s.density = parseDoubleArg("--density", requireFlag(args, "--density"));
    if (s.density <= 0.0) printUsageAndExit("--density must be positive");

    s.temperature = parseDoubleArg("--temperature", requireFlag(args, "--temperature"));
    if (s.temperature < 0.0) printUsageAndExit("--temperature must be non-negative");

    s.cutoff = parseDoubleArg("--cutoff", requireFlag(args, "--cutoff"));
    if (s.cutoff <= 0.0) printUsageAndExit("--cutoff must be positive");

    s.skin = parseDoubleArg("--skin", requireFlag(args, "--skin"));
    if (s.skin < 0.0) printUsageAndExit("--skin must be non-negative");

    s.dt = parseDoubleArg("--dt", requireFlag(args, "--dt"));
    if (s.dt <= 0.0) printUsageAndExit("--dt must be positive");

    s.equilibrateSteps = parseIntArg("--equilibrate", requireFlag(args, "--equilibrate"));
    if (s.equilibrateSteps < 0) printUsageAndExit("--equilibrate must be non-negative");

    s.productionSteps = parseIntArg("--production", requireFlag(args, "--production"));
    if (s.productionSteps < 0) printUsageAndExit("--production must be non-negative");

    s.thermostat = parseThermostat(requireFlag(args, "--thermostat"));
    s.seed = parseSeedArg(requireFlag(args, "--seed"));

    const std::optional<std::string> truncationArg = optionalFlag(args, "--truncation");
    if (truncationArg) s.truncation = parseTruncation(*truncationArg);

    s.rdfOut = optionalFlag(args, "--rdf-out");
    s.msdOut = optionalFlag(args, "--msd-out");
    s.thermoOut = optionalFlag(args, "--thermo-out");
    s.trajOut = optionalFlag(args, "--traj-out");

    return s;
}

// Fixed physical constants for the thermostats. These are not exposed as
// flags -- the flag list is closed by the task brief -- but are standard
// choices consistent with this project's test suite (tests/test_thermostat.cpp).
constexpr double kLangevinFriction = 1.0;
constexpr double kNoseHooverTau = 0.5;
constexpr int kNoseHooverChainLength = 3;

// RDF binning: 200 bins out to just under half the smallest box length,
// the hard limit Rdf::accumulate enforces.
constexpr std::size_t kRdfBins = 200;
constexpr double kRdfRMaxFraction = 0.98;  // of 0.5 * min box length

// Bounds the MSD/VACF ring buffer: enough lags to see the diffusive regime
// without the memory cost growing unbounded on a very long production run.
constexpr std::size_t kMaxLagCap = 2000;

// Carries everything a step needs to record: thermo row, RDF/MSD/VACF
// accumulation (production only), and trajectory frames.
class StepSampler {
public:
    StepSampler(moldyn::System& sys, const moldyn::LennardJones& lj, double& lastEnergy,
                double& lastVirial, std::size_t dof, std::ofstream* thermoFile,
                std::ofstream* trajFile, moldyn::Rdf* rdf, moldyn::Msd* msd, moldyn::Vacf* vacf)
        : sys_(sys),
          lj_(lj),
          lastEnergy_(lastEnergy),
          lastVirial_(lastVirial),
          dof_(dof),
          thermoFile_(thermoFile),
          trajFile_(trajFile),
          rdf_(rdf),
          msd_(msd),
          vacf_(vacf) {}

    void record(long step, double t, const char* phase, bool sampleAnalysis) {
        if (thermoFile_ != nullptr) {
            const double ke = sys_.kineticEnergy();
            const double temperature = sys_.temperature(dof_);
            const double volume = sys_.box().volume();
            const double numberDensity = static_cast<double>(sys_.size()) / volume;
            const double pressure = numberDensity * temperature + lastVirial_ / (3.0 * volume) +
                                     lj_.longRangeCorrectionPressure(sys_);
            const double potentialEnergy = lastEnergy_ + lj_.longRangeCorrection(sys_);
            const double totalEnergy = potentialEnergy + ke;

            (*thermoFile_) << step << ',' << t << ',' << phase << ',' << potentialEnergy << ','
                           << ke << ',' << totalEnergy << ',' << temperature << ',' << pressure
                           << '\n';
        }

        if (!sampleAnalysis) return;

        if (rdf_ != nullptr) rdf_->accumulate(sys_);
        if (msd_ != nullptr) msd_->push(sys_);
        if (vacf_ != nullptr) vacf_->push(sys_);

        if (trajFile_ != nullptr) {
            const std::size_t n = sys_.size();
            (*trajFile_) << n << '\n';
            (*trajFile_) << "step=" << step << " time=" << t << " phase=" << phase
                         << " git_sha=" << gitSha() << '\n';
            for (std::size_t atomIndex = 0; atomIndex < n; ++atomIndex) {
                const moldyn::Vec3 r = sys_.position(atomIndex);
                (*trajFile_) << "Ar " << r.x << ' ' << r.y << ' ' << r.z << '\n';
            }
        }
    }

private:
    moldyn::System& sys_;
    const moldyn::LennardJones& lj_;
    double& lastEnergy_;
    double& lastVirial_;
    std::size_t dof_;
    std::ofstream* thermoFile_;
    std::ofstream* trajFile_;
    moldyn::Rdf* rdf_;
    moldyn::Msd* msd_;
    moldyn::Vacf* vacf_;
};

int runSimulation(const SimArgs& args, const std::string& commandLine) {
    moldyn::System sys = moldyn::fccLattice(args.cells, args.density, args.temperature, args.seed);
    const std::size_t dof = 3 * sys.size() - 3;

    moldyn::LennardJones lj(1.0, 1.0, args.cutoff, args.truncation);
    moldyn::VerletList verlet(sys.box(), args.cutoff, args.skin);
    verlet.build(sys);

    double lastEnergy = 0.0;
    double lastVirial = 0.0;
    auto computeForces = [&](moldyn::System& s) {
        verlet.rebuildIfNeeded(s);
        const moldyn::EnergyVirial ev = lj.computeForces(s, verlet);
        lastEnergy = ev.energy;
        lastVirial = ev.virial;
        return ev.energy;
    };
    computeForces(sys);  // prime forces before the first integrator step

    std::optional<std::ofstream> thermoFile;
    if (args.thermoOut) {
        thermoFile = openCsvOrExit(*args.thermoOut, commandLine);
        // Both potential_energy and pressure include the analytic LJ
        // long-range (tail) correction beyond the cutoff -- see
        // LennardJones::longRangeCorrection[Pressure]. That term is
        // identically zero under linear-force-shift truncation (the
        // potential is defined to vanish at rc), so the flag below is
        // accurate either way: this line is the provenance record of the
        // fact for anyone reading this CSV later.
        const bool tailCorrected = (args.truncation == moldyn::Truncation::Truncated);
        (*thermoFile) << "# truncation: "
                       << (tailCorrected ? "truncated" : "linear-force-shift") << "\n";
        (*thermoFile) << "# pressure_includes_tail_correction: " << (tailCorrected ? "yes" : "no")
                       << "\n";
        (*thermoFile) << "# potential_energy_includes_tail_correction: "
                       << (tailCorrected ? "yes" : "no") << "\n";
        (*thermoFile) << "step,time,phase,potential_energy,kinetic_energy,total_energy,"
                         "temperature,pressure\n";
    }

    std::optional<std::ofstream> trajFile;
    if (args.trajOut) {
        trajFile.emplace(*args.trajOut);
        if (!trajFile->is_open()) {
            std::cerr << "moldyn_run: cannot open '" << *args.trajOut << "' for writing\n";
            return 1;
        }
    }

    std::optional<moldyn::Rdf> rdf;
    if (args.rdfOut) {
        const moldyn::Vec3 lengths = sys.box().lengths();
        const double minLength = std::min({lengths.x, lengths.y, lengths.z});
        const double rMax = kRdfRMaxFraction * 0.5 * minLength;
        rdf.emplace(rMax, kRdfBins);
    }

    std::size_t maxLag = 0;
    std::optional<moldyn::Msd> msd;
    std::optional<moldyn::Vacf> vacf;
    if (args.msdOut && args.productionSteps > 0) {
        maxLag = std::max<std::size_t>(
            1, std::min<std::size_t>(static_cast<std::size_t>(args.productionSteps) / 2, kMaxLagCap));
        msd.emplace(sys.size(), maxLag);
        vacf.emplace(sys.size(), maxLag);
    }

    StepSampler sampler(sys, lj, lastEnergy, lastVirial, dof, thermoFile ? &(*thermoFile) : nullptr,
                         trajFile ? &(*trajFile) : nullptr, rdf ? &(*rdf) : nullptr,
                         msd ? &(*msd) : nullptr, vacf ? &(*vacf) : nullptr);

    long globalStep = 0;
    sampler.record(globalStep, 0.0, "start", false);

    // Equilibration always uses Langevin, regardless of the production
    // thermostat, purely to thermalise the lattice start configuration
    // towards the target temperature before sampling begins.
    if (args.equilibrateSteps > 0) {
        moldyn::LangevinBAOAB equil(args.dt, kLangevinFriction, args.temperature, args.seed);
        for (long step = 0; step < args.equilibrateSteps; ++step) {
            equil.step(sys, computeForces);
            ++globalStep;
            sampler.record(globalStep, static_cast<double>(globalStep) * args.dt, "equilibrate",
                            false);
        }
    }

    // Production, with the thermostat the caller selected. Distinct seed
    // offsets keep each stochastic component's random stream independent
    // of the equilibration Langevin thermostat and of each other.
    if (args.productionSteps > 0) {
        switch (args.thermostat) {
            case ThermostatKind::None: {
                moldyn::VelocityVerlet integ(args.dt);
                for (long step = 0; step < args.productionSteps; ++step) {
                    integ.step(sys, computeForces);
                    ++globalStep;
                    sampler.record(globalStep, static_cast<double>(globalStep) * args.dt,
                                    "production", true);
                }
                break;
            }
            case ThermostatKind::Langevin: {
                moldyn::LangevinBAOAB integ(args.dt, kLangevinFriction, args.temperature,
                                             args.seed + 1);
                for (long step = 0; step < args.productionSteps; ++step) {
                    integ.step(sys, computeForces);
                    ++globalStep;
                    sampler.record(globalStep, static_cast<double>(globalStep) * args.dt,
                                    "production", true);
                }
                break;
            }
            case ThermostatKind::NoseHoover: {
                moldyn::NoseHooverChain integ(args.dt, args.temperature, kNoseHooverTau,
                                               kNoseHooverChainLength, dof);
                for (long step = 0; step < args.productionSteps; ++step) {
                    integ.step(sys, computeForces);
                    ++globalStep;
                    sampler.record(globalStep, static_cast<double>(globalStep) * args.dt,
                                    "production", true);
                }
                break;
            }
        }
    }

    if (args.rdfOut) {
        std::ofstream file = openCsvOrExit(*args.rdfOut, commandLine);
        file << "# frames: " << rdf->frames() << "\n";
        file << "r,g\n";
        const std::vector<double> centres = rdf->binCentres();
        const std::vector<double> g = rdf->g();
        for (std::size_t k = 0; k < centres.size(); ++k) {
            file << centres[k] << ',' << g[k] << '\n';
        }
    }

    if (args.msdOut && msd && vacf) {
        std::ofstream file = openCsvOrExit(*args.msdOut, commandLine);

        // Fit window for the Einstein-relation diffusion coefficient: skip
        // the ballistic short-lag regime and the noisy long-lag tail.
        std::size_t fitFrom = maxLag / 4;
        std::size_t fitTo = maxLag - maxLag / 4 + 1;
        if (fitTo <= fitFrom + 1) {
            fitFrom = 0;
            fitTo = maxLag + 1;
        }
        const double dMsd = msd->diffusionCoefficient(args.dt, fitFrom, fitTo);
        const double dVacf = vacf->diffusionCoefficient(args.dt);

        file << "# D_msd (Einstein, fit lags [" << fitFrom << "," << fitTo << ")): " << dMsd
             << "\n";
        file << "# D_vacf (Green-Kubo): " << dVacf << "\n";
        file << "lag,time,msd,vacf\n";
        const std::vector<double> msdValues = msd->msd();
        const std::vector<double> vacfValues = vacf->vacf();
        for (std::size_t lag = 0; lag < msdValues.size(); ++lag) {
            file << lag << ',' << static_cast<double>(lag) * args.dt << ',' << msdValues[lag]
                 << ',' << vacfValues[lag] << '\n';
        }
    }

    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string commandLine = reconstructCommandLine(argc, argv);
    const std::map<std::string, std::string> args = parseArgs(argc, argv);

    const auto nistIt = args.find("--nist-config");
    if (nistIt != args.end()) {
        return runNistTable(nistIt->second);
    }

    const SimArgs simArgs = parseSimArgs(args);
    return runSimulation(simArgs, commandLine);
}
