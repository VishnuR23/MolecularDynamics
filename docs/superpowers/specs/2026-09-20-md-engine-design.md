# Design: A Validated Molecular Dynamics Engine

**Date:** 2026-09-20
**Repo:** github.com/VishnuR23/MolecularDynamics
**Status:** Approved, pre-implementation

---

## 1. Motivation

The repository currently contains 653 lines across one commit. Every headline claim in
its README is a stub:

- `train.py` optimises against `torch.randn` noise with a freshly random label each step.
- `gnn.py` is an MLP; it accepts `edge_index` and ignores it. There is no graph.
- `predict.py` never loads a model; it returns `-4.0 - 3.0/mean_distance`.
- The "diffusion ligand generator" selects one of three hardcoded SMILES strings by hashing noise.
- `run_mini_md()` — the only function that names molecular dynamics — is a no-op.
- `coulomb()` applies the SI constant 8.99e9 to Angstrom distances, and is never called.
- The README cites a LICENSE file that does not exist.

The stated goal is a repository that stands up to review by D. E. Shaw Research. DESR
builds Anton, Desmond, and the DES-Amber force fields; their reviewers assess integrator
accuracy, energy conservation, and sampling quality professionally. A repository that
names molecular dynamics but does not perform it is a negative signal.

This design replaces the project with a real molecular dynamics engine whose every claim
is backed by a number a reviewer can check against published literature.

## 2. Goal

Ship a molecular dynamics engine in C++ that:

1. Reproduces the NIST Standard Reference Simulation benchmarks to their published precision.
2. Reproduces published physical results for liquid argon and SPC/E water.
3. Is fast, with the performance work measured and reported honestly.
4. Regenerates every committed figure and table from a single command, under CI.

Non-goal: machine learning of any kind. Non-goal: drug discovery. Non-goal: breadth. The
project is judged on whether the physics is right, not on how many features it lists.

## 3. Why NIST reference data is the spine of this project

NIST publishes exact atomic configurations together with the energies and virials those
configurations must produce. This converts "the code runs" into "the code agrees with the
national reference standard to five significant figures", which is checkable by a reviewer
in about a minute.

### Lennard-Jones benchmark (reduced units, sigma = epsilon = m = 1)

Four configurations, three truncation schemes. Source:
`nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-reference-calculations-cuboid-cell`

| Config | N | L | Scheme | Upair* | Wpair* | ULRC* |
|---|---|---|---|---|---|---|
| 1 | 800 | 10.0 | LRC, rc*=3.0 | -4.3515E+03 | -5.6867E+02 | -1.9849E+02 |
| 1 | 800 | 10.0 | LRC, rc*=4.0 | -4.4675E+03 | -1.2639E+03 | -8.3769E+01 |
| 1 | 800 | 10.0 | LFS, rc*=3.0 | -3.8709E+03 | 3.1754E+02 | 0.0 |
| 2 | 200 | 8.0 | LRC, rc*=3.0 | -6.9000E+02 | -5.6846E+02 | -2.4230E+01 |
| 2 | 200 | 8.0 | LRC, rc*=4.0 | -7.0460E+02 | -6.5599E+02 | -1.0226E+01 |
| 2 | 200 | 8.0 | LFS, rc*=3.0 | -6.2012E+02 | -4.4533E+02 | 0.0 |
| 3 | 400 | 10.0 | LRC, rc*=3.0 | -1.1467E+03 | -1.1649E+03 | -4.9622E+01 |
| 3 | 400 | 10.0 | LRC, rc*=4.0 | -1.1754E+03 | -1.3371E+03 | -2.0942E+01 |
| 3 | 400 | 10.0 | LFS, rc*=3.0 | -1.0210E+03 | -9.3578E+02 | 0.0 |
| 4 | 30 | 8.0 | LRC, rc*=3.0 | -1.6790E+01 | -4.6249E+01 | -5.4517E-01 |
| 4 | 30 | 8.0 | LRC, rc*=4.0 | -1.7060E+01 | -4.7869E+01 | -2.3008E-01 |
| 4 | 30 | 8.0 | LFS, rc*=3.0 | -1.5001E+01 | -4.3096E+01 | 0.0 |

LRC = truncated with analytic long-range correction. LFS = linear force shift at the cutoff.

### SPC/E water benchmark (energies as E/kB in Kelvin)

Source: `nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10a-cutoff`

| Quantity | Config 1 | Config 2 | Config 3 | Config 4 |
|---|---|---|---|---|
| M molecules | 100 | 200 | 300 | 750 |
| L (A) | 20.0 | 20.0 | 20.0 | 30.0 |
| Edisp | 9.95387E+04 | 1.93712E+05 | 3.54344E+05 | 4.48593E+05 |
| ELRC | -8.23715E+02 | -3.29486E+03 | -7.41343E+03 | -1.37286E+04 |
| Ereal | -5.58889E+05 | -1.19295E+06 | -1.96297E+06 | -3.57226E+06 |
| Efourier | 6.27009E+03 | 6.03495E+03 | 5.24461E+03 | 7.58785E+03 |
| Eself | -2.84469E+06 | -5.68938E+06 | -8.53407E+06 | -1.42235E+07 |
| Eintra | 2.80999E+06 | 5.61998E+06 | 8.42998E+06 | 1.41483E+07 |
| **Etotal** | **-4.88604E+05** | **-1.06590E+06** | **-1.71488E+06** | **-3.20501E+06** |

Model and method parameters fixed by NIST:

- q = 0.42380 e on each hydrogen, -2q on oxygen
- sigma = 0.316555789 nm, epsilon/kB = 78.19743111 K, single LJ site on oxygen
- r_OH = 1 A, angle H-O-H = 109.47 degrees, rigid
- Ewald damping alpha = 5.6 / min(Lx,Ly,Lz)
- kmax = 5, including only k with k^2 < 27
- real-space cutoff rcut = 10 A; dispersion truncated at rcut with analytic LRC
- periodic, tin-foil (conducting) boundary conditions
- erfc via the Numerical Recipes ERFCC routine

Both configuration sets are committed under `data/nist/` with provenance recorded.

## 4. Architecture

```
src/moldyn/          C++ core. No Python dependency. Builds standalone.
  core/              vec3, box (minimum image), system (structure-of-arrays), units
  forcefield/        lennard_jones, ewald, spce
  neighbor/          cell_list, verlet_list
  integrate/         velocity_verlet, langevin (BAOAB), nose_hoover_chain
  constrain/         settle
  analyze/           rdf, msd, vacf, block averaging with error bars
  io/                xyz, nist_config, trajectory
tests/               C++ validation suite
bindings/            pybind11 wrapper. Thin. No physics.
experiments/         one script per published result
results/             committed figures and CSVs, regenerated by CI
data/nist/           reference configurations with provenance
```

Design constraints:

- The C++ core never depends on Python. `cmake --build` and `ctest` work with no Python present.
- Physics lives only in C++. The bindings expose it; they do not reimplement it.
- Each module answers: what it does, how it is used, what it depends on. Force fields do
  not know about integrators; integrators do not know about neighbor lists.

### Units

One internal unit system: Angstrom, picosecond, atomic mass unit, kJ/mol, elementary
charge. Lennard-Jones reduced units are not a separate code path — they are the ordinary
path with sigma = epsilon = m = 1. This keeps a single tested kernel.

### Determinism

Every stochastic component takes an explicit seed. Given a seed, a run is bit-for-bit
reproducible on the same binary. Committed result files record the git commit and host
description that produced them.

## 5. Validation suite

This is the deliverable that carries the project. Tests are C++, run under `ctest`, and
gate CI.

| Test | Assertion |
|---|---|
| NIST LJ energy and virial | All 4 configs x 3 schemes match to NIST's published 5 digits |
| NIST SPC/E Ewald | All 6 energy components match on all 4 configs |
| Force vs energy | Analytic force equals central finite difference of the energy to 1e-8 |
| Energy conservation order | Drift measured across timesteps; fitted log-log slope is 2.0 +/- 0.15 |
| Time reversibility | Forward N steps, negate velocities, N steps returns to start within tolerance |
| Cell list correctness | Forces identical to brute-force O(N^2) to machine precision |
| Verlet list correctness | Same, with skin and rebuild triggering exercised |
| Minimum image | Correct across all 27 periodic images, including exact box-edge cases |
| Equipartition | Mean kinetic energy equals (3/2)NkT within statistical error |
| Maxwell-Boltzmann | Langevin velocity distribution passes a Kolmogorov-Smirnov test |
| Nose-Hoover conserved quantity | Extended-system Hamiltonian conserved to the same order |
| SETTLE | Bond lengths and angle held to 1e-10 over a full trajectory |

The energy-conservation-order test is called out deliberately: verifying that drift scales
as dt^2 rather than merely that it is "small" catches integrator errors that no
single-timestep check detects.

## 6. Experiments

Each is one command, writes a CSV and a figure to `results/`, and states its comparison
target.

| ID | Experiment | Compared against |
|---|---|---|
| 01 | NIST reference table | Our values vs NIST, printed side by side with deviations |
| 02 | Energy conservation | Drift vs timestep, log-log, fitted slope |
| 03 | Argon structure | RDF at Rahman's 1964 state point: 864 atoms, 94.4 K, 1.374 g/cm^3 |
| 04 | Argon diffusion | D from MSD and from VACF; must agree with each other and with experiment |
| 05 | LJ equation of state | Pressure and energy across state points vs published LJ EOS data |
| 06 | Water structure | SPC/E O-O, O-H, H-H RDFs; O-O first peak near 2.75 A |
| 07 | Water diffusion | SPC/E self-diffusion at 298 K, literature range 2.2-2.8e-5 cm^2/s |
| 08 | Performance | O(N^2) vs cell list scaling, thread scaling, ns/day |

Experiment 04 computing D two independent ways is a deliberate internal consistency check:
MSD and VACF are different estimators of the same quantity, and disagreement indicates a bug.

## 7. Performance

Progression, each step measured against the last:

1. Brute-force O(N^2) baseline. Correct, slow, and the reference for every later test.
2. Linked cell list, reducing to O(N).
3. Verlet neighbor list with skin distance and a displacement-triggered rebuild.
4. Structure-of-arrays layout in the force kernel.
5. ARM NEON intrinsics in the inner loop, with a scalar fallback compiled elsewhere.
6. Thread parallelism over cells using `std::thread`.

`std::thread` is chosen over OpenMP deliberately: Apple's system clang does not ship
libomp, and a build that fails on the author's own machine is not acceptable in a
portfolio repository.

Reported as a speedup table and a scaling plot, with the hardware stated (Apple M1,
8 cores). Absolute performance claims are always accompanied by the machine.

## 8. Test framework and dependencies

- C++17, CMake >= 3.16.
- doctest, vendored as a single header. Chosen over GoogleTest via FetchContent so that
  configuring the project requires no network access.
- pybind11 for bindings, optional: `-DMOLDYN_BUILD_PYTHON=OFF` skips it entirely.
- Python side, used only for plotting and driving experiments: numpy, matplotlib, pytest.
  Managed with `uv`.

Removed entirely: torch, torchvision, torchaudio, rdkit-pypi, shap, captum, streamlit,
plotly, networkx, typer, rich, pandas, scipy.

## 9. Phasing

Each phase leaves the repository in a shippable state.

**Phase 0 — Foundation.** Delete the stub tree. CMake skeleton, doctest vendored, CI on
Linux and macOS, MIT LICENSE (the file the old README already claimed), `.gitignore`,
NIST data committed with provenance.

**Phase 1 — Lennard-Jones.** vec3, box, system, LJ kernel with all three truncation
schemes, cell and Verlet lists, velocity-Verlet, Langevin, Nose-Hoover chain, RDF, MSD,
VACF. Exit criterion: NIST LJ tests green, experiments 01-05 producing committed figures.

**Phase 2 — Electrostatics and water.** Ewald real, reciprocal, self and intramolecular
terms; SPC/E; SETTLE. Exit criterion: NIST SPC/E tests green on all four configurations,
experiments 06-07 producing committed figures.

**Phase 3 — Performance.** The progression in section 7, with experiment 08.

**Phase 4 — Presentation.** README rewritten to lead with the results table and figures.
Every claim traceable to a committed artifact.

Commits are small and frequent — roughly one per file or logical unit — so that the
history reads as engineering rather than as a single drop.

## 10. Risks

**Ewald agreement to five digits (Phase 2).** The self-energy and intramolecular
correction terms are the usual source of subtle error, and NIST specifies the `erfc`
implementation because the choice is numerically visible. Mitigation: implement and test
each of the six energy components against its NIST column independently rather than
checking only the total, so a discrepancy localises immediately. If Phase 2 stalls,
Phase 1 is independently shippable and already satisfies the project goal.

**Sampling time for transport properties (experiments 04, 07).** Diffusion coefficients
need long trajectories to converge. Mitigation: report block-averaged values with error
bars and state the trajectory length; a converged error bar is a stronger result than a
bare number.

**Overstating performance.** Mitigation: every timing states hardware, compiler, flags,
and system size. No speedup is quoted without its baseline.

## 11. Definition of done

- `cmake --build build && ctest` passes on Linux and macOS in CI.
- `make reproduce` regenerates every figure and CSV in `results/`.
- The README's first screen is a results table comparing our numbers to published values.
- No claim in the repository lacks a committed artifact supporting it.
