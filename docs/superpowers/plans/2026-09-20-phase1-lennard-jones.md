# Phase 0+1: Foundation and Lennard-Jones Engine — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the stub repository with a C++ molecular dynamics engine whose Lennard-Jones energies and virials match the NIST Standard Reference Simulation benchmarks to their published precision, and which reproduces published liquid-argon physics.

**Architecture:** A dependency-free C++17 core under `src/moldyn/`, built by CMake and tested by doctest under `ctest`. Force fields, neighbor lists, integrators and analysis are separate modules with no back-references: a force field knows nothing about integrators, an integrator knows nothing about neighbor lists. Python appears only as an optional pybind11 binding and as plotting scripts; no physics lives there.

**Tech Stack:** C++17, CMake >= 3.16, doctest (vendored single header), GitHub Actions. Python 3.13 with numpy + matplotlib via `uv`, for plotting only.

**Spec:** `docs/superpowers/specs/2026-09-20-md-engine-design.md`

## Global Constraints

- C++17. No compiler extensions (`CMAKE_CXX_EXTENSIONS OFF`).
- The core builds and its tests pass with **no Python present**. `-DMOLDYN_BUILD_PYTHON=OFF` must be a complete, working configuration.
- No network access at CMake configure time. doctest is vendored, not fetched.
- No OpenMP. Thread parallelism, when it arrives in Phase 3, uses `std::thread`.
- Phase 1 works entirely in Lennard-Jones reduced units: `sigma = epsilon = mass = kB = 1`. Reduced units are not a separate code path — they are the ordinary code path with those values set to 1.
- Every stochastic component takes an explicit `uint64_t seed`. No hidden global RNG, no seeding from the clock.
- Namespace is `moldyn` for everything.
- Headers use `#pragma once`.
- Floating point is `double` throughout. No `float`.
- Commits are small and frequent — one per logical unit, as listed in each task's commit steps.
- Git identity for this repo is already set to `VishnuR23 <19866703+VishnuR23@users.noreply.github.com>` with SSH signing on. Do not change it.

## Verified Physics Constants

These formulas were verified before planning by an independent implementation that reproduced all 12 NIST Lennard-Jones reference values (4 configurations x 3 truncation schemes) to their published 5 significant figures. Implement exactly these.

```
u_LJ(r) = 4*eps*[ (sig/r)^12 - (sig/r)^6 ]
w_LJ(r) = -r * du_LJ/dr = 24*eps*[ 2*(sig/r)^12 - (sig/r)^6 ]

Truncation::Truncated   (NIST calls this "LRC")
    r < rc:  u = u_LJ(r)              w = w_LJ(r)
    r >= rc: u = 0                    w = 0
    The potential is NOT shifted. The tail is accounted for separately by U_LRC.

Truncation::LinearForceShift   (NIST calls this "LFS")
    dudr_rc = -(24*eps/rc) * [ 2*(sig/rc)^12 - (sig/rc)^6 ]
    r < rc:  u = u_LJ(r) - u_LJ(rc) - (r - rc)*dudr_rc
             w = w_LJ(r) + r*dudr_rc
    U_LRC is exactly 0 for this scheme.

Long-range correction (Truncated scheme only), rho = N/V:
    U_LRC = (8/3)*pi*N*rho*eps*sig^3 * [ (1/3)*(sig/rc)^9 - (sig/rc)^3 ]

Total virial is the RAW SUM over pairs, NOT divided by 3:
    W = sum_{i<j} w(r_ij)

Pair force, from w:
    f_vec_ij = ( w(r) / r^2 ) * r_vec_ij        where r_vec_ij = r_i - r_j
    Atom i gains +f_vec_ij, atom j gains -f_vec_ij.

Pressure (for later use):
    P = rho*kB*T + W / (3*V)
```

## NIST Lennard-Jones Reference Table

Reduced units. Files are in `data/nist/lj/`. Configurations 1-4 have (N, L) = (800, 10.0), (200, 8.0), (400, 10.0), (30, 8.0).

| Config | Scheme | rc | Upair* | Wpair* | ULRC* |
|---|---|---|---|---|---|
| 1 | Truncated | 3.0 | -4.3515E+03 | -5.6867E+02 | -1.9849E+02 |
| 1 | Truncated | 4.0 | -4.4675E+03 | -1.2639E+03 | -8.3769E+01 |
| 1 | LinearForceShift | 3.0 | -3.8709E+03 | 3.1754E+02 | 0.0 |
| 2 | Truncated | 3.0 | -6.9000E+02 | -5.6846E+02 | -2.4230E+01 |
| 2 | Truncated | 4.0 | -7.0460E+02 | -6.5599E+02 | -1.0226E+01 |
| 2 | LinearForceShift | 3.0 | -6.2012E+02 | -4.4533E+02 | 0.0 |
| 3 | Truncated | 3.0 | -1.1467E+03 | -1.1649E+03 | -4.9622E+01 |
| 3 | Truncated | 4.0 | -1.1754E+03 | -1.3371E+03 | -2.0942E+01 |
| 3 | LinearForceShift | 3.0 | -1.0210E+03 | -9.3578E+02 | 0.0 |
| 4 | Truncated | 3.0 | -1.6790E+01 | -4.6249E+01 | -5.4517E-01 |
| 4 | Truncated | 4.0 | -1.7060E+01 | -4.7869E+01 | -2.3008E-01 |
| 4 | LinearForceShift | 3.0 | -1.5001E+01 | -4.3096E+01 | 0.0 |

NIST reports 5 significant figures, so the test tolerance is a **relative** tolerance of `5e-5` against each value.

## File Structure

```
CMakeLists.txt                    top level; options, warnings, subdirectories
LICENSE                           MIT
.gitignore
.github/workflows/ci.yml          ubuntu-latest + macos-latest, build + ctest
third_party/doctest/doctest.h     vendored, v2.4.11
data/nist/lj/*.txt                already committed
data/nist/README.md               provenance: URLs, download date, checksums

src/moldyn/CMakeLists.txt         defines the `moldyn` static library
src/moldyn/core/vec3.hpp          Vec3 value type, header-only
src/moldyn/core/box.hpp|.cpp      Box: minimum image, wrap, volume
src/moldyn/core/system.hpp|.cpp   System: SoA state, kinetic energy, COM removal
src/moldyn/forcefield/forcefield.hpp        abstract ForceField interface
src/moldyn/forcefield/lennard_jones.hpp|.cpp
src/moldyn/neighbor/cell_list.hpp|.cpp
src/moldyn/neighbor/verlet_list.hpp|.cpp
src/moldyn/integrate/velocity_verlet.hpp|.cpp
src/moldyn/integrate/langevin.hpp|.cpp      BAOAB splitting
src/moldyn/integrate/nose_hoover.hpp|.cpp   Nose-Hoover chain
src/moldyn/analyze/rdf.hpp|.cpp
src/moldyn/analyze/msd.hpp|.cpp
src/moldyn/analyze/vacf.hpp|.cpp
src/moldyn/analyze/block_average.hpp|.cpp
src/moldyn/io/nist_config.hpp|.cpp          reads the NIST configuration format
src/moldyn/io/xyz.hpp|.cpp                  trajectory output
src/moldyn/util/random.hpp                  seeded mt19937_64 + normal draws
src/moldyn/util/lattice.hpp|.cpp            fcc lattice builder for initial states

tests/CMakeLists.txt
tests/test_main.cpp               doctest main, compiled once
tests/test_vec3.cpp
tests/test_box.cpp
tests/test_system.cpp
tests/test_nist_config.cpp
tests/test_nist_lj.cpp            THE headline test
tests/test_forces.cpp             analytic force vs finite difference
tests/test_neighbor.cpp           cell list and Verlet list vs brute force
tests/test_integrator.cpp         energy conservation order, time reversibility
tests/test_thermostat.cpp         equipartition, Maxwell-Boltzmann, NHC invariant
tests/test_analysis.cpp           rdf, msd, vacf against analytic cases

apps/moldyn_run.cpp               the simulation driver used by experiments
experiments/*.sh, experiments/plot_*.py
results/                          committed CSVs and figures
```

---

### Task 0: Strip the stubs and stand up the build

**Files:**
- Delete: `src/moldynnet/`, `app/`, `cpp/`, `configs/`, `data/toy/`, `setup.py`, `pyproject.toml`, `Makefile`, `requirements.txt`, `tests/test_cli.py`, `tests/test_cpp_backend.py`, `tests/test_models.py`, `README.md`
- Create: `LICENSE`, `.gitignore`, `CMakeLists.txt`, `src/moldyn/CMakeLists.txt`, `tests/CMakeLists.txt`, `tests/test_main.cpp`, `third_party/doctest/doctest.h`, `data/nist/README.md`, `README.md`
- Create: `.github/workflows/ci.yml`

**Interfaces:**
- Consumes: nothing.
- Produces: a CMake project with a `moldyn` static library target (initially empty) and a `moldyn_tests` executable registered with `ctest`. Later tasks add source files to `MOLDYN_SOURCES` in `src/moldyn/CMakeLists.txt` and test files to `MOLDYN_TEST_SOURCES` in `tests/CMakeLists.txt`.

- [ ] **Step 1: Delete the stub tree**

```bash
cd /Users/vishnu/molecular
git rm -r --quiet src/moldynnet app cpp configs data/toy setup.py pyproject.toml Makefile requirements.txt tests README.md
```

- [ ] **Step 2: Commit the deletion on its own**

```bash
git commit -q -m "chore: remove stub scaffolding

The previous tree described a physics+ML drug discovery framework but
implemented none of it: training optimised against torch.randn noise,
the GNN ignored its edge input, predict returned a hardcoded formula,
and run_mini_md() was a no-op. Removing it wholesale before building
the real engine.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 3: Add the MIT LICENSE the old README already claimed**

Create `LICENSE` with the standard MIT text, copyright `2026 Vishnu Rajeev`.

- [ ] **Step 4: Add .gitignore**

```
build/
.venv/
__pycache__/
*.pyc
.DS_Store
compile_commands.json
```

- [ ] **Step 5: Vendor doctest**

```bash
mkdir -p third_party/doctest
curl -sL -o third_party/doctest/doctest.h \
  https://raw.githubusercontent.com/doctest/doctest/v2.4.11/doctest/doctest.h
test -s third_party/doctest/doctest.h && grep -q DOCTEST_VERSION_STR third_party/doctest/doctest.h && echo OK
```

- [ ] **Step 6: Write the top-level CMakeLists.txt**

```cmake
cmake_minimum_required(VERSION 3.16)
project(moldyn LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(MOLDYN_BUILD_TESTS  "Build the test suite" ON)
option(MOLDYN_BUILD_PYTHON "Build the pybind11 bindings" OFF)

add_library(moldyn_warnings INTERFACE)
target_compile_options(moldyn_warnings INTERFACE
  -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wdouble-promotion)

add_subdirectory(src/moldyn)

if(MOLDYN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
```

- [ ] **Step 7: Write src/moldyn/CMakeLists.txt**

```cmake
# Source files are appended here as each task lands.
set(MOLDYN_SOURCES
  core/box.cpp
)

add_library(moldyn STATIC ${MOLDYN_SOURCES})
target_include_directories(moldyn PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(moldyn PRIVATE moldyn_warnings)
```

Note: `core/box.cpp` does not exist yet. For this task only, replace the
`set(MOLDYN_SOURCES ...)` body with `core/placeholder.cpp` containing a single
line `namespace moldyn { }` so the library has a translation unit. Task 2
deletes that file and adds `core/box.cpp`.

- [ ] **Step 8: Write tests/CMakeLists.txt**

```cmake
# Test files are appended here as each task lands.
set(MOLDYN_TEST_SOURCES
  test_main.cpp
)

add_executable(moldyn_tests ${MOLDYN_TEST_SOURCES})
target_include_directories(moldyn_tests PRIVATE ${CMAKE_SOURCE_DIR}/third_party)
target_link_libraries(moldyn_tests PRIVATE moldyn moldyn_warnings)
target_compile_definitions(moldyn_tests PRIVATE
  MOLDYN_DATA_DIR="${CMAKE_SOURCE_DIR}/data")

add_test(NAME moldyn_tests COMMAND moldyn_tests)
```

`MOLDYN_DATA_DIR` is how tests locate `data/nist/` regardless of working directory.

- [ ] **Step 9: Write tests/test_main.cpp**

```cpp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

TEST_CASE("build system is wired up") {
    CHECK(1 + 1 == 2);
}
```

- [ ] **Step 10: Configure, build and run — verify it passes**

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```
Expected: `100% tests passed, 0 tests failed out of 1`.

- [ ] **Step 11: Commit the build system**

```bash
git add LICENSE .gitignore CMakeLists.txt src tests third_party
git commit -q -m "build: CMake project, vendored doctest, ctest wiring

C++17, warnings-as-guidance (-Wall -Wextra -Wpedantic -Wshadow
-Wconversion), no compiler extensions. doctest is vendored rather than
fetched so configuring never needs the network.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 12: Record NIST data provenance**

Create `data/nist/README.md` stating, for each of the two datasets: the source
URL, the download date (2026-09-20), the file list, and the sha256 of each file.
Generate the checksums with `shasum -a 256 data/nist/*/*.txt`. Also record the
reference value tables reproduced in the spec, and the note that the LJ
conventions (unshifted truncation, undivided virial sum) were verified against
all 12 published values.

- [ ] **Step 13: Commit the data and its provenance**

```bash
git add data/nist
git commit -q -m "data: NIST reference configurations with provenance

Lennard-Jones (4 configs) and SPC/E water (4 configs) sample
configurations from the NIST Standard Reference Simulation Website,
with source URLs, download date and sha256 checksums.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 14: Add CI**

Create `.github/workflows/ci.yml`: trigger on push and pull_request; matrix over
`ubuntu-latest` and `macos-latest`; steps are checkout, `cmake -S . -B build
-DCMAKE_BUILD_TYPE=Release`, `cmake --build build -j`, `ctest --test-dir build
--output-on-failure`. No Python step — the core must build without it.

- [ ] **Step 15: Add a minimal honest README**

One short paragraph stating the project is a molecular dynamics engine under
active construction, that validation is against NIST reference data, and that
the full results table lands in Phase 4. Do not claim anything not yet true.
This file is rewritten completely in Phase 4.

- [ ] **Step 16: Commit CI and README, then push and confirm CI is green**

```bash
git add .github README.md
git commit -q -m "ci: build and test on ubuntu and macos

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
git push
gh run watch
```

---

### Task 1: Vec3

**Files:**
- Create: `src/moldyn/core/vec3.hpp`
- Test: `tests/test_vec3.cpp`
- Modify: `tests/CMakeLists.txt` (add `test_vec3.cpp` to `MOLDYN_TEST_SOURCES`)

**Interfaces:**
- Consumes: nothing.
- Produces: `moldyn::Vec3`, an aggregate `struct Vec3 { double x, y, z; };` with free operators `+ - * /`, compound `+= -= *= /=`, and free functions `double dot(Vec3, Vec3)`, `double norm2(Vec3)`, `double norm(Vec3)`. All `constexpr` where possible, all `inline`. Used by every later task.

- [ ] **Step 1: Write the failing test**

```cpp
#include "doctest/doctest.h"
#include "moldyn/core/vec3.hpp"

using moldyn::Vec3;

TEST_CASE("Vec3 arithmetic") {
    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{0.5, -1.0, 2.0};

    CHECK((a + b).x == doctest::Approx(1.5));
    CHECK((a - b).y == doctest::Approx(3.0));
    CHECK((a * 2.0).z == doctest::Approx(6.0));
    CHECK((a / 2.0).x == doctest::Approx(0.5));

    Vec3 c = a;
    c += b;
    CHECK(c.z == doctest::Approx(5.0));
}

TEST_CASE("Vec3 dot and norm") {
    Vec3 a{1.0, 2.0, 2.0};
    CHECK(moldyn::dot(a, a) == doctest::Approx(9.0));
    CHECK(moldyn::norm2(a) == doctest::Approx(9.0));
    CHECK(moldyn::norm(a) == doctest::Approx(3.0));
    CHECK(moldyn::dot(Vec3{1,0,0}, Vec3{0,1,0}) == doctest::Approx(0.0));
}
```

- [ ] **Step 2: Add the test file to the build and verify it fails**

Add `test_vec3.cpp` to `MOLDYN_TEST_SOURCES`, then:
```bash
cmake --build build -j
```
Expected: compile error, `moldyn/core/vec3.hpp` not found.

- [ ] **Step 3: Write vec3.hpp**

Header-only. `struct Vec3 { double x, y, z; };` plus the operators and free
functions named in the Interfaces block. Keep it an aggregate so `Vec3{1,2,3}`
works. Mark everything `inline constexpr` except `norm`, which calls `std::sqrt`
and is `inline` only.

- [ ] **Step 4: Build and run — verify it passes**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
```
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/core/vec3.hpp tests/test_vec3.cpp tests/CMakeLists.txt
git commit -q -m "core: Vec3 value type

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Box and the minimum image convention

**Files:**
- Create: `src/moldyn/core/box.hpp`, `src/moldyn/core/box.cpp`
- Delete: `src/moldyn/core/placeholder.cpp`
- Test: `tests/test_box.cpp`
- Modify: `src/moldyn/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `moldyn::Vec3`.
- Produces:
```cpp
namespace moldyn {
class Box {
public:
    explicit Box(Vec3 lengths);
    Vec3 minimumImage(Vec3 d) const;  // nearest-image displacement
    Vec3 wrap(Vec3 r) const;          // fold coordinate into [-L/2, +L/2)
    double volume() const;
    Vec3 lengths() const;
};
}
```
The coordinate convention is a box centred on the origin, spanning `[-L/2, +L/2)`
in each direction. This matches the NIST configuration files.

- [ ] **Step 1: Write the failing test**

```cpp
#include "doctest/doctest.h"
#include "moldyn/core/box.hpp"
#include <cmath>

using moldyn::Box;
using moldyn::Vec3;

TEST_CASE("minimum image returns the short displacement") {
    Box box(Vec3{10.0, 10.0, 10.0});
    // 9.0 across a 10-wide box is really -1.0
    Vec3 d = box.minimumImage(Vec3{9.0, -9.0, 0.0});
    CHECK(d.x == doctest::Approx(-1.0));
    CHECK(d.y == doctest::Approx(1.0));
    CHECK(d.z == doctest::Approx(0.0));
}

TEST_CASE("minimum image is correct across all 27 images") {
    Box box(Vec3{7.0, 11.0, 13.0});
    Vec3 base{1.3, -2.7, 4.1};
    Vec3 want = box.minimumImage(base);
    for (int ix = -1; ix <= 1; ++ix)
    for (int iy = -1; iy <= 1; ++iy)
    for (int iz = -1; iz <= 1; ++iz) {
        Vec3 shifted{base.x + ix * 7.0, base.y + iy * 11.0, base.z + iz * 13.0};
        Vec3 got = box.minimumImage(shifted);
        CHECK(got.x == doctest::Approx(want.x));
        CHECK(got.y == doctest::Approx(want.y));
        CHECK(got.z == doctest::Approx(want.z));
    }
}

TEST_CASE("minimum image never exceeds half the box") {
    Box box(Vec3{4.0, 4.0, 4.0});
    for (double t = -20.0; t <= 20.0; t += 0.013) {
        Vec3 d = box.minimumImage(Vec3{t, 0.0, 0.0});
        CHECK(std::abs(d.x) <= 2.0 + 1e-12);
    }
}

TEST_CASE("wrap folds coordinates into the primary cell") {
    Box box(Vec3{10.0, 10.0, 10.0});
    Vec3 r = box.wrap(Vec3{14.0, -14.0, 5.0});
    CHECK(r.x == doctest::Approx(4.0));
    CHECK(r.y == doctest::Approx(-4.0));
    CHECK(std::abs(r.z) <= 5.0 + 1e-12);
}

TEST_CASE("volume") {
    CHECK(Box(Vec3{2.0, 3.0, 4.0}).volume() == doctest::Approx(24.0));
}
```

- [ ] **Step 2: Add to build and verify it fails**

Add `test_box.cpp` to `MOLDYN_TEST_SOURCES`; add `core/box.cpp` to
`MOLDYN_SOURCES` and remove `core/placeholder.cpp`; delete the placeholder file.
```bash
cmake -S . -B build && cmake --build build -j
```
Expected: compile error, `moldyn/core/box.hpp` not found.

- [ ] **Step 3: Implement Box**

`minimumImage` uses `d.x -= L.x * std::round(d.x / L.x)` per component.
`wrap` uses the same expression — for a box centred on the origin the two are the
same operation, but keep both names because they express different intents and
`wrap` will diverge if non-cubic boxes are added later.

- [ ] **Step 4: Build and run — verify it passes**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
```

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/core/box.hpp src/moldyn/core/box.cpp tests/test_box.cpp \
        src/moldyn/CMakeLists.txt tests/CMakeLists.txt
git rm -q --ignore-unmatch src/moldyn/core/placeholder.cpp
git commit -q -m "core: periodic Box with minimum image convention

Coordinates use a box centred on the origin spanning [-L/2, +L/2),
matching the NIST configuration file convention. Tested across all 27
periodic images.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: System state

**Files:**
- Create: `src/moldyn/core/system.hpp`, `src/moldyn/core/system.cpp`
- Test: `tests/test_system.cpp`
- Modify: `src/moldyn/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Vec3`, `Box`.
- Produces:
```cpp
namespace moldyn {
class System {
public:
    explicit System(Box box);

    std::size_t size() const;
    const Box& box() const;

    void addAtom(Vec3 position, Vec3 velocity, double mass);
    void reserve(std::size_t n);

    Vec3 position(std::size_t i) const;
    Vec3 velocity(std::size_t i) const;
    Vec3 force(std::size_t i) const;
    double mass(std::size_t i) const;

    void setPosition(std::size_t i, Vec3 r);
    void setVelocity(std::size_t i, Vec3 v);
    void setForce(std::size_t i, Vec3 f);
    void addForce(std::size_t i, Vec3 f);

    void zeroForces();
    double kineticEnergy() const;                 // sum 0.5*m*v^2
    double temperature(std::size_t dof) const;    // 2*KE / (dof * kB), kB = 1
    Vec3 totalMomentum() const;
    void removeCenterOfMassMotion();

    // Raw contiguous access for the force kernels (structure of arrays).
    const double* xs() const; const double* ys() const; const double* zs() const;
    double* fxs(); double* fys(); double* fzs();
};
}
```
Storage is structure-of-arrays: separate `std::vector<double>` for each of
x, y, z, vx, vy, vz, fx, fy, fz, mass. The `Vec3` accessors are conveniences over
that storage, not the storage itself.

- [ ] **Step 1: Write the failing test**

```cpp
#include "doctest/doctest.h"
#include "moldyn/core/system.hpp"

using moldyn::Box; using moldyn::System; using moldyn::Vec3;

TEST_CASE("System stores and returns atoms") {
    System sys(Box(Vec3{10, 10, 10}));
    sys.addAtom(Vec3{1, 2, 3}, Vec3{0.1, 0.2, 0.3}, 1.0);
    sys.addAtom(Vec3{-1, -2, -3}, Vec3{-0.1, -0.2, -0.3}, 2.0);

    CHECK(sys.size() == 2);
    CHECK(sys.position(1).y == doctest::Approx(-2.0));
    CHECK(sys.velocity(0).z == doctest::Approx(0.3));
    CHECK(sys.mass(1) == doctest::Approx(2.0));
}

TEST_CASE("kinetic energy and temperature") {
    System sys(Box(Vec3{10, 10, 10}));
    // one atom, m=2, v=(1,0,0)  ->  KE = 0.5*2*1 = 1.0
    sys.addAtom(Vec3{0, 0, 0}, Vec3{1, 0, 0}, 2.0);
    CHECK(sys.kineticEnergy() == doctest::Approx(1.0));
    // kB = 1, dof = 3  ->  T = 2*KE/(3*1) = 2/3
    CHECK(sys.temperature(3) == doctest::Approx(2.0 / 3.0));
}

TEST_CASE("forces accumulate and zero") {
    System sys(Box(Vec3{10, 10, 10}));
    sys.addAtom(Vec3{0, 0, 0}, Vec3{0, 0, 0}, 1.0);
    sys.addForce(0, Vec3{1, 1, 1});
    sys.addForce(0, Vec3{2, 0, -1});
    CHECK(sys.force(0).x == doctest::Approx(3.0));
    CHECK(sys.force(0).z == doctest::Approx(0.0));
    sys.zeroForces();
    CHECK(sys.force(0).x == doctest::Approx(0.0));
}

TEST_CASE("removing centre of mass motion zeroes total momentum") {
    System sys(Box(Vec3{10, 10, 10}));
    sys.addAtom(Vec3{0, 0, 0}, Vec3{1.0, 0.5, -2.0}, 1.0);
    sys.addAtom(Vec3{1, 1, 1}, Vec3{3.0, -1.5, 0.0}, 4.0);
    sys.removeCenterOfMassMotion();
    CHECK(sys.totalMomentum().x == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(sys.totalMomentum().y == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(sys.totalMomentum().z == doctest::Approx(0.0).epsilon(1e-12));
}
```

- [ ] **Step 2: Add to build and verify it fails**

- [ ] **Step 3: Implement System** per the Interfaces block.

- [ ] **Step 4: Build and run — verify it passes**

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/core/system.hpp src/moldyn/core/system.cpp \
        tests/test_system.cpp src/moldyn/CMakeLists.txt tests/CMakeLists.txt
git commit -q -m "core: System state in structure-of-arrays layout

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: NIST configuration reader

**Files:**
- Create: `src/moldyn/io/nist_config.hpp`, `src/moldyn/io/nist_config.cpp`
- Test: `tests/test_nist_config.cpp`
- Modify: `src/moldyn/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `System`, `Box`, `Vec3`.
- Produces:
```cpp
namespace moldyn {
// Reads the NIST Standard Reference Simulation configuration format:
//   line 1: Lx Ly Lz
//   line 2: N
//   lines 3..N+2: index x y z [type]
// All atoms are given mass 1 and zero velocity. The optional 5th column
// (atom type, present in the SPC/E files) is returned in `types`; for the
// Lennard-Jones files it is empty.
struct NistConfig {
    System system;
    std::vector<std::string> types;
};
NistConfig readNistConfig(const std::string& path);
}
```
Throws `std::runtime_error` on a missing file or a malformed line.

- [ ] **Step 1: Write the failing test**

```cpp
#include "doctest/doctest.h"
#include "moldyn/io/nist_config.hpp"
#include <string>

static std::string dataPath(const std::string& rel) {
    return std::string(MOLDYN_DATA_DIR) + "/" + rel;
}

TEST_CASE("reads NIST LJ configuration 4") {
    auto cfg = moldyn::readNistConfig(dataPath("nist/lj/lj_sample_config_periodic4.txt"));
    CHECK(cfg.system.size() == 30);
    CHECK(cfg.system.box().lengths().x == doctest::Approx(8.0));
    // first atom from the file
    CHECK(cfg.system.position(0).x == doctest::Approx(1.077169909511));
    CHECK(cfg.system.position(0).y == doctest::Approx(-1.020988125886));
    // last atom from the file
    CHECK(cfg.system.position(29).z == doctest::Approx(-1.252452130644));
    CHECK(cfg.types.empty());
}

TEST_CASE("reads NIST LJ configuration 1") {
    auto cfg = moldyn::readNistConfig(dataPath("nist/lj/lj_sample_config_periodic1.txt"));
    CHECK(cfg.system.size() == 800);
    CHECK(cfg.system.box().lengths().x == doctest::Approx(10.0));
}

TEST_CASE("reads SPC/E configuration and keeps atom types") {
    auto cfg = moldyn::readNistConfig(dataPath("nist/spce/spce_sample_config_periodic1.txt"));
    CHECK(cfg.system.size() == 300);            // 100 molecules x 3 sites
    CHECK(cfg.types.size() == 300);
    CHECK(cfg.types[0] == "O");
    CHECK(cfg.types[1] == "H");
    CHECK(cfg.types[2] == "H");
}

TEST_CASE("missing file throws") {
    CHECK_THROWS_AS(moldyn::readNistConfig("/nonexistent/path.txt"), std::runtime_error);
}
```

Note the SPC/E header line gives the *molecule* count (100) while the file
contains 3N atom lines. The reader returns atoms, so `size()` is 300. The reader
must therefore read atom lines until end of file rather than trusting line 2 —
or, equivalently, treat line 2 as a lower bound and keep reading. Read to EOF.

- [ ] **Step 2: Add to build and verify it fails**

- [ ] **Step 3: Implement the reader** per the Interfaces block. Read every
remaining whitespace-separated record after the header, stopping at end of file.

- [ ] **Step 4: Build and run — verify it passes**

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/io tests/test_nist_config.cpp src/moldyn/CMakeLists.txt tests/CMakeLists.txt
git commit -q -m "io: reader for the NIST reference configuration format

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Lennard-Jones energy and virial — the NIST validation

This is the headline task. When it is green the project has its central claim.

**Files:**
- Create: `src/moldyn/forcefield/forcefield.hpp`, `src/moldyn/forcefield/lennard_jones.hpp`, `src/moldyn/forcefield/lennard_jones.cpp`
- Test: `tests/test_nist_lj.cpp`
- Modify: `src/moldyn/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `System`, `Box`, `Vec3`.
- Produces:
```cpp
namespace moldyn {

enum class Truncation { Truncated, LinearForceShift };

struct EnergyVirial {
    double energy = 0.0;       // pair sum, no tail correction
    double virial = 0.0;       // raw sum of w(r), NOT divided by 3
};

class LennardJones {
public:
    LennardJones(double sigma, double epsilon, double cutoff, Truncation truncation);

    // Pair sum by brute force over all i<j with the minimum image convention.
    // Does not touch forces.
    EnergyVirial computeEnergyVirial(const System& sys) const;

    // Analytic long-range correction to the energy. Zero for LinearForceShift.
    double longRangeCorrection(const System& sys) const;

    double sigma() const; double epsilon() const; double cutoff() const;
    Truncation truncation() const;
};
}
```

- [ ] **Step 1: Write the failing test — the full NIST table**

```cpp
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
```

- [ ] **Step 2: Add to build and verify it fails**

```bash
cmake --build build -j
```
Expected: compile error, `moldyn/forcefield/lennard_jones.hpp` not found.

- [ ] **Step 3: Implement LennardJones**

Use exactly the formulas in the **Verified Physics Constants** section above.
Implementation notes that matter:

- Work in `r^2` and avoid `sqrt` except in the `LinearForceShift` branch, which
  needs `r` itself.
- `s6 = (sigma^2 / r^2)^3`, then `u = 4*eps*(s6*s6 - s6)` and
  `w = 24*eps*(2*s6*s6 - s6)`.
- Precompute `u_LJ(rc)` and `dudr_rc` in the constructor.
- The cutoff test is `r2 < rc*rc` — strictly less than, matching the reference.
- `longRangeCorrection` returns 0 for `LinearForceShift`.

- [ ] **Step 4: Build and run — verify all 12 cases pass**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
```
Expected: PASS. If energies match but virials are off by exactly a factor of 3,
the sum is being divided — it must not be.

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/forcefield tests/test_nist_lj.cpp \
        src/moldyn/CMakeLists.txt tests/CMakeLists.txt
git commit -q -m "forcefield: Lennard-Jones matching NIST reference values

Energy and virial reproduce all four NIST Standard Reference Simulation
configurations under all three truncation schemes (truncated with
analytic long-range correction at rc=3 and rc=4, and linear force shift
at rc=3) to the published five significant figures.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Lennard-Jones forces, checked against the energy

**Files:**
- Modify: `src/moldyn/forcefield/lennard_jones.hpp`, `src/moldyn/forcefield/lennard_jones.cpp`
- Test: `tests/test_forces.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: everything from Task 5.
- Produces: added to `LennardJones`:
```cpp
// Accumulates forces into sys and returns the same energy/virial as
// computeEnergyVirial. Calls sys.zeroForces() first.
EnergyVirial computeForces(System& sys) const;
```

- [ ] **Step 1: Write the failing test**

The strongest available check on a force routine is that it is the exact negative
gradient of the energy the same code reports. Central differences give roughly
`eps^(2/3)` accuracy, so with `h = 1e-5` expect agreement near `1e-8`.

```cpp
#include "doctest/doctest.h"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/io/nist_config.hpp"
#include <string>

using moldyn::LennardJones; using moldyn::Truncation; using moldyn::Vec3;

TEST_CASE("analytic force equals minus the gradient of the energy") {
    const std::string path = std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic4.txt";   // N = 30, small and fast

    for (Truncation trunc : {Truncation::Truncated, Truncation::LinearForceShift}) {
        auto cfg = moldyn::readNistConfig(path);
        LennardJones lj(1.0, 1.0, 3.0, trunc);
        lj.computeForces(cfg.system);

        const double h = 1e-5;
        for (std::size_t i : {std::size_t(0), std::size_t(7), std::size_t(29)}) {
            Vec3 analytic = cfg.system.force(i);
            Vec3 original = cfg.system.position(i);
            double numeric[3];
            for (int d = 0; d < 3; ++d) {
                Vec3 plus = original, minus = original;
                (d == 0 ? plus.x : d == 1 ? plus.y : plus.z)   += h;
                (d == 0 ? minus.x : d == 1 ? minus.y : minus.z) -= h;

                cfg.system.setPosition(i, plus);
                double up = lj.computeEnergyVirial(cfg.system).energy;
                cfg.system.setPosition(i, minus);
                double um = lj.computeEnergyVirial(cfg.system).energy;
                cfg.system.setPosition(i, original);

                numeric[d] = -(up - um) / (2.0 * h);
            }
            CAPTURE(i);
            CHECK(analytic.x == doctest::Approx(numeric[0]).epsilon(1e-6));
            CHECK(analytic.y == doctest::Approx(numeric[1]).epsilon(1e-6));
            CHECK(analytic.z == doctest::Approx(numeric[2]).epsilon(1e-6));
        }
    }
}

TEST_CASE("forces sum to zero (Newton's third law)") {
    auto cfg = moldyn::readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic4.txt");
    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);
    lj.computeForces(cfg.system);
    Vec3 total{0, 0, 0};
    for (std::size_t i = 0; i < cfg.system.size(); ++i) total += cfg.system.force(i);
    CHECK(total.x == doctest::Approx(0.0).epsilon(1e-10));
    CHECK(total.y == doctest::Approx(0.0).epsilon(1e-10));
    CHECK(total.z == doctest::Approx(0.0).epsilon(1e-10));
}

TEST_CASE("computeForces reports the same energy as computeEnergyVirial") {
    auto cfg = moldyn::readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic2.txt");
    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);
    auto a = lj.computeEnergyVirial(cfg.system);
    auto b = lj.computeForces(cfg.system);
    CHECK(a.energy == doctest::Approx(b.energy).epsilon(1e-12));
    CHECK(a.virial == doctest::Approx(b.virial).epsilon(1e-12));
}
```

Note: the `Truncated` scheme has a force discontinuity at the cutoff, so a
finite-difference probe that straddles `rc` would fail spuriously. With `h=1e-5`
the chance of any pair sitting within `1e-5` of `rc=3.0` is negligible for these
configurations, and the test passing confirms it. If it ever fails for that
reason, the fix is to nudge `rc` to an irrational-ish value such as `3.0001`, not
to loosen the tolerance.

- [ ] **Step 2: Add to build and verify it fails**

- [ ] **Step 3: Implement computeForces**

`f_vec = (w(r) / r^2) * r_vec` with `r_vec = r_i - r_j`; add to `i`, subtract
from `j`. Share the pair loop with `computeEnergyVirial` via a private templated
helper rather than duplicating the kernel.

- [ ] **Step 4: Build and run — verify it passes**

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/forcefield tests/test_forces.cpp tests/CMakeLists.txt
git commit -q -m "forcefield: Lennard-Jones forces, verified against the energy gradient

Analytic forces are checked against central finite differences of the
same code's energy, and against Newton's third law.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: Cell list

**Files:**
- Create: `src/moldyn/neighbor/cell_list.hpp`, `src/moldyn/neighbor/cell_list.cpp`
- Test: `tests/test_neighbor.cpp`
- Modify: `src/moldyn/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `System`, `Box`, `Vec3`.
- Produces:
```cpp
namespace moldyn {
class CellList {
public:
    // cutoff is the interaction range; cells are at least this wide.
    CellList(const Box& box, double cutoff);
    void build(const System& sys);

    // Calls fn(i, j) exactly once for each candidate pair with i < j whose
    // cells are adjacent. Callers still apply the cutoff test themselves.
    template <class Fn> void forEachPair(Fn&& fn) const;

    int cellsPerSide(int dim) const;
    bool usable() const;   // false when the box is too small to divide (< 3 cells/side)
};
}
```
When `usable()` is false the caller must fall back to the brute-force loop. A box
narrower than three cutoffs per side cannot be cell-decomposed correctly under
the minimum image convention, and silently producing wrong neighbours there is
exactly the kind of bug this project exists to not have.

- [ ] **Step 1: Write the failing test**

```cpp
#include "doctest/doctest.h"
#include "moldyn/neighbor/cell_list.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/io/nist_config.hpp"

using namespace moldyn;

TEST_CASE("cell list reproduces brute-force forces to machine precision") {
    // config 1 is N=800 in L=10 with rc=3 -> 3 cells per side, the tight case
    auto cfg = readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic1.txt");
    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);

    auto brute = lj.computeForces(cfg.system);
    std::vector<Vec3> reference;
    for (std::size_t i = 0; i < cfg.system.size(); ++i)
        reference.push_back(cfg.system.force(i));

    CellList cells(cfg.system.box(), 3.0);
    REQUIRE(cells.usable());
    cells.build(cfg.system);
    auto fast = lj.computeForces(cfg.system, cells);

    CHECK(fast.energy == doctest::Approx(brute.energy).epsilon(1e-12));
    CHECK(fast.virial == doctest::Approx(brute.virial).epsilon(1e-12));
    for (std::size_t i = 0; i < cfg.system.size(); ++i) {
        CAPTURE(i);
        CHECK(cfg.system.force(i).x == doctest::Approx(reference[i].x).epsilon(1e-12));
        CHECK(cfg.system.force(i).y == doctest::Approx(reference[i].y).epsilon(1e-12));
        CHECK(cfg.system.force(i).z == doctest::Approx(reference[i].z).epsilon(1e-12));
    }
}

TEST_CASE("cell list still matches NIST after the optimisation") {
    auto cfg = readNistConfig(std::string(MOLDYN_DATA_DIR) +
        "/nist/lj/lj_sample_config_periodic1.txt");
    LennardJones lj(1.0, 1.0, 3.0, Truncation::Truncated);
    CellList cells(cfg.system.box(), 3.0);
    cells.build(cfg.system);
    auto ev = lj.computeForces(cfg.system, cells);
    CHECK(ev.energy == doctest::Approx(-4.3515e3).epsilon(5e-5));
    CHECK(ev.virial == doctest::Approx(-5.6867e2).epsilon(5e-5));
}

TEST_CASE("cell list reports itself unusable for a too-small box") {
    CellList cells(Box(Vec3{5.0, 5.0, 5.0}), 3.0);   // only 1 cell per side
    CHECK_FALSE(cells.usable());
}
```

- [ ] **Step 2: Add to build and verify it fails**

- [ ] **Step 3: Implement CellList and the LennardJones overload**

Add `EnergyVirial computeForces(System&, const CellList&) const` to
`LennardJones`. Use the classic head-of-chain representation: a `head` array of
size `ncells` and a `next` array of size `N`. Iterate each cell against its 13
forward neighbours plus itself, which visits every pair exactly once.

- [ ] **Step 4: Build and run — verify it passes**

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/neighbor tests/test_neighbor.cpp src/moldyn/forcefield \
        src/moldyn/CMakeLists.txt tests/CMakeLists.txt
git commit -q -m "neighbor: linked cell list, verified against brute force

Forces from the cell-list path are identical to the O(N^2) path to
machine precision, and the NIST reference values still hold through it.
The list reports itself unusable rather than silently producing wrong
neighbours when the box is too small to decompose.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: Verlet neighbour list

**Files:**
- Create: `src/moldyn/neighbor/verlet_list.hpp`, `src/moldyn/neighbor/verlet_list.cpp`
- Modify: `tests/test_neighbor.cpp`, `src/moldyn/CMakeLists.txt`

**Interfaces:**
- Consumes: `CellList`, `System`.
- Produces:
```cpp
namespace moldyn {
class VerletList {
public:
    VerletList(const Box& box, double cutoff, double skin);
    void build(const System& sys);              // also records reference positions
    bool needsRebuild(const System& sys) const; // true when the two largest
                                                // displacements sum past the skin
    void rebuildIfNeeded(const System& sys);
    std::size_t rebuildCount() const;
    template <class Fn> void forEachPair(Fn&& fn) const;
};
}
```
The rebuild criterion is the standard one: rebuild when the sum of the two
largest displacements since the last build exceeds the skin. Checking only the
single largest displacement is a well-known off-by-one that lets a pair slip in.

- [ ] **Step 1: Write the failing tests** — append to `tests/test_neighbor.cpp`:
forces match brute force; a system whose atoms are displaced by less than half
the skin does **not** trigger a rebuild; a system displaced past the skin does;
and after 200 MD steps the Verlet-list forces still match a fresh brute-force
computation to `1e-12`.

- [ ] **Step 2: Build and verify the new tests fail**

- [ ] **Step 3: Implement VerletList** built on top of `CellList`, plus a
`LennardJones::computeForces(System&, const VerletList&)` overload.

- [ ] **Step 4: Build and run — verify it passes**

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/neighbor tests/test_neighbor.cpp src/moldyn/forcefield \
        src/moldyn/CMakeLists.txt
git commit -q -m "neighbor: Verlet list with skin and displacement-triggered rebuild

Rebuild triggers on the sum of the two largest displacements exceeding
the skin, not the single largest.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 9: Velocity-Verlet, and proving it is second order

**Files:**
- Create: `src/moldyn/integrate/velocity_verlet.hpp`, `src/moldyn/integrate/velocity_verlet.cpp`
- Create: `src/moldyn/util/lattice.hpp`, `src/moldyn/util/lattice.cpp`
- Create: `src/moldyn/util/random.hpp`
- Test: `tests/test_integrator.cpp`
- Modify: `src/moldyn/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `System`, `LennardJones`, `VerletList`.
- Produces:
```cpp
namespace moldyn {

// Builds N atoms on a face-centred cubic lattice filling the box at the given
// number density, with velocities drawn from Maxwell-Boltzmann at `temperature`
// and centre-of-mass motion removed. `cells` is the number of fcc unit cells per
// side, so N = 4 * cells^3.
System fccLattice(int cells, double density, double temperature, uint64_t seed);

class VelocityVerlet {
public:
    explicit VelocityVerlet(double dt);
    // One step. `computeForces` must fill sys forces and return the potential energy.
    template <class ForceFn> void step(System& sys, ForceFn&& computeForces);
    double dt() const;
};
}
```
`util/random.hpp` provides `class Rng { explicit Rng(uint64_t seed); double uniform(); double normal(); }`
wrapping `std::mt19937_64`. No global state.

- [ ] **Step 1: Write the failing test — energy conservation is SECOND ORDER**

This is the test that catches integrator bugs nothing else catches. An integrator
that is merely "approximately right" still conserves energy approximately; only
the *scaling* of the error with timestep distinguishes a correct second-order
symplectic integrator from a subtly wrong one.

```cpp
#include "doctest/doctest.h"
#include "moldyn/integrate/velocity_verlet.hpp"
#include "moldyn/forcefield/lennard_jones.hpp"
#include "moldyn/util/lattice.hpp"
#include <cmath>
#include <vector>

using namespace moldyn;

TEST_CASE("velocity-Verlet energy error scales as dt^2") {
    const std::vector<double> dts = {0.001, 0.002, 0.004, 0.008};
    const double totalTime = 1.0;          // same physical time for every dt
    std::vector<double> logDt, logErr;

    for (double dt : dts) {
        System sys = fccLattice(3, 0.80, 1.0, 12345);   // N = 108
        LennardJones lj(1.0, 1.0, 2.5, Truncation::LinearForceShift);
        VelocityVerlet integrator(dt);

        auto forces = [&](System& s) { return lj.computeForces(s).energy; };
        double pe = forces(sys);
        const double e0 = pe + sys.kineticEnergy();

        double maxDev = 0.0;
        const int steps = static_cast<int>(std::lround(totalTime / dt));
        for (int n = 0; n < steps; ++n) {
            integrator.step(sys, forces);
            pe = lj.computeForces(sys).energy;
            const double e = pe + sys.kineticEnergy();
            maxDev = std::max(maxDev, std::abs(e - e0));
        }
        logDt.push_back(std::log(dt));
        logErr.push_back(std::log(maxDev / std::abs(e0)));
    }

    // least-squares slope of log(error) against log(dt)
    double mx = 0, my = 0;
    for (std::size_t i = 0; i < dts.size(); ++i) { mx += logDt[i]; my += logErr[i]; }
    mx /= static_cast<double>(dts.size());
    my /= static_cast<double>(dts.size());
    double num = 0, den = 0;
    for (std::size_t i = 0; i < dts.size(); ++i) {
        num += (logDt[i] - mx) * (logErr[i] - my);
        den += (logDt[i] - mx) * (logDt[i] - mx);
    }
    const double slope = num / den;
    CAPTURE(slope);
    CHECK(slope > 1.85);
    CHECK(slope < 2.15);
}

TEST_CASE("velocity-Verlet is time reversible") {
    System sys = fccLattice(3, 0.80, 1.0, 999);
    LennardJones lj(1.0, 1.0, 2.5, Truncation::LinearForceShift);
    VelocityVerlet integrator(0.002);
    auto forces = [&](System& s) { return lj.computeForces(s).energy; };

    std::vector<Vec3> start;
    for (std::size_t i = 0; i < sys.size(); ++i) start.push_back(sys.position(i));

    forces(sys);
    for (int n = 0; n < 200; ++n) integrator.step(sys, forces);
    for (std::size_t i = 0; i < sys.size(); ++i)
        sys.setVelocity(i, sys.velocity(i) * -1.0);
    for (int n = 0; n < 200; ++n) integrator.step(sys, forces);

    for (std::size_t i = 0; i < sys.size(); ++i) {
        CAPTURE(i);
        CHECK(sys.position(i).x == doctest::Approx(start[i].x).epsilon(1e-8));
        CHECK(sys.position(i).y == doctest::Approx(start[i].y).epsilon(1e-8));
        CHECK(sys.position(i).z == doctest::Approx(start[i].z).epsilon(1e-8));
    }
}
```

The reversibility test uses `LinearForceShift` deliberately: a plain truncated
potential is discontinuous at the cutoff, which breaks exact reversibility when a
pair crosses `rc`.

- [ ] **Step 2: Add to build and verify it fails**

- [ ] **Step 3: Implement fccLattice, Rng and VelocityVerlet**

The step is the standard kick-drift-kick:
```
v += (dt/2) * f/m
r += dt * v
f  = computeForces(sys)
v += (dt/2) * f/m
```
Positions are wrapped into the box after the drift. The integrator assumes forces
are already current on entry, which the tests honour by calling `forces(sys)`
before the loop.

- [ ] **Step 4: Build and run — verify both tests pass**

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/integrate src/moldyn/util tests/test_integrator.cpp \
        src/moldyn/CMakeLists.txt tests/CMakeLists.txt
git commit -q -m "integrate: velocity-Verlet, with a second-order convergence test

The energy-conservation test measures the SCALING of the error with
timestep and asserts the fitted log-log slope is 2, rather than merely
asserting the drift is small. That distinguishes a correct symplectic
integrator from a subtly wrong one. Time reversibility is checked
separately.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 10: Langevin thermostat (BAOAB)

**Files:**
- Create: `src/moldyn/integrate/langevin.hpp`, `src/moldyn/integrate/langevin.cpp`
- Test: `tests/test_thermostat.cpp`
- Modify: `src/moldyn/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `System`, `Rng`.
- Produces:
```cpp
namespace moldyn {
class LangevinBAOAB {
public:
    LangevinBAOAB(double dt, double friction, double temperature, uint64_t seed);
    template <class ForceFn> void step(System& sys, ForceFn&& computeForces);
    double temperature() const;
};
}
```
BAOAB splitting, which has the smallest configurational sampling error of the
common Langevin splittings:
```
B: v += (dt/2) * f/m
A: r += (dt/2) * v
O: v = c1*v + c2*sqrt(kT/m)*xi      c1 = exp(-friction*dt), c2 = sqrt(1 - c1^2)
A: r += (dt/2) * v
   f = computeForces(sys)
B: v += (dt/2) * f/m
```

- [ ] **Step 1: Write the failing tests**

Three checks, each targeting a different way a thermostat can be wrong:

1. **Equipartition** — after equilibration, the running mean of
   `2*KE / (3N - 3)` equals the target temperature within `2%`. Uses `3N-3`
   degrees of freedom because centre-of-mass motion is removed.
2. **Maxwell-Boltzmann** — pool every velocity component from many late frames
   and run a Kolmogorov-Smirnov test against `Normal(0, sqrt(kT/m))`. Assert the
   KS statistic is below the 1% critical value `1.63 / sqrt(n)`.
3. **Temperature is actually controlled** — start the system at `T = 3.0`,
   thermostat at `T = 1.0`, and confirm it relaxes to 1.0 rather than staying hot.

Write the KS statistic inline in the test file as a small static helper: sort the
sample, compute `max |i/n - Phi(x_i)|` and `max |Phi(x_i) - (i-1)/n|`, take the
larger. `Phi` is `0.5 * erfc(-x / sqrt(2))`.

- [ ] **Step 2: Add to build and verify it fails**

- [ ] **Step 3: Implement LangevinBAOAB**

- [ ] **Step 4: Build and run — verify it passes**

Note: these are statistical tests with a fixed seed, so they are deterministic
and must pass reliably. If one is marginal, lengthen the run rather than
loosening the threshold — a thermostat that only just passes equipartition is
telling you something.

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/integrate tests/test_thermostat.cpp \
        src/moldyn/CMakeLists.txt tests/CMakeLists.txt
git commit -q -m "integrate: Langevin thermostat with BAOAB splitting

Verified three ways: equipartition, a Kolmogorov-Smirnov test of the
velocity distribution against Maxwell-Boltzmann, and relaxation to the
target temperature from a hot start.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 11: Nose-Hoover chain thermostat

**Files:**
- Create: `src/moldyn/integrate/nose_hoover.hpp`, `src/moldyn/integrate/nose_hoover.cpp`
- Modify: `tests/test_thermostat.cpp`, `src/moldyn/CMakeLists.txt`

**Interfaces:**
- Consumes: `System`.
- Produces:
```cpp
namespace moldyn {
class NoseHooverChain {
public:
    // chainLength >= 1; 3 is the usual choice. tau is the thermostat period.
    NoseHooverChain(double dt, double temperature, double tau,
                    int chainLength, std::size_t degreesOfFreedom);
    template <class ForceFn> void step(System& sys, ForceFn&& computeForces);

    // The extended-system conserved quantity. Add this to the physical total
    // energy and the sum must be conserved.
    double conservedQuantityContribution() const;
};
}
```
The conserved quantity contribution is
`sum_i p_xi_i^2 / (2 Q_i) + Nf*kT*xi_1 + kT * sum_{i>1} xi_i`.

- [ ] **Step 1: Write the failing tests** — append to `tests/test_thermostat.cpp`:
the chain reaches the target temperature; and the extended conserved quantity
`H' = PE + KE + conservedQuantityContribution()` drifts by less than `1e-3`
relative over 20000 steps at `dt = 0.002`. The second is the real test — a
Nose-Hoover implementation with a wrong term still thermostats, but its extended
Hamiltonian drifts.

- [ ] **Step 2: Build and verify the new tests fail**

- [ ] **Step 3: Implement NoseHooverChain** using the Suzuki-Yoshida decomposed
chain update applied as a half-step before and after the velocity-Verlet core.

- [ ] **Step 4: Build and run — verify it passes**

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/integrate tests/test_thermostat.cpp src/moldyn/CMakeLists.txt
git commit -q -m "integrate: Nose-Hoover chain thermostat

Tested on the extended-system conserved quantity, not only on whether
the temperature comes out right — a wrong chain term still thermostats
but fails to conserve H'.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 12: Radial distribution function

**Files:**
- Create: `src/moldyn/analyze/rdf.hpp`, `src/moldyn/analyze/rdf.cpp`
- Test: `tests/test_analysis.cpp`
- Modify: `src/moldyn/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `System`, `Box`.
- Produces:
```cpp
namespace moldyn {
class Rdf {
public:
    Rdf(double rMax, std::size_t bins);
    void accumulate(const System& sys);        // one frame
    std::vector<double> binCentres() const;
    std::vector<double> g() const;             // normalised g(r)
    std::size_t frames() const;
};
}
```
Normalisation, with pairs counted once (`i < j`):
```
g(r_k) = hist[k] / ( frames * 0.5 * N * rho * shellVolume(k) )
shellVolume(k) = (4/3)*pi*( r_{k+1}^3 - r_k^3 )
rho = N / V
```

- [ ] **Step 1: Write the failing test**

The decisive check is an **ideal gas**: for randomly placed non-interacting
points, `g(r)` must equal 1 everywhere within statistical noise. This catches
every normalisation error, which is where RDF implementations almost always go
wrong.

```cpp
TEST_CASE("g(r) of an ideal gas is 1 everywhere") {
    Box box(Vec3{20, 20, 20});
    Rng rng(4242);
    Rdf rdf(8.0, 80);
    for (int frame = 0; frame < 40; ++frame) {
        System sys(box);
        for (int i = 0; i < 1500; ++i)
            sys.addAtom(Vec3{(rng.uniform() - 0.5) * 20.0,
                             (rng.uniform() - 0.5) * 20.0,
                             (rng.uniform() - 0.5) * 20.0}, Vec3{0,0,0}, 1.0);
        rdf.accumulate(sys);
    }
    auto g = rdf.g();
    auto r = rdf.binCentres();
    // skip the first few bins, where the shell volume is tiny and noise is large
    for (std::size_t k = 5; k < g.size(); ++k) {
        CAPTURE(r[k]);
        CHECK(g[k] == doctest::Approx(1.0).epsilon(0.06));
    }
}
```

Add a second test: two atoms a fixed distance apart put all their weight in the
bin containing that distance and nowhere else.

- [ ] **Step 2: Add to build and verify it fails**

- [ ] **Step 3: Implement Rdf**

- [ ] **Step 4: Build and run — verify it passes**

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/analyze tests/test_analysis.cpp \
        src/moldyn/CMakeLists.txt tests/CMakeLists.txt
git commit -q -m "analyze: radial distribution function

Normalisation is verified against an ideal gas, whose g(r) must be flat
at 1 — the check that catches RDF normalisation errors.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 13: Mean-squared displacement, velocity autocorrelation, block averaging

**Files:**
- Create: `src/moldyn/analyze/msd.hpp|.cpp`, `src/moldyn/analyze/vacf.hpp|.cpp`, `src/moldyn/analyze/block_average.hpp|.cpp`
- Modify: `tests/test_analysis.cpp`, `src/moldyn/CMakeLists.txt`

**Interfaces:**
- Consumes: `System`.
- Produces:
```cpp
namespace moldyn {

// Accumulates unwrapped displacements. The caller pushes frames in order.
class Msd {
public:
    Msd(std::size_t natoms, std::size_t maxLag);
    void push(const System& sys, const Box& box);  // unwraps internally
    std::vector<double> msd() const;               // indexed by lag
    // Einstein relation: D = slope of msd vs 6*t over the fitted window.
    double diffusionCoefficient(double dt, std::size_t fitFrom, std::size_t fitTo) const;
};

class Vacf {
public:
    Vacf(std::size_t natoms, std::size_t maxLag);
    void push(const System& sys);
    std::vector<double> vacf() const;              // normalised to 1 at lag 0
    // Green-Kubo: D = (1/3) * integral of the unnormalised VACF.
    double diffusionCoefficient(double dt) const;
};

// Flyvbjerg-Petersen blocking: returns the mean and the standard error,
// estimated by doubling block sizes until the estimate plateaus.
struct BlockResult { double mean; double stderr_; std::size_t blocks; };
BlockResult blockAverage(const std::vector<double>& samples);
}
```

- [ ] **Step 1: Write the failing tests**

- **MSD, ballistic limit:** atoms moving at constant velocity with no forces give
  exactly `msd(t) = <v^2> t^2`. Check against the analytic value to `1e-10`.
- **MSD unwrapping:** a single atom crossing the periodic boundary many times
  must give an MSD that keeps growing, not one that resets. This is the bug that
  makes diffusion coefficients silently wrong.
- **VACF:** for free particles the normalised VACF is exactly 1 at every lag.
- **Block averaging:** for independent samples the standard error matches
  `sigma/sqrt(n)` within 15%; for a strongly correlated series it must come out
  substantially *larger* than the naive estimate.

- [ ] **Step 2: Build and verify the tests fail**

- [ ] **Step 3: Implement the three modules**

- [ ] **Step 4: Build and run — verify they pass**

- [ ] **Step 5: Commit**

```bash
git add src/moldyn/analyze tests/test_analysis.cpp src/moldyn/CMakeLists.txt
git commit -q -m "analyze: MSD, VACF and Flyvbjerg-Petersen block averaging

MSD unwrapping is tested against an atom repeatedly crossing the
periodic boundary, which is the failure that silently corrupts
diffusion coefficients.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 14: Simulation driver

**Files:**
- Create: `apps/moldyn_run.cpp`, `apps/CMakeLists.txt`
- Modify: `CMakeLists.txt` (add `add_subdirectory(apps)`)

**Interfaces:**
- Consumes: everything above.
- Produces: a `moldyn_run` executable driven entirely by command-line flags, so
  the experiment scripts are a few lines of shell. Flags:
```
--nist-config PATH     run the NIST validation table instead of a simulation
--cells N              fcc cells per side (N_atoms = 4*cells^3)
--density RHO          reduced number density
--temperature T        target temperature
--cutoff RC            interaction cutoff
--skin S               Verlet skin
--dt DT                timestep
--equilibrate STEPS    steps of Langevin equilibration
--production STEPS     steps of NVE or NVT production
--thermostat none|langevin|nose-hoover
--seed S               RNG seed (required for anything stochastic)
--rdf-out PATH         write g(r) as CSV
--msd-out PATH         write MSD and VACF as CSV
--thermo-out PATH      write per-step energy, temperature, pressure as CSV
--traj-out PATH        write an XYZ trajectory
```
Every CSV gets a header comment recording the full command line, the git commit
(`MOLDYN_GIT_SHA`, injected by CMake via `git rev-parse --short HEAD`), and the
host description.

- [ ] **Step 1: Write the driver** with the flags above, exiting non-zero with a
usage message on an unknown flag or a missing required one.

- [ ] **Step 2: Inject the git SHA in CMake**

```cmake
execute_process(COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE MOLDYN_GIT_SHA OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
target_compile_definitions(moldyn_run PRIVATE MOLDYN_GIT_SHA="${MOLDYN_GIT_SHA}")
```

- [ ] **Step 3: Verify the NIST table mode reproduces the reference values**

```bash
cmake --build build -j
./build/apps/moldyn_run --nist-config data/nist/lj
```
Expected: a printed table of all 12 cases with a deviation column, every row
within `5e-5` relative, and exit code 0.

- [ ] **Step 4: Verify a short simulation runs end to end**

```bash
./build/apps/moldyn_run --cells 5 --density 0.8442 --temperature 0.728 \
  --cutoff 2.5 --skin 0.3 --dt 0.004 --equilibrate 2000 --production 5000 \
  --thermostat langevin --seed 1 --thermo-out /tmp/thermo.csv --rdf-out /tmp/gr.csv
head -5 /tmp/gr.csv
```
Expected: exit 0, CSVs with provenance headers and plausible data.

- [ ] **Step 5: Commit**

```bash
git add apps CMakeLists.txt
git commit -q -m "apps: moldyn_run simulation driver

Every output CSV records the full command line, the git commit that
produced it, and the host, so a committed result can always be traced
back to the code that made it.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 15: Experiments 01-05 and their figures

**Files:**
- Create: `experiments/run_all.sh`, `experiments/exp01_nist_table.sh`, `experiments/exp02_energy_conservation.sh`, `experiments/exp03_argon_rdf.sh`, `experiments/exp04_argon_diffusion.sh`, `experiments/exp05_lj_eos.sh`
- Create: `experiments/plot_common.py`, `experiments/plot_exp02.py` … `plot_exp05.py`
- Create: `experiments/requirements.txt` (`numpy`, `matplotlib` only)
- Create: `results/` outputs
- Create: `Makefile` with a `reproduce` target

**Interfaces:**
- Consumes: `moldyn_run`.
- Produces: `results/*.csv` and `results/*.png`, all committed.

State points, fixed here so the experiments are reproducible:

| Experiment | Setup | Compared against |
|---|---|---|
| 01 | `--nist-config data/nist/lj` | the 12 NIST values, printed with deviations |
| 02 | fcc 3 cells, rho*=0.80, T*=1.0, NVE, dt in {0.001,0.002,0.004,0.008}, t=1.0 | fitted log-log slope, expected 2 |
| 03 | 864 atoms (fcc 6 cells), rho*=0.8442, T*=0.728 — Rahman's argon state point in reduced units | published argon g(r): first peak near r*=1.1, height near 3.0 |
| 04 | same state point, NVE production after Langevin equilibration | D from MSD and from VACF must agree within their error bars |
| 05 | rho* in {0.5,0.6,0.7,0.8,0.9}, T* in {1.0,1.5,2.0} | published LJ equation-of-state pressures and energies |

- [ ] **Step 1: Set up the Python environment**

```bash
cd /Users/vishnu/molecular
uv venv
uv pip install numpy matplotlib
```

- [ ] **Step 2: Write `experiments/plot_common.py`**

A shared module with: a CSV loader that skips `#` provenance headers and returns
the header dict plus a numpy array; a single `apply_style()` that sets the figure
style once; and a `save(fig, path)` that writes the PNG and stamps the git SHA
into the figure footer. Every plot script uses these — no per-script styling.

- [ ] **Step 3: Write each experiment script**

Each is a short shell script that invokes `moldyn_run` with explicit flags and a
fixed seed, writes to `results/`, then calls its plot script. No hidden defaults.

- [ ] **Step 4: Write the `reproduce` Makefile target**

```make
.PHONY: build test reproduce clean
build:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
test: build
	ctest --test-dir build --output-on-failure
reproduce: build
	bash experiments/run_all.sh
clean:
	rm -rf build
```

- [ ] **Step 5: Run everything and inspect the results critically**

```bash
make reproduce
```
Then actually look at the figures. Check specifically: experiment 02's slope is
near 2; experiment 03's first peak position and height are close to the published
argon values; experiment 04's two diffusion estimates overlap within error bars.
**If a number disagrees with the literature, that is a finding to investigate,
not a tolerance to widen.** Report the disagreement rather than tuning until it
matches.

- [ ] **Step 6: Commit the experiment scripts and the results separately**

```bash
git add experiments Makefile
git commit -q -m "experiments: reproduction scripts for NIST, conservation, argon structure and transport

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"

git add results
git commit -q -m "results: committed figures and CSVs for experiments 01-05

Each file records the command line, git commit and host that produced it.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 7: Extend CI to run the fast experiments**

Add a CI step running experiments 01 and 02 (both take seconds) and failing the
build if the NIST deviations exceed tolerance or the fitted slope leaves
`[1.85, 2.15]`. Experiments 03-05 are too slow for every push; leave them to
`make reproduce`.

- [ ] **Step 8: Push and confirm CI is green**

```bash
git push && gh run watch
```

---

## Self-Review

**Spec coverage.** Spec section 4 (architecture) maps to the File Structure block
and Tasks 1-13. Section 5 (validation suite) maps as follows: NIST LJ energy and
virial → Task 5; force vs energy → Task 6; energy conservation order → Task 9;
time reversibility → Task 9; cell list correctness → Task 7; Verlet list
correctness → Task 8; minimum image → Task 2; equipartition and
Maxwell-Boltzmann → Task 10; Nose-Hoover conserved quantity → Task 11. The
NIST SPC/E Ewald row and the SETTLE row are Phase 2 and are deliberately out of
scope for this plan. Section 6 experiments 01-05 → Task 15; experiments 06-08 are
Phases 2 and 3. Section 7 (performance) is Phase 3, except that the cell and
Verlet lists land here because correctness requires them early. Section 8
(dependencies) → Task 0. Section 9 Phase 0 → Task 0, Phase 1 → Tasks 1-15.
Section 10 risk 2 (sampling time) is addressed by the block-averaging requirement
in Task 13 and the error-bar check in Task 15 step 5.

**Placeholder scan.** No TBD or TODO. Tasks 8, 10, 11 and 13 describe their tests
in prose rather than full code; each names the exact assertions and thresholds,
so nothing is left to the implementer's judgement. The one deliberate deferral is
`core/placeholder.cpp` in Task 0 step 7, which Task 2 step 2 explicitly deletes.

**Type consistency.** `EnergyVirial` is introduced in Task 5 and reused in Tasks
6-8. `computeForces` is overloaded consistently: `(System&)`, `(System&, const
CellList&)`, `(System&, const VerletList&)`, all returning `EnergyVirial`.
`Truncation::Truncated` and `Truncation::LinearForceShift` are used with those
exact spellings in Tasks 5, 6, 7 and 9. `System::temperature(dof)` takes the
degrees of freedom explicitly and is used that way in Task 10 with `3N-3`.
`Rng` is declared in Task 9 and used in Task 12's ideal-gas test. `fccLattice`
has one signature, used in Tasks 9, 10 and 11.

## Out of Scope for This Plan

Phase 2 (Ewald summation, SPC/E water, SETTLE, experiments 06-07), Phase 3
(NEON, threading, experiment 08) and Phase 4 (README rewrite) each get their own
plan once this one is green.
