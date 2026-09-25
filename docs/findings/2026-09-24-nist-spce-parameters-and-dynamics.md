# NIST's SPC/E parameters reproduce static energies exactly, and do not conserve energy

**Date:** 2026-09-24
**Code:** `tests/test_rigid_water.cpp`, `src/moldyn/forcefield/{ewald,spce_forcefield}.cpp`
**Status:** understood, not a defect. It changes how one test is set up.

## Summary

The engine reproduces all 24 NIST SPC/E energy components and all 4 totals to
their published precision. Those same parameters, used for dynamics, put a
floor under energy conservation that no integrator or constraint can get below.

Measured with the full SPC/E interaction, the energy error stops scaling as
`dt^2` and the fitted log-log slope comes out near **1.2**. With a continuous
force field and the identical constraint code, it is **1.997**, with successive
error ratios of 4.00, 3.99 and 3.98.

The constraint is not at fault. The force field is discontinuous.

## Why

NIST's reference conditions are chosen so that a static energy on a fixed
configuration is unambiguous and reproducible. Two of them break continuity:

**Dispersion is truncated and unshifted.** At `rc = 10 Å` the Lennard-Jones
potential is still at `-0.314 K`, so an oxygen pair crossing the cutoff steps
the energy by that much.

**The Ewald screening is loose.** `alpha = 5.6 / min(L)` with `L = 20 Å` gives
`alpha = 0.28 Å^-1`, so `alpha * rc = 2.8` and

```
erfc(2.8) = 7.50e-05
```

7.5e-5 of the real-space Coulomb term is still present *at the cutoff*. A
charged pair crossing it steps the energy by about `0.450 K`, or `9.2e-07`
relative to a total of `4.9e5 K`.

For a static energy these are simply part of the defined convention, and both
sides of the comparison use them. For dynamics they are steps in the potential.

## Why it floors the slope rather than just adding noise

The number of cutoff crossings depends on elapsed **physical time**, not on the
timestep. A convergence test holds physical time fixed and varies `dt`, so every
run sees roughly the same number of crossings and inherits roughly the same
error from them. That is a `dt`-independent floor sitting under a `dt^2` signal,
and it flattens the fit from below — worst at the smallest timestep, where the
genuine integration error is smallest.

The data show exactly that shape. With the full interaction:

| dt (fs) | relative error | ratio |
|---|---|---|
| 0.5 | 7.18e-05 | |
| 1.0 | 9.95e-05 | 1.39 |
| 2.0 | 2.56e-04 | 2.58 |
| 4.0 | 8.25e-04 | 3.22 |

The ratio climbs toward 4 as the `dt^2` term grows past the floor. With a
continuous force field the floor is gone and the ratios are flat at 4:

| dt (fs) | relative error | ratio |
|---|---|---|
| 0.5 | 1.43e-05 | |
| 1.0 | 5.74e-05 | 4.00 |
| 2.0 | 2.29e-04 | 3.99 |
| 4.0 | 9.12e-04 | 3.98 |

## What was ruled out first

Two plausible explanations were tested and rejected before landing on this one.

**Constraint solver tolerance.** SHAKE converges to a relative constraint
violation of `1e-12` by default. Tightening it to `1e-13` and `1e-14` changed
the measured errors *not at all* — identical to every printed digit. (At
`1e-15` it stops converging, which is the double-precision limit for this
residual.) So the floor is not solver convergence.

**Round-off in the velocity correction.** The position reset feeds back into the
half-step velocity as `(r_constrained - r_unconstrained) / dt`, and dividing by
`dt` amplifies round-off as `dt` shrinks — the right shape for the observed
floor. It is far too small: position round-off of `eps * 20 Å ≈ 4.4e-15 Å`
gives a relative energy error near `6e-13`, eight orders below what is observed.

## Consequence

`tests/test_rigid_water.cpp`'s convergence test uses oxygen-oxygen dispersion
with the force shifted to zero at the cutoff, and no electrostatics. That is
deliberate and documented in the test itself: the test exists to measure the
constraint, so the force field must not contribute a floor that swamps it.

Energy-conserving SPC/E dynamics would want a shifted or switched dispersion
term and a tighter Ewald screening — a larger `alpha`, a longer real-space
cutoff, or both. Those are different parameters from the ones that reproduce
the NIST reference, so an engine that wants to do both needs to carry both
conventions and be explicit about which is in use. This one currently carries
NIST's.
