# Phase 2: Ewald summation and SPC/E water

**Spec:** `docs/superpowers/specs/2026-09-20-md-engine-design.md` (section 9, Phase 2)
**Status:** conventions verified, ready to implement

## Why this order

The headline validation is a **static energy calculation** on configurations
already committed under `data/nist/spce/`. It costs milliseconds, like the
Lennard-Jones table. Water *simulation* — structure, diffusion — is the
expensive part and is scoped last and small, deliberately.

## Conventions, verified before planning

`tools/spce_reference.py` is an independent, deliberately naive Python
implementation. It reproduces all 28 NIST values. Everything below is what it
took to do that, so the C++ has an exact target rather than a description.

```
SPC/E          q_H = +0.42380 e, q_O = -2 q_H
               sigma = 3.16555789 A, eps/kB = 78.19743111 K, LJ on oxygen only
               r_OH = 1.0 A, angle H-O-H = 109.47 deg, rigid
Ewald          alpha = 5.6 / min(Lx,Ly,Lz)
               kmax = 5, including only k with |n|^2 < 27
               rcut = 10 A, tin-foil (conducting) boundary conditions
Dispersion     truncated at rcut, unshifted, with the analytic tail correction
Coulomb factor e^2/(4 pi eps0 kB) = 167100.956632 K*A per e^2
               (from NIST's own constants; energies are reported as E/kB in K)
```

Term by term, with the sign conventions that reproduce NIST:

```
Edisp     4 eps [(sig/r)^12 - (sig/r)^6] over OXYGEN PAIRS ONLY, r < rcut
ELRC      (8/3) pi M rho eps sig^3 [ (1/3)(sig/rcut)^9 - (sig/rcut)^3 ],  rho = M/V
Ereal     sum over INTERMOLECULAR pairs, r < rcut, of  C q_i q_j erfc(alpha r)/r
Efourier  C (2 pi / V) sum_{k != 0} exp(-k^2/4 alpha^2)/k^2 * |S(k)|^2,
          S(k) = sum_j q_j exp(i k . r_j)
Eself     -C (alpha/sqrt(pi)) sum_j q_j^2
Eintra    -C sum_molecules sum_{kappa<lambda} q_k q_l erf(alpha r_kl)/r_kl
```

### The convention NIST does not state

**Intramolecular distances need the minimum image.** Some molecules in the
sample files are stored split across the periodic boundary — oxygen on one
side, hydrogens wrapped to the other — so raw separations reach 19.9 A where
the bond is 1.0 A. In configuration 1, 28 of 300 intramolecular distances are
affected. Using raw separations leaves `Eintra` **4.8% low** while every other
term still matches to 1e-7, which is exactly the kind of error that looks like
success.

## Tolerances, derived not chosen

NIST prints 6 significant figures, so each published value carries a rounding
half-width. Components are compared against **their own** half-width, the same
way Phase 1 compares the Lennard-Jones table.

`Etotal` needs a different bound. It is a small residual of large cancelling
terms — at configuration 4, `-3.2e6` from terms of `+/-1.42e7` — so component
rounding amplifies. The total's bound is therefore **derived** by propagating
the components' half-widths through the cancellation:

| config | sum of component half-widths | \|Etotal\| | derived bound | reference achieves | margin |
|---|---|---|---|---|---|
| 1 | 10.6 | 4.886e5 | 2.16e-05 | 1.68e-06 | 12.9x |
| 2 | 15.5 | 1.066e6 | 1.46e-05 | 4.76e-06 | 3.1x |
| 3 | 15.5 | 1.715e6 | 9.04e-06 | 6.93e-06 | 1.3x |
| 4 | 105.6 | 3.205e6 | 3.29e-05 | 1.71e-05 | 1.9x |

Configuration 3 sits at 1.3x, so this is a real constraint rather than a
formality. Do not round these bounds up.

## Tasks

### A. Ewald and SPC/E — the headline, cheap

`src/moldyn/forcefield/ewald.{hpp,cpp}`, `spce.{hpp,cpp}`; charges and molecule
identity on `System`; `tests/test_nist_spce.cpp`.

Test each of the six components against its own rounding half-width on all four
configurations, and the total against its derived bound — **24 component
assertions plus 4 totals, not 4 total-only assertions.** A discrepancy must
localise to a term. That is the whole reason NIST publishes the breakdown.

Forces follow from the same kernel and are checked against the finite
difference of the energy, as in Phase 1.

Exit criterion: 28/28 green, and `tools/spce_reference.py` agrees with the C++
to the precision both share.

### B. SETTLE — rigid water

`src/moldyn/constrain/settle.{hpp,cpp}`, `tests/test_settle.cpp`.

Analytic constraint, not iterative. Tests: bond lengths and angle held to 1e-10
over a full trajectory; the constraint is a rotation so kinetic energy is not
silently injected; and energy conservation still shows its dt^2 slope with
constraints active.

### C. Water in motion — scoped small

One state point, SPC/E at 298 K, one short NVT run. Report the O-O radial
distribution function against the published first-peak position (~2.75 A).

**Diffusion is deliberately excluded.** Phase 1 established that converged
transport coefficients need long trajectories, and the rho*=0.9 investigation
showed what happens when they are not converged. A structure comparison is
cheap and honest; a diffusion number here would be neither.

## Out of scope

PME (Ewald is enough to validate the physics; PME is a performance technique
and belongs with Phase 3), flexible water, and any biomolecular force field.
