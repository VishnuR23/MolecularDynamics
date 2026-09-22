# Experiment 04: the Einstein/Green-Kubo disagreement was centre-of-mass drift

**Date:** 2026-09-21
**Experiment:** `experiments/exp04_argon_diffusion.sh`
**Data:** `results/exp04_argon_diffusion.csv`, `results/exp04_msd_*_seed*.csv`
**Status:** diagnosed, corrected, and now passing. A second, independent
claim — agreement with Rahman's 1964 laboratory measurement — is recorded
at the end.

## Summary

Experiment 04 measures the self-diffusion coefficient of Lennard-Jones
argon two ways from the same trajectories: Einstein (slope of the mean
squared displacement) and Green-Kubo (integral of the velocity
autocorrelation function). They are two routes to the same number, so
they must agree.

They did not. `D_vacf` read systematically higher than `D_msd` at both
state points, by a similar absolute amount (+0.0026 and +0.0027), with
the same sign in 14 of 16 replicas. At the Rahman state point the gap
exceeded the two estimates' combined standard error and the experiment
failed. This was reported honestly and left unexplained.

The cause was uncorrected centre-of-mass motion. It is now fixed, and
the disagreement is gone: the offset has dropped by a factor of six and
changed sign, its residual is inside the error bar at both state points,
and the per-replica sign is now random. The predicted size of the effect
matched the observed size.

## The mechanism

`System::removeCenterOfMassMotion()` was called in exactly one place,
`util/lattice.cpp:62`, when the lattice is built. Every experiment then
ran thousands of steps of Langevin dynamics before production began.

Langevin does not conserve momentum. Each particle receives independent
noise, so the centre-of-mass momentum performs its own Ornstein-Uhlenbeck
random walk and thermalises like any other degree of freedom, reaching

    <v_com^2> = 3 k T / N

By the time NVE production started, the whole box was drifting at a
typical speed of sqrt(3T/N) — about 0.05 in reduced units at N = 864.
NVE then conserves that drift exactly, so it persists, unchanged, through
the entire production run.

A uniform drift contaminates both estimators, but not equally:

- **MSD** gains `v_com^2 * t^2`. Fitting a straight line to that over the
  window t* in [1.25, 3.755) adds `2 * <v_com^2> * t_mid` to the slope,
  hence `<v_com^2> * t_mid / 3` to D, with `t_mid` = 2.5.
- **VACF** gains a constant, non-decaying floor `<v_com^2>`, because a
  constant drift never decorrelates. Integrating that over the
  Green-Kubo window t* in [0, 5] adds `<v_com^2> * t_GK / 3` to D, with
  `t_GK` = maxLag*dt = 5.

The Green-Kubo window is twice the Einstein fit midpoint, so the VACF
estimator picks up twice the bias the MSD estimator does, and the
*difference* between them is biased by

    <v_com^2> * (t_GK - t_mid) / 3  =  3T/N * 2.5 / 3

which is +0.00228 at the Rahman point and +0.00211 at Verlet — against
the +0.00263 and +0.00269 actually observed. Same sign, same magnitude,
at both state points.

## The fix

Three changes, all in the physics rather than in a tolerance:

1. `apps/moldyn_run.cpp` re-zeroes the centre-of-mass momentum **once**,
   at the equilibration/production boundary. Velocity-Verlet and the
   Nosé-Hoover chain both conserve total momentum (the chain only ever
   rescales every velocity by one common factor), so once it is zero it
   stays zero for the whole production run.
2. It is **not** re-zeroed during Langevin production. Continually
   projecting out motion that the thermostat is supposed to be generating
   would perturb the dynamics being sampled.
3. The degrees of freedom used to report temperature now follow the
   dynamics: `3N` while Langevin is driving (the centre of mass is
   thermalised, so all 3N degrees of freedom carry kT/2), `3N-3` under
   the momentum-conserving integrators (those three are frozen at zero).
   `tests/test_thermostat.cpp` was corrected the same way. Using `3N-3`
   under Langevin had been inflating the reported temperature by
   `3N/(3N-3)` = 324/321 = 1.00935 at N = 108; the equipartition test's
   measured mean temperature moved from 1.00753 to **0.99743** as a
   result, which is the bias disappearing rather than a new one arriving.

## Result

n = 8 independent replicas per state point, 5000 steps Langevin
equilibration + 2000 steps NVE production at N = 864. "Before" is the
previously committed result; "after" is this tree.

| | before | after |
|---|---|---|
| **Rahman** (rho*=0.8177, T*=0.7880) | | |
| D_msd | 0.04535 ± 0.00062 | 0.04535 ± 0.00077 |
| D_vacf | 0.04798 ± 0.00110 | 0.04492 ± 0.00091 |
| D_vacf − D_msd | **+0.00263** | **−0.00043** |
| combined standard error | 0.00172 | 0.00168 |
| offset / combined SE | **1.53x — FAIL** | **0.25x — PASS** |
| **Verlet** (rho*=0.8442, T*=0.7280) | | |
| D_msd | 0.03567 ± 0.00160 | 0.03163 ± 0.00061 |
| D_vacf | 0.03836 ± 0.00241 | 0.03188 ± 0.00063 |
| D_vacf − D_msd | **+0.00269** | **+0.00024** |
| combined standard error | 0.00401 | 0.00124 |
| offset / combined SE | 0.67x — PASS | **0.20x — PASS** |
| **sign of D_vacf − D_msd** | positive in **14 of 16** replicas | positive in **6 of 16** |

The offset shrank by 6.1x at Rahman and 11x at Verlet, and the
per-replica sign went from strongly one-sided to indistinguishable from a
coin flip — which is what "the systematic part is gone and only scatter
remains" looks like.

Quantitatively, against the prediction above:

| quantity | predicted shift | observed shift |
|---|---|---|
| Rahman, D_vacf − D_msd | −0.00228 | −0.00306 |
| Verlet, D_vacf − D_msd | −0.00211 | −0.00245 |

Both within 0.0008 of prediction, with no parameter fitted.

## What the numbers do *not* show, stated plainly

The individual estimators moved less tidily than their difference.
`D_msd` was predicted to fall by 0.00228 (Rahman) and 0.00211 (Verlet);
it actually fell by 0.00000 and 0.00403. That asymmetry is worth naming
rather than glossing, because an 11% shift in Verlet's `D_msd` is the
sort of thing a reader comparing against the previous committed results
will notice.

Nothing else about the experiment changed: same N, cutoff, skin, dt,
5000/2000 equilibration/production steps, same eight seeds per state
point, same Einstein fit window (lags [250, 751)) and same Green-Kubo
window. The only physics change is the centre-of-mass re-zeroing. (One
non-physics change does re-roll the dice: `LangevinBAOAB` now uses
`std::fma` in its deterministic kick and drift, matching
`VelocityVerlet`. That alters trajectories in the last bits, which in a
chaotic system means completely different trajectories — but it cannot
shift a statistical mean.)

So the scatter is sampling noise, and it is the right size for it. Three
sources contribute to the uncertainty on an observed before/after shift:
the standard error of the before mean, of the after mean, and the spread
of the bias itself — `v_com^2` is a chi-squared variate with 3 degrees of
freedom, so it has a relative standard deviation of sqrt(2/3) = 82% per
replica, 29% after averaging 8. Combined:

| | predicted D_msd shift | observed | combined sd | discrepancy |
|---|---|---|---|---|
| Rahman | −0.00228 | −0.00000 | 0.00119 | 1.9 sigma |
| Verlet | −0.00211 | −0.00403 | 0.00182 | 1.1 sigma |

The two state points fall on *opposite* sides of the prediction, and
their mean shift, −0.00202, matches the predicted −0.00219 to 8%. So the
asymmetry is two draws from the expected distribution, not a second
effect.

This also explains why the *difference* is the better-determined
quantity, and why the prediction was framed on it: `D_msd` and `D_vacf`
carry the full run-to-run scatter of D itself, and the before and after
runs are independent samples of that scatter. Within a single trajectory,
though, both estimators see the same atoms and the same drift, so that
physical scatter largely cancels in their difference. The difference is
what the experiment's tolerance tests, and it is what matched.

One further consistency check falls out for free: removing a positive,
replica-varying bias should also *reduce* the replica-to-replica scatter.
Verlet's standard error on `D_msd` fell from 0.00160 to 0.00061 and on
`D_vacf` from 0.00241 to 0.00063. Rahman's changed little (0.00062 to
0.00077). Standard errors estimated from n = 8 are themselves uncertain
by about 27%, so this is supporting evidence, not proof.

## The second claim: agreement with a 1964 laboratory measurement

Rahman (1964), *Phys. Rev.* **136**, A405 reports a *measured*
self-diffusion coefficient for liquid argon at 94.4 K and 1.374 g/cm^3:

    D = 2.43e-5 cm^2/s

With the same Lennard-Jones constants this repository uses throughout
(sigma = 3.405 Å, eps/k_B = 119.8 K, M = 39.948 g/mol), the reduced unit
of diffusion is

    sigma * sqrt(eps/m) = 5.3767e-4 cm^2/s

so the experimental value in reduced units is

    D*_exp = 2.43e-5 / 5.3767e-4 = 0.0452

Our Einstein estimate at that same state point is **0.04535 ± 0.00077**
(n = 8). The agreement is **0.33%**, well inside our own error bar
(0.19 standard errors).

This is the strongest physics result in the repository: a diffusion
coefficient measured in a laboratory in 1964, reproduced from nothing but
a pair potential and Newton's equations, to a third of a percent. It was
previously computed and never stated. It is now a column (`D_exp_star`),
a status column (`experiment_status`), an assertion in
`exp04_argon_diffusion.sh`, and a line on
`results/exp04_argon_diffusion.png`.

The assertion's tolerance is **10% relative**, chosen before the
centre-of-mass fix was run and deliberately not tightened afterwards to
flatter the result. That is the level at which agreement of this kind is
meaningful: sigma and eps are themselves fitted parameters good to about
1%, D near the triple point is steeply sensitive to both, and Rahman's
number is a single 1964 measurement quoted to three figures. A tighter
bound would assert precision the comparison does not have.

## Status of the experiment

Both state points now pass both claims. Experiment 04 no longer
contributes to `make reproduce`'s non-zero exit; the remaining failure is
experiment 05's rho* = 0.900 row, written up in
[2026-09-21-rho090-outlier.md](2026-09-21-rho090-outlier.md).

## Reproducing every number above

```
make build
bash experiments/exp04_argon_diffusion.sh
```

Means, standard errors, the offset and both verdicts are written into the
header of `results/exp04_argon_diffusion.csv` (`result_rahman`,
`result_verlet`, `result_experiment`); the per-replica values are its
rows. The "before" column of the table comes from the previously
committed `results/exp04_argon_diffusion.csv` at commit `81f55ab`.
