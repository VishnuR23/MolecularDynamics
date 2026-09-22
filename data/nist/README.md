# NIST reference data provenance

Two datasets from the NIST Standard Reference Simulation Website, used to
validate the Lennard-Jones and SPC/E water force field implementations
against published values.

## Lennard-Jones fluid (`lj/`)

- Source: https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-reference-calculations-cuboid-cell
- Downloaded: 2026-09-20

| File | sha256 |
|---|---|
| `lj_sample_config_periodic1.txt` | `6fa4bf9d4cb84c07d4766e4c1143c5fdf0c560c47a695336c6d06b19e425d295` |
| `lj_sample_config_periodic2.txt` | `d6714977c486280577f2fe69fdfd52188b47a1bdace653edf47ad44b29f37c9a` |
| `lj_sample_config_periodic3.txt` | `c8e13fef51b21cea5f958ae911c0f45fe3f7ef4f92578d68569258bbb29113ec` |
| `lj_sample_config_periodic4.txt` | `aa823314ac4abf597cbe9faf9b9caab267a78ddab1433779a6d7da72f59e1883` |

### Reference values (reduced units, sigma = epsilon = m = 1)

Four configurations, three truncation schemes: LRC = truncated with analytic
long-range correction, LFS = linear force shift at the cutoff.

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

Our LJ conventions (unshifted truncation, undivided virial sum) were verified
against all 12 published values above.

## Lennard-Jones fluid equation of state

Unlike `lj/` and `spce/` above, this is not a downloaded reference-config
file (there is none to download -- NIST's page presents this table
directly in HTML). Nothing is stored under `data/nist/` for it; the
quoted table below is copied here verbatim as the provenance record.

Used by experiment 05 (`experiments/exp05_lj_eos.sh`) to compare our
simulated pressure and potential energy against NIST's own published values,
rather than a digitised plot from a paper, and by experiment 05b
(`experiments/exp05b_rho090_diagnostic.sh`) for the D* comparison described
below.

- Source: https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-properties,
  MD results table at https://mmlapps.nist.gov/srs/LJ_PURE/md.htm
- Retrieved: 2026-09-21

This page tabulates only one MD isotherm with both U* and p* reported
together (its EOS-TMMC page has a much finer density grid but pressure
only, no energy, so it cannot be used for a direct U* comparison). We use
that one isotherm exactly as published, at NIST's own state points, rather
than picking our own density grid and interpolating.

Stated methodology (quoted from the page): method NVE molecular dynamics,
integrator velocity-Verlet, time step 0.005 t*, N = 500, truncation "3σ +
standard long range corrections", equilibration >50 t*, production 100 t*.
"Standard long range corrections are included in the tabulated values of
energy and pressure. ... The reported pressure was calculated using the
virial expression." This is exactly our `--truncation truncated` convention
with the analytic LJ tail correction applied to both potential energy and
pressure (moldyn_run commit c277bb2), so the two are directly comparable
with no unit or convention conversion.

The standard uncertainty in the last digit is quoted in parentheses,
obtained from block averages of length 2-4 time units (67% confidence).

### Reference values (reduced units, sigma = epsilon = m = 1, rc* = 3.0)

| T* | rho* | U* | p* | D* |
|---|---|---|---|---|
| 0.851(1) | 0.776 | -5.517(1) | 0.030(6) | 0.061 |
| 0.853(1) | 0.780 | -5.533(1) | 0.072(5) | 0.063 |
| 0.852(1) | 0.820 | -5.803(1) | 0.573(5) | 0.048 |
| 0.851(1) | 0.840 | -5.909(1) | 0.910(6) | 0.042 |
| 0.849(1) | 0.860 | -6.027(1) | 1.282(5) | 0.035 |
| 0.851(1) | 0.900 | -6.234(1) | 2.544(6) | 0.027 |

We run our own simulation at each of these six (T*, rho*) pairs with
N = 500 (fcc, 5 cells/side), rc* = 3.0, dt = 0.005, matching NIST's own
protocol as closely as our CLI allows: Langevin equilibration to each row's
own target T* (not a single nominal 0.85 for all six -- NIST's tabulated
T* values are themselves each a distinct row, 0.849-0.853), then **NVT**
production (Nose-Hoover chain, pinned at that same target T*). Production
was originally NVE, matching NIST's stated method; a review found NVE let
our achieved T* drift per-density (0.8256-0.8713) against NIST's tight
0.849-0.853 band, confounding part of the U*/p* comparison with a
different-state-point effect rather than engine error. NVT removes that
confound at the cost of no longer matching NIST's own ensemble exactly;
at N=500 that ensemble difference is far smaller than the temperature
offset it replaces. `our_T_star` in the results CSV is the achieved mean
production temperature, so this is a checked, reported quantity, not an
unexamined one. NIST's own equilibration method is unspecified beyond
">50 t*". See `experiments/exp05_lj_eos.sh` for the exact commands and
`results/exp05_lj_eos.csv` for our measured values and their deviation
from this table: `our_T_star`, `our_U_star`, `our_p_star`, the matching
`nist_*` columns copied from the rows above, the signed deviations
`delta_T_star`, `delta_U_star`, `delta_p_star` (all `ours - NIST`), and a
per-row `status`. Five rows pass; `rho*=0.900` does not, and the next
section is why.

`experiments/exp05b_rho090_diagnostic.sh` characterises what the
simulation actually reaches at rho*=0.900, across three independent
seeds. The answer is that it never melts: it stays in the fcc lattice it
was started from. All three seeds give D_msd between 1.5e-5 and 3.7e-5
against NIST's tabulated **D*=0.027** for the same (rho*, T*) -- about
1000x smaller, and indistinguishable from zero. All three fail the
experiment's a-priori "not frozen" floor of D_msd > 0.005.

Four signatures, all read from the committed trajectories:

- **The MSD is flat.** It rises ballistically to 0.085 by t*=1 and then
  stays inside 0.0790-0.0853 all the way to t*=10. A liquid's MSD grows
  linearly forever; a bounded one is atoms vibrating about fixed sites.
- **The vibration amplitude is a melting solid's.** MSD = 0.0805 at
  t*=10 is 0.164 sigma of rms displacement per axis against a
  nearest-neighbour distance of 1.098 sigma: a Lindemann ratio of
  **0.149**, where the textbook criterion for a crystal at its melting
  point is about 0.15.
- **The VACF backscatters and dies** (minimum -0.303 at t*=0.20, down to
  0.0009 by t*=10). Backscattering alone does not separate a liquid from
  a solid -- an atom in a cage and an atom on a lattice site both
  reverse -- but with a flat MSD it means the reversals never add up to
  transport.
- **g(r) keeps order the reference liquid has lost.** The principal
  peaks decay monotonically, 3.125 -> 1.531 -> 1.318; an earlier version
  of this note called that "an oscillation regrowing rather than damping
  out", which the committed data never supported. The accurate statement
  is stronger. Against experiment 03's Rahman g(r) (identical code and
  normalisation, lower density rho*=0.8177, a genuine liquid), |g-1| at
  the third principal peak is **0.318 against 0.087 -- an envelope 3.7x
  further from the bulk value**, and where the reference has settled to
  0.98-1.02 by r*≈4.5-5, ours is still 1.102 at r*=4.02. There is also a
  **split second-shell shoulder** our g(r) has and the reference does
  not: two minima at r*=2.367 (g=0.846) and r*=2.789 (g=0.812) with a
  sub-unity maximum at r*=2.588 (g=0.960) between them, where the
  reference has one smooth trough at r*=2.533 (g=0.853). That splitting
  is the fcc lattice's distinct neighbour distances resolving where a
  liquid's would be smeared into one.

Both signs of the EOS disagreement follow: an ordered lattice sits deeper
in the pair potential (U* too negative by 0.390) and relieves the
short-range repulsion that dominates the virial at this density (p* too
low by 2.452).

This supersedes an earlier n=1 reading of the same state point, which
measured D_msd=0.0203 and concluded "partially ordered but still mobile".
About 0.0046 of that was not diffusion at all but uncorrected
centre-of-mass drift, which adds `<v_com^2> * t_mid / 3` to any Einstein
estimate; with the drift fixed the same seed gives 1.5e-5. The
diagnostic now runs three seeds rather than one, and the threshold was
left untouched.

One limitation, stated rather than glossed: the Green-Kubo window is
capped at t*=10 by `kMaxLagCap=2000` (`apps/moldyn_run.cpp`), 10% of the
100 t* trajectory. It does not affect the conclusion here -- both
estimators agree on "indistinguishable from zero" -- but measuring a real
D at this density would need more. Raising `maxLag` is not a config
change: `Msd`/`Vacf` accumulation is O(frames x maxLag x natoms), and a
t*=50 window at the current sampling rate is roughly 1e11 operations.
Extending the correlation window properly needs a frame-sampling stride,
which is new driver scope, recorded as future work for Phase 3's
performance pass.

See `results/exp05b_rho090_diagnostic.csv`/`.png` and
[`docs/findings/2026-09-21-rho090-outlier.md`](../../docs/findings/2026-09-21-rho090-outlier.md)
for the full discussion, and
[`docs/findings/2026-09-21-exp04-einstein-vs-green-kubo.md`](../../docs/findings/2026-09-21-exp04-einstein-vs-green-kubo.md)
for the separate experiment 04 finding.

## SPC/E water (`spce/`)

- Source: https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10a-cutoff
- Downloaded: 2026-09-20

| File | sha256 |
|---|---|
| `spce_sample_config_periodic1.txt` | `a2b01a73b918444cd9a6a47895999924bb87b45b9315cf836710d7d5e2d840d6` |
| `spce_sample_config_periodic2.txt` | `a52db516d5c8450f08cc71d220bfda26b4936596ba17210dd04297421da8238b` |
| `spce_sample_config_periodic3.txt` | `5d89c3e07199c2fffff11076b35f5b40d3dce20e2dbad128719216634b0bcf92` |
| `spce_sample_config_periodic4.txt` | `095a4ee6a64f281ebc0b77374e410e7dc2338e6bd0963643f74001b2946e68ee` |

### Reference values (energies as E/kB in Kelvin)

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
