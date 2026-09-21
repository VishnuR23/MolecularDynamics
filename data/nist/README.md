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
