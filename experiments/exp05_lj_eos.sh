#!/usr/bin/env bash
# Experiment 05: Lennard-Jones equation of state vs. NIST's own tabulated
# values -- not a digitised plot from a paper. Source and exact quoted
# numbers are in data/nist/README.md ("Lennard-Jones fluid equation of
# state"), retrieved 2026-09-21 from
# https://mmlapps.nist.gov/srs/LJ_PURE/md.htm. That page tabulates one MD
# isotherm with BOTH U* and p* (its EOS-TMMC page has a much finer density
# grid but pressure only, so cannot be used for a direct U* comparison):
#
#   NVE molecular dynamics, velocity-Verlet, dt=0.005 t*, N=500,
#   truncation "3 sigma + standard long range corrections", equilibration
#   >50 t*, production 100 t*. Six densities, each at ITS OWN measured T*
#   (0.849-0.853, not a single fixed 0.85 -- NIST's own NVE production
#   simply drifts to whatever T results from their initial velocities).
#
# CORRECTED PROTOCOL (post-review): the first version of this script ran
# NVE production at a single target T*=0.85 for every density. NVE does
# not pin T -- it only conserves E -- so our own achieved T* drifted
# per-density, from 0.8256 to 0.8713, while NIST's six rows sit in the
# tight range 0.849-0.853. A reviewer spot-check found the sign of every
# passing row's delta-U* tracked the sign of its own T* offset, with a
# magnitude consistent with dU* ~ Cv* * dT* -- i.e. part of what looked
# like engine deviation was actually a different-state-point comparison,
# not error. Fixed at the source: production now runs NVT (Nose-Hoover
# chain, already validated elsewhere in this project -- its extended
# Hamiltonian conserves to 7e-5, tests/test_thermostat.cpp), pinned at
# EACH ROW'S OWN NIST-quoted T*, not a shared nominal 0.85. The achieved
# T* is still measured and reported per row (our_T_star column) and
# asserted close to target, so this is verified, not assumed. At N=500 the
# NVT vs NVE ensemble difference in U and p is far smaller than the ~0.03
# temperature offset the old NVE-at-0.85 protocol carried.
#
# N=500 (fcc, 5 cells/side -- 4*5^3=500 exactly), rc*=3.0, dt=0.005,
# equilibrate 50 t* (Langevin, to reach each row's own target T*),
# production 100 t* (NVT, Nose-Hoover, pinned at that same T*), plain
# truncation with the analytic tail correction on both potential energy
# and pressure (moldyn_run's default --truncation truncated; this IS
# NIST's "standard long range corrections" convention, and moldyn_run's
# thermo-out already includes the LRC pressure term as of commit c277bb2
# -- the script checks that flag explicitly below so a convention mismatch
# fails loudly instead of silently producing a ~100% pressure error).
#
# Tolerances (picked before looking at any NVT result -- the U*/p*
# tolerances are unchanged from the original, pre-registered choice;
# T_TOL is new, for the temperature-pinning check this fix adds):
#   |delta U*| <= 0.1  (about 1.5-2% of the ~-6 magnitude here)
#   |delta p*| <= 0.3  (absolute -- p* itself ranges 0.03 to 2.5 here, a
#                        relative bound is meaningless near p*=0.03)
#   |delta T*| <= 0.02 (about 2.4% of T*~0.85 -- confirms the Nose-Hoover
#                        chain actually landed on NIST's own T*, so the
#                        U*/p* comparison above is now apples to apples)
# No plotting here -- see plot_exp05.py.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT/build/apps/moldyn_run"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

GIT_SHA="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
HOST="$(uname -s) $(uname -m)"

CELLS=5   # N = 4*5^3 = 500, matching NIST's N exactly
N_ATOMS=500
CUTOFF=3.0
SKIN=0.3
DT=0.005
EQUILIBRATE=10000   # 50 t*
PRODUCTION=20000    # 100 t*

DELTA_U_TOL=0.1
DELTA_P_TOL=0.3
DELTA_T_TOL=0.02

# rho*, NIST T*, NIST U*, NIST p*, seed
STATE_POINTS="
0.776 0.851 -5.517 0.030 5001
0.780 0.853 -5.533 0.072 5002
0.820 0.852 -5.803 0.573 5003
0.840 0.851 -5.909 0.910 5004
0.860 0.849 -6.027 1.282 5005
0.900 0.851 -6.234 2.544 5006
"

SUMMARY_CSV="$RESULTS/exp05_lj_eos.csv"
ROWS="$(mktemp)"
# Per-density thermo-out is a full per-step trajectory (30000 rows each);
# it is the source of the U*/p*/T* averages but is not itself committed
# (six of them would add ~16MB of raw per-step logs to results/ for no
# benefit once the mean is taken -- the summary CSV's header has every
# fixed parameter and each row's own seed and target T*, so the exact
# command that reproduces any one of them is fully specified without
# storing its output). Written to a scratch directory cleaned up on exit.
SCRATCH="$(mktemp -d)"
trap 'rm -f "$ROWS"; rm -rf "$SCRATCH"' EXIT

while read -r rho nist_t nist_u nist_p seed; do
  [ -z "$rho" ] && continue
  thermo_csv="$SCRATCH/exp05_thermo_rho${rho}.csv"
  cmd="$BIN --cells $CELLS --density $rho --temperature $nist_t --cutoff $CUTOFF \
--skin $SKIN --dt $DT --equilibrate $EQUILIBRATE --production $PRODUCTION \
--thermostat nose-hoover --seed $seed --thermo-out $thermo_csv"
  echo "exp05[rho*=$rho]: $cmd"
  eval "$cmd"

  if ! grep -q '^# pressure_includes_tail_correction: yes' "$thermo_csv"; then
    echo "exp05[rho*=$rho]: FATAL -- thermo-out is not tail-corrected; would mismatch NIST's sLRC convention" >&2
    exit 1
  fi

  read -r mean_u mean_p mean_t <<EOF
$(awk -F, '
    /^#/ { next }
    !seen_header { seen_header=1; next }
    $3=="production" { n++; u+=$4; p+=$8; t+=$7 }
    END { printf "%.6f %.6f %.6f", u/n/N, p/n, t/n }
  ' N="$N_ATOMS" "$thermo_csv")
EOF

  echo "exp05[rho*=$rho]: our U*=$mean_u p*=$mean_p T*=$mean_t (NIST U*=$nist_u p*=$nist_p T*=$nist_t, target=$nist_t)"
  echo "$rho,$seed,$mean_t,$mean_u,$mean_p,$nist_t,$nist_u,$nist_p" >> "$ROWS"
done <<< "$STATE_POINTS"

{
  echo "# source: https://mmlapps.nist.gov/srs/LJ_PURE/md.htm (NIST SRSW LJ fluid MD results)"
  echo "# retrieved: 2026-09-21"
  echo "# reference_convention: 3 sigma cutoff + standard long range corrections (matches our --truncation truncated)"
  echo "# git_sha: $GIT_SHA"
  echo "# host: $HOST"
  echo "# units: reduced Lennard-Jones units (sigma = epsilon = m = 1)"
  echo "# cells: $CELLS (N=$N_ATOMS, matching NIST's N exactly)"
  echo "# cutoff: $CUTOFF"
  echo "# skin: $SKIN"
  echo "# dt: $DT"
  echo "# equilibrate_steps: $EQUILIBRATE (Langevin, target = this row's nist_T_star)"
  echo "# production_steps: $PRODUCTION (NVT, Nose-Hoover chain, pinned at this row's nist_T_star)"
  echo "# protocol_correction: production changed from NVE at a shared T*=0.85 to NVT (Nose-Hoover)"
  echo "#       pinned at each row's own NIST-quoted T*, after a review found the old NVE T* drift"
  echo "#       (0.8256-0.8713) confounded the U*/p* comparison -- see task-15 report"
  echo "# tolerance: |delta U*| <= $DELTA_U_TOL, |delta p*| <= $DELTA_P_TOL, |delta T*| <= $DELTA_T_TOL (all absolute, picked a priori)"
  echo "# note: raw per-step thermo-out (30000 rows/density) is not committed -- every flag above"
  echo "#       is fixed except --density, --temperature and --seed (this row's rho_star, nist_T_star"
  echo "#       and seed columns), so 'moldyn_run --cells $CELLS --density RHO --temperature NIST_T"
  echo "#       --cutoff $CUTOFF --skin $SKIN --dt $DT --equilibrate $EQUILIBRATE --production $PRODUCTION"
  echo "#       --thermostat nose-hoover --seed SEED --thermo-out FILE' reproduces it exactly"
  echo "rho_star,seed,our_T_star,our_U_star,our_p_star,nist_T_star,nist_U_star,nist_p_star"
  cat "$ROWS"
} > "$SUMMARY_CSV"

echo "exp05: wrote $SUMMARY_CSV"

FAIL=0
while IFS=, read -r rho seed our_t our_u our_p nist_t nist_u nist_p; do
  result=$(awk -v ou="$our_u" -v nu="$nist_u" -v op="$our_p" -v np="$nist_p" -v ot="$our_t" -v nt="$nist_t" \
            -v ut="$DELTA_U_TOL" -v pt="$DELTA_P_TOL" -v tt="$DELTA_T_TOL" '
    BEGIN {
      du = (ou-nu<0) ? nu-ou : ou-nu
      dp = (op-np<0) ? np-op : op-np
      dtem = (ot-nt<0) ? nt-ot : ot-nt
      up_ok = (du<=ut && dp<=pt) ? "yes" : "no"
      t_ok = (dtem<=tt) ? "yes" : "no"
      printf "%s %s %.6f %.6f %.6f", up_ok, t_ok, du, dp, dtem
    }')
  read -r up_ok t_ok du dp dtem <<< "$result"

  if [ "$t_ok" != "yes" ]; then
    echo "exp05[rho*=$rho]: FAIL -- Nose-Hoover did not pin T* (our=$our_t, target=$nist_t, |delta|=$dtem > $DELTA_T_TOL)" >&2
    FAIL=1
  elif [ "$up_ok" != "yes" ]; then
    echo "exp05[rho*=$rho]: FAIL -- our (U*=$our_u, p*=$our_p) vs NIST (U*=$nist_u, p*=$nist_p) outside tolerance (T* pinned OK, |dT*|=$dtem)" >&2
    FAIL=1
  else
    echo "exp05[rho*=$rho]: PASS -- our (U*=$our_u, p*=$our_p, T*=$our_t) vs NIST (U*=$nist_u, p*=$nist_p, T*=$nist_t) within tolerance"
  fi
done < "$ROWS"

exit $FAIL
