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
#   >50 t*, production 100 t*. Six densities along T*~0.85.
#
# We reproduce NIST's own protocol as closely as our CLI allows: N=500
# (fcc, 5 cells/side -- 4*5^3=500 exactly), rc*=3.0, dt=0.005, equilibrate
# 50 t* (Langevin, to reach target T*), production 100 t* (NVE), plain
# truncation with the analytic tail correction on both potential energy
# and pressure (moldyn_run's default --truncation truncated; this IS
# NIST's "standard long range corrections" convention, and moldyn_run's
# thermo-out already includes the LRC pressure term as of commit c277bb2
# -- the script checks that flag explicitly below so a convention mismatch
# fails loudly instead of silently producing a ~100% pressure error).
#
# Tolerances (picked before looking at any result): |delta U*| <= 0.1
# (about 1.5-2% of the ~-6 magnitude here), |delta p*| <= 0.3 (an absolute
# bound, since p* itself ranges from 0.03 to 2.5 here and a relative bound
# is meaningless near p*=0.03). No plotting here -- see plot_exp05.py.
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
TARGET_T=0.85

DELTA_U_TOL=0.1
DELTA_P_TOL=0.3

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
# it is the source of the U*/p* averages but is not itself committed (six
# of them would add ~16MB of raw per-step logs to results/ for no benefit
# once the mean is taken -- the summary CSV's header has every fixed
# parameter and each row's own seed, so the exact command that reproduces
# any one of them is fully specified without storing its output). Written
# to a scratch directory that is cleaned up on exit.
SCRATCH="$(mktemp -d)"
trap 'rm -f "$ROWS"; rm -rf "$SCRATCH"' EXIT

while read -r rho nist_t nist_u nist_p seed; do
  [ -z "$rho" ] && continue
  thermo_csv="$SCRATCH/exp05_thermo_rho${rho}.csv"
  cmd="$BIN --cells $CELLS --density $rho --temperature $TARGET_T --cutoff $CUTOFF \
--skin $SKIN --dt $DT --equilibrate $EQUILIBRATE --production $PRODUCTION \
--thermostat none --seed $seed --thermo-out $thermo_csv"
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

  echo "exp05[rho*=$rho]: our U*=$mean_u p*=$mean_p T*=$mean_t (NIST U*=$nist_u p*=$nist_p T*=$nist_t)"
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
  echo "# equilibrate_steps: $EQUILIBRATE (Langevin, target T*=$TARGET_T)"
  echo "# production_steps: $PRODUCTION (NVE)"
  echo "# tolerance: |delta U*| <= $DELTA_U_TOL, |delta p*| <= $DELTA_P_TOL (absolute, picked a priori)"
  echo "# note: raw per-step thermo-out (30000 rows/density) is not committed -- every flag above"
  echo "#       is fixed except --density and --seed (this row's rho_star and seed columns), so"
  echo "#       'moldyn_run --cells $CELLS --density RHO --temperature $TARGET_T --cutoff $CUTOFF --skin $SKIN --dt $DT --equilibrate $EQUILIBRATE --production $PRODUCTION --thermostat none --seed SEED --thermo-out FILE' reproduces it exactly"
  echo "rho_star,seed,our_T_star,our_U_star,our_p_star,nist_T_star,nist_U_star,nist_p_star"
  cat "$ROWS"
} > "$SUMMARY_CSV"

echo "exp05: wrote $SUMMARY_CSV"

FAIL=0
while IFS=, read -r rho seed our_t our_u our_p nist_t nist_u nist_p; do
  ok=$(awk -v ou="$our_u" -v nu="$nist_u" -v op="$our_p" -v np="$nist_p" \
            -v ut="$DELTA_U_TOL" -v pt="$DELTA_P_TOL" '
    BEGIN {
      du = (ou-nu<0) ? nu-ou : ou-nu
      dp = (op-np<0) ? np-op : op-np
      print (du<=ut && dp<=pt) ? "yes" : "no"
    }')
  if [ "$ok" != "yes" ]; then
    echo "exp05[rho*=$rho]: FAIL -- our (U*=$our_u, p*=$our_p) vs NIST (U*=$nist_u, p*=$nist_p) outside tolerance" >&2
    FAIL=1
  else
    echo "exp05[rho*=$rho]: PASS -- our (U*=$our_u, p*=$our_p) vs NIST (U*=$nist_u, p*=$nist_p) within tolerance"
  fi
done < "$ROWS"

exit $FAIL
