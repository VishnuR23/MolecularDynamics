#!/usr/bin/env bash
# Experiment 05b: diagnostic evidence for the rho*=0.9 outlier found in
# experiment 05. exp05_lj_eos.sh's report claims that state point is not a
# frozen/non-ergodic crystal -- this script is what that claim rests on,
# committed and reproducible rather than asserted from an ad hoc run.
#
# Same physical state point as exp05's rho*=0.9 row, same N=500/cutoff/dt,
# equilibrated the same way (Langevin, target = NIST's own quoted T* for
# this row, 0.851). Production here uses NVE, not the NVT (Nose-Hoover)
# exp05 itself now uses: measuring a transport coefficient (D, from MSD
# and VACF) requires real, unthermostatted Newtonian dynamics -- a
# thermostatted production would bias the velocity autocorrelation function
# and hence D. This is the same choice experiment 04 makes (Langevin
# equilibration, NVE production) for the same reason.
#
# Writes g(r) and MSD/VACF to results/ (this is the fix: the original
# investigation ran these as one-off, uncommitted diagnostic commands).
# Asserts D_msd stays above a floor that distinguishes "diffusing, just
# slower than a fully-equilibrated liquid" from "frozen, effectively
# D=0" -- 0.005, chosen a priori as well below any liquid D measured
# anywhere else in this project (experiment 04's smallest is ~0.029,
# NIST's own tabulated D* at this exact state point is 0.027 -- see
# data/nist/README.md) and well above true zero. No plotting here -- see
# plot_exp05b.py.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT/build/apps/moldyn_run"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

GIT_SHA="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
HOST="$(uname -s) $(uname -m)"

CELLS=5              # N = 4*5^3 = 500, matching exp05 and NIST
RHO=0.900
TARGET_T=0.851        # NIST's own quoted T* for this exact row
CUTOFF=3.0
SKIN=0.3
DT=0.005
EQUILIBRATE=10000     # 50 t*, matching exp05
PRODUCTION=20000      # 100 t*, matching exp05
SEED=5006             # same seed as exp05's rho*=0.9 row

D_FLOOR=0.005          # "not frozen" floor -- see header
NIST_D_STAR=0.027      # data/nist/README.md, same state point

RDF_CSV="$RESULTS/exp05b_rdf_rho0.900.csv"
MSD_CSV="$RESULTS/exp05b_msd_rho0.900.csv"
SUMMARY_CSV="$RESULTS/exp05b_rho090_diagnostic.csv"

cmd="$BIN --cells $CELLS --density $RHO --temperature $TARGET_T --cutoff $CUTOFF \
--skin $SKIN --dt $DT --equilibrate $EQUILIBRATE --production $PRODUCTION \
--thermostat none --seed $SEED --rdf-out $RDF_CSV --msd-out $MSD_CSV"
echo "exp05b: $cmd"
eval "$cmd"

D_MSD=$(grep '^# D_msd' "$MSD_CSV" | sed -E 's/.*: *//')
D_VACF=$(grep '^# D_vacf' "$MSD_CSV" | sed -E 's/.*: *//')

# First RDF peak, same method as exp03: max g(r) with r* in a window
# containing the first coordination shell.
read -r PEAK_R PEAK_G <<EOF
$(awk -F, -v lo=0.85 -v hi=1.40 '
    /^#/ { next }
    !seen_header { seen_header=1; next }
    { r=$1; g=$2; if (r>=lo && r<=hi && g>maxg) { maxg=g; maxr=r } }
    END { printf "%.6f %.6f", maxr, maxg }
  ' "$RDF_CSV")
EOF

# g(r) at the largest tabulated r: how close to 1 (the ideal-gas/bulk
# value a liquid's g(r) relaxes to) it has settled by the edge of the
# accessible range -- a crystal would still show strong, sharp deviation
# from 1 out at r* close to L/2, a liquid decays to ~1 well before it.
TAIL_G=$(awk -F, '
    /^#/ { next }
    !seen_header { seen_header=1; next }
    { last=$2 }
    END { printf "%.4f", last }
  ' "$RDF_CSV")

DELTA_D=$(awk -v d="$D_MSD" -v n="$NIST_D_STAR" 'BEGIN { printf "%.6f", d-n }')

{
  echo "# purpose: evidence for experiment 05's rho*=0.9 outlier investigation -- see task-15 report"
  echo "# command: $cmd"
  echo "# git_sha: $GIT_SHA"
  echo "# host: $HOST"
  echo "# units: reduced Lennard-Jones units (sigma = epsilon = m = 1)"
  echo "# state_point: rho*=$RHO, T*=$TARGET_T (NIST's own quoted T* for this row)"
  echo "# nist_D_star: $NIST_D_STAR (data/nist/README.md, same rho*/T*, fully-equilibrated liquid)"
  echo "# d_floor: $D_FLOOR (a priori 'not frozen' threshold, see script header)"
  echo "# tolerance: D_msd > $D_FLOOR"
  echo "rho_star,T_star,seed,D_msd,D_vacf,nist_D_star,delta_D_msd,rdf_first_peak_r,rdf_first_peak_g,rdf_tail_g"
  echo "$RHO,$TARGET_T,$SEED,$D_MSD,$D_VACF,$NIST_D_STAR,$DELTA_D,$PEAK_R,$PEAK_G,$TAIL_G"
} > "$SUMMARY_CSV"

echo "exp05b: D_msd=$D_MSD  D_vacf=$D_VACF  (NIST D*=$NIST_D_STAR, delta=$DELTA_D)"
echo "exp05b: RDF first peak r*=$PEAK_R g=$PEAK_G, tail g(r_max)=$TAIL_G"
echo "exp05b: wrote $SUMMARY_CSV, $MSD_CSV, $RDF_CSV"

PASS=$(awk -v d="$D_MSD" -v floor="$D_FLOOR" 'BEGIN { print (d>floor) ? "yes" : "no" }')
if [ "$PASS" != "yes" ]; then
  echo "exp05b: FAIL -- D_msd=$D_MSD does not clear the not-frozen floor $D_FLOOR" >&2
  exit 1
fi

echo "exp05b: PASS -- D_msd=$D_MSD clears the not-frozen floor $D_FLOOR (suppressed vs. NIST's $NIST_D_STAR, but not zero)"
