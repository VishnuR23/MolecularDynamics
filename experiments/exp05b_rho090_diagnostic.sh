#!/usr/bin/env bash
# Experiment 05b: what state is the system actually in at rho*=0.900,
# T*=0.851 -- the one density where experiment 05 disagrees with NIST?
#
# NIST tabulates that row as a liquid, with a self-diffusion coefficient
# D* = 0.027. The question this experiment asks is whether our simulation
# reaches that liquid at all. The a-priori test, unchanged since it was
# written: D_msd must clear 0.005, a floor chosen to sit well below any
# liquid D measured anywhere in this project (experiment 04's smallest is
# ~0.031, NIST's own D* here is 0.027) and well above true zero. Clearing
# it means "diffusing, even if more slowly than NIST's liquid"; failing
# it means "frozen".
#
# THREE INDEPENDENT SEEDS, not one. This experiment previously ran a
# single trajectory, and that turned out not to be enough: the same seed
# gave D_msd = 0.020 (mobile) under one revision of the engine and
# D_msd = 1.5e-5 (frozen) under the next, which differed only by a
# centre-of-mass correction and a rounding detail. A state point whose
# answer flips between two last-bit-different trajectories is a state
# point one trajectory cannot characterise. Replicas are the cheapest
# honest fix; the threshold itself is untouched.
#
# Same physical state point as exp05's rho*=0.9 row, same N=500/cutoff/dt,
# equilibrated the same way (Langevin, target = NIST's own quoted T* for
# this row, 0.851). Production uses NVE, not the NVT (Nose-Hoover) exp05
# itself uses: measuring a transport coefficient (D, from MSD and VACF)
# requires real, unthermostatted Newtonian dynamics -- a thermostatted
# production would bias the velocity autocorrelation function and hence D.
# This is the same choice experiment 04 makes, for the same reason.
#
# Writes g(r) and MSD/VACF per seed to results/, so every claim made about
# this state point is backed by a committed, re-readable trajectory rather
# than an ad hoc run. Findings:
# docs/findings/2026-09-21-rho090-outlier.md. No plotting here -- see
# plot_exp05b.py.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT/build/apps/moldyn_run"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

# `git describe --always --dirty`, not `rev-parse --short HEAD`: a run
# from a modified working tree is stamped `<sha>-dirty`, so it cannot
# masquerade as the clean commit it was derived from.
GIT_SHA="$(git -C "$ROOT" describe --always --dirty 2>/dev/null || echo unknown)"
HOST="$(uname -s) $(uname -m)"

CELLS=5              # N = 4*5^3 = 500, matching exp05 and NIST
RHO=0.900
TARGET_T=0.851        # NIST's own quoted T* for this exact row
CUTOFF=3.0
SKIN=0.3
DT=0.005
EQUILIBRATE=10000     # 50 t*, matching exp05
PRODUCTION=20000      # 100 t*, matching exp05

# 5006 is exp05's own seed for this row, so one replica is exactly the
# trajectory the failing EOS point came from; the other two are
# independent.
SEEDS=(5006 5007 5008)
PRIMARY_SEED=5006

D_FLOOR=0.005          # "not frozen" floor -- see header
NIST_D_STAR=0.027      # data/nist/README.md, same state point

SUMMARY_CSV="$RESULTS/exp05b_rho090_diagnostic.csv"
ROWS="$(mktemp)"
trap 'rm -f "$ROWS"' EXIT

FAIL=0

for SEED in "${SEEDS[@]}"; do
  RDF_CSV="$RESULTS/exp05b_rdf_rho0.900_seed${SEED}.csv"
  MSD_CSV="$RESULTS/exp05b_msd_rho0.900_seed${SEED}.csv"

  cmd="$BIN --cells $CELLS --density $RHO --temperature $TARGET_T --cutoff $CUTOFF \
--skin $SKIN --dt $DT --equilibrate $EQUILIBRATE --production $PRODUCTION \
--thermostat none --seed $SEED --rdf-out $RDF_CSV --msd-out $MSD_CSV"
  echo "exp05b[seed=$SEED]: $cmd"
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

  # g(r) at the largest tabulated r: how close to 1 (the bulk value a
  # liquid's g(r) relaxes to) the structure has settled by the edge of
  # the accessible range. A solid still shows strong deviation from 1 out
  # near L/2; a liquid is back to ~1 well before that.
  TAIL_G=$(awk -F, '
      /^#/ { next }
      !seen_header { seen_header=1; next }
      { last=$2 }
      END { printf "%.4f", last }
    ' "$RDF_CSV")

  DELTA_D=$(awk -v d="$D_MSD" -v n="$NIST_D_STAR" 'BEGIN { printf "%.6f", d-n }')
  STATUS=$(awk -v d="$D_MSD" -v floor="$D_FLOOR" 'BEGIN { print (d>floor) ? "PASS" : "FAIL" }')
  [ "$STATUS" = "PASS" ] || FAIL=1

  echo "exp05b[seed=$SEED]: D_msd=$D_MSD D_vacf=$D_VACF (NIST D*=$NIST_D_STAR, delta=$DELTA_D) first peak r*=$PEAK_R g=$PEAK_G tail g=$TAIL_G -> $STATUS"
  echo "$RHO,$TARGET_T,$SEED,$D_MSD,$D_VACF,$NIST_D_STAR,$DELTA_D,$PEAK_R,$PEAK_G,$TAIL_G,$STATUS" >> "$ROWS"
done

N_PASS=$(awk -F, '$11=="PASS" { n++ } END { print n+0 }' "$ROWS")
N_SEEDS=${#SEEDS[@]}

{
  echo "# purpose: characterise the state experiment 05 actually reaches at its one"
  echo "#       failing density -- see docs/findings/2026-09-21-rho090-outlier.md"
  echo "# git_sha: $GIT_SHA"
  echo "# host: $HOST"
  echo "# units: reduced Lennard-Jones units (sigma = epsilon = m = 1)"
  echo "# state_point: rho*=$RHO, T*=$TARGET_T (NIST's own quoted T* for this row)"
  echo "# command_template: $BIN --cells $CELLS --density $RHO --temperature $TARGET_T"
  echo "#       --cutoff $CUTOFF --skin $SKIN --dt $DT --equilibrate $EQUILIBRATE"
  echo "#       --production $PRODUCTION --thermostat none --seed SEED"
  echo "#       --rdf-out results/exp05b_rdf_rho0.900_seedSEED.csv"
  echo "#       --msd-out results/exp05b_msd_rho0.900_seedSEED.csv"
  echo "# seeds: ${SEEDS[*]} (independent replicas; $PRIMARY_SEED is exp05's own seed for this row)"
  echo "# primary_seed: $PRIMARY_SEED (the trajectory plotted in the figure)"
  echo "# nist_D_star: $NIST_D_STAR (data/nist/README.md, same rho*/T*, fully-equilibrated liquid)"
  echo "# d_floor: $D_FLOOR (a priori 'not frozen' threshold, see script header -- unchanged)"
  echo "# tolerance: D_msd > $D_FLOOR, every seed"
  echo "# outcome: $N_PASS of $N_SEEDS seeds cleared the floor"
  echo "rho_star,T_star,seed,D_msd,D_vacf,nist_D_star,delta_D_msd,rdf_first_peak_r,rdf_first_peak_g,rdf_tail_g,status"
  cat "$ROWS"
} > "$SUMMARY_CSV"

echo "exp05b: wrote $SUMMARY_CSV"
echo "exp05b: $N_PASS of $N_SEEDS seeds cleared the D_msd > $D_FLOOR floor"

if [ "$FAIL" -ne 0 ]; then
  echo "exp05b: FAIL -- at least one seed did not clear the not-frozen floor $D_FLOOR; see docs/findings/2026-09-21-rho090-outlier.md" >&2
  exit 1
fi

echo "exp05b: PASS -- every seed cleared the not-frozen floor $D_FLOOR"
