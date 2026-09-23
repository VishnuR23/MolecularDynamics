#!/usr/bin/env bash
# Experiment 04: the argon self-diffusion coefficient at the same two state
# points as experiment 03 (Rahman 1964 and Verlet 1967, in reduced units --
# see exp03_argon_rdf.sh for the derivation), NVE production after Langevin
# equilibration, estimated two independent ways: Einstein (slope of MSD)
# and Green-Kubo (integral of the VACF).
#
# moldyn_run computes both from a single trajectory (see its --msd-out
# header), but a single trajectory gives no error bar. So each state point
# is run as N_REPLICAS independent trajectories (distinct seeds), and the
# error bar on D is the standard error of the mean across those replicas --
# this is the Flyvbjerg-Petersen block-averaging idea (src/moldyn/analyze/
# block_average.*) taken to its coarsest, fully-decorrelated block: whole
# independent runs, rather than blocks cut from one run. Asserts the two
# estimates (D_msd, D_vacf) agree within their combined standard errors, at
# both state points. A diffusion coefficient without an error bar is not a
# result.
#
# TWO separate claims are asserted here, and they are independent:
#
#   1. INTERNAL CONSISTENCY -- the Einstein and Green-Kubo estimates of the
#      same D, from the same trajectories, agree within their combined
#      standard errors. This is a check on our own estimators.
#      OUTCOME: the Verlet point passes; the Rahman point FAILS, and has
#      from the first run. See docs/findings/2026-09-21-exp04-einstein-vs-
#      green-kubo.md for the investigation, including what the
#      centre-of-mass fix did and did not explain. This script's non-zero
#      exit is that finding, not a broken engine.
#
#   2. AGREEMENT WITH EXPERIMENT -- mean(D_msd) at the Rahman state point
#      against Rahman's own 1964 laboratory measurement of liquid argon,
#      D = 2.43e-5 cm^2/s at 94.4 K and 1.374 g/cm^3. Reduced with the
#      same constants exp03_argon_rdf.sh uses (sigma = 3.405 A,
#      eps/kB = 119.8 K, M = 39.948 g/mol), the diffusion unit is
#      sigma*sqrt(eps/m) = 5.3767e-4 cm^2/s, so
#      D*_exp = 2.43e-5 / 5.3767e-4 = 0.0452.
#      Spec section 6 asks that D agree "with each other AND with
#      experiment"; only the first half used to be tested.
#      Tolerance: 10% relative, chosen a priori and before the
#      centre-of-mass fix was run. That is the level at which LJ argon
#      transport agreement is meaningful: sigma and eps are themselves
#      fitted parameters good to ~1%, D near the triple point is steeply
#      sensitive to both, and Rahman's number is a single 1964 measurement
#      quoted to three figures. A tighter bound would be asserting
#      precision the comparison does not have; a looser one would not be
#      a test.
#
# No plotting here -- see plot_exp04.py.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT/build/apps/moldyn_run"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

. "$(dirname "${BASH_SOURCE[0]}")/git_sha.sh"
GIT_SHA="$(moldyn_git_sha "$ROOT")"
HOST="$(uname -s) $(uname -m)"

CELLS=6
CUTOFF=2.5
SKIN=0.3
DT=0.005
EQUILIBRATE=5000
# Msd::push/Vacf::push are O(maxLag * natoms) per production step (they
# accumulate every valid time origin at every lag, every step), and
# moldyn_run sets maxLag = min(production/2, 2000). At N=864 that made a
# single 20000-step production run take ~2.5 minutes; with 8 replicas
# total that is too slow for a script meant to be re-run. 2000 steps
# (maxLag=1000, fit window lags [250,751) i.e. t* in [1.25, 3.76]) is
# short enough to keep exp04 to a couple of minutes total while still
# giving each lag in the fit window O(10^6) origin-pairs of statistics --
# the same MSD/VACF methodology already
# validated against the exact Ornstein-Uhlenbeck result in
# tests/test_analysis.cpp, just run for less physical time per replica,
# which is why this experiment averages over N_REPLICAS independent
# replicas rather than relying on one long trajectory.
PRODUCTION=2000
# 8, not 4: an early run with 4 replicas/point found the two D estimates
# agreeing at the Verlet point but disagreeing by ~2x their combined
# standard error at the Rahman point. With only 4 replicas the standard
# error estimate itself is noisy (relative uncertainty ~1/sqrt(2*(n-1)) =
# 41%), so that could be a real short-trajectory Einstein/Green-Kubo bias
# or just n=4 sampling noise -- not enough replicas to tell. 8 halves that
# uncertainty and is still only ~2-3 minutes total.
N_REPLICAS=8

# Rahman 1964's measured self-diffusion coefficient for liquid argon,
# reduced to LJ units -- see the derivation in this script's header.
D_EXP_STAR=0.0452
D_EXP_TOL_FRAC=0.10

SUMMARY_CSV="$RESULTS/exp04_argon_diffusion.csv"
ROWS="$(mktemp)"
trap 'rm -f "$ROWS"' EXIT

run_replica() {
  local label="$1" rho="$2" temp="$3" seed="$4"
  local msd_csv="$RESULTS/exp04_msd_${label}_seed${seed}.csv"
  local cmd="$BIN --cells $CELLS --density $rho --temperature $temp --cutoff $CUTOFF \
--skin $SKIN --dt $DT --equilibrate $EQUILIBRATE --production $PRODUCTION \
--thermostat none --seed $seed --msd-out $msd_csv"
  echo "exp04[$label seed=$seed]: $cmd"
  eval "$cmd"

  local d_msd d_vacf
  d_msd=$(grep '^# D_msd' "$msd_csv" | sed -E 's/.*: *//')
  d_vacf=$(grep '^# D_vacf' "$msd_csv" | sed -E 's/.*: *//')
  echo "exp04[$label seed=$seed]: D_msd=$d_msd  D_vacf=$d_vacf"
  echo "$label,$rho,$temp,$seed,$d_msd,$d_vacf" >> "$ROWS"
}

RAHMAN_SEEDS=(4001 4002 4003 4004 4005 4006 4007 4008)
VERLET_SEEDS=(4101 4102 4103 4104 4105 4106 4107 4108)

for s in "${RAHMAN_SEEDS[@]}"; do
  run_replica rahman 0.8177 0.7880 "$s"
done
for s in "${VERLET_SEEDS[@]}"; do
  run_replica verlet 0.8442 0.7280 "$s"
done

# --- statistics, computed before the CSV is written so the CSV can carry
# --- each row's verdict (a reader of the data alone should learn the
# --- outcome, not just the tolerance).
stats_for() {
  local label="$1"
  awk -F, -v label="$label" '
    $1==label { n++; sum_msd+=$5; sum_vacf+=$6; msd[n]=$5; vacf[n]=$6 }
    END {
      if (n<2) { exit 1 }
      mean_msd=sum_msd/n; mean_vacf=sum_vacf/n
      for (i=1;i<=n;i++) { ss_msd+=(msd[i]-mean_msd)^2; ss_vacf+=(vacf[i]-mean_vacf)^2 }
      se_msd=sqrt(ss_msd/(n-1))/sqrt(n)
      se_vacf=sqrt(ss_vacf/(n-1))/sqrt(n)
      diff = mean_vacf-mean_msd
      adiff = (diff<0) ? -diff : diff
      combined = se_msd+se_vacf
      printf "%.6f %.6f %.6f %.6f %.6f %.6f %s", \
             mean_msd, se_msd, mean_vacf, se_vacf, diff, combined, \
             (adiff<=combined ? "PASS" : "FAIL")
    }
  ' "$ROWS"
}

read -r R_MEAN_MSD R_SE_MSD R_MEAN_VACF R_SE_VACF R_DIFF R_COMBINED R_STATUS <<EOF
$(stats_for rahman)
EOF
read -r V_MEAN_MSD V_SE_MSD V_MEAN_VACF V_SE_VACF V_DIFF V_COMBINED V_STATUS <<EOF
$(stats_for verlet)
EOF

# Rahman's laboratory value, compared against the Einstein estimate at his
# own state point. Only the Rahman row has an experimental counterpart --
# Verlet 1967 is a simulation paper and publishes no measured D.
read -r EXP_DELTA EXP_PCT EXP_STATUS <<EOF
$(awk -v d="$R_MEAN_MSD" -v e="$D_EXP_STAR" -v tol="$D_EXP_TOL_FRAC" 'BEGIN {
    delta = d - e
    adelta = (delta<0) ? -delta : delta
    pct = 100.0 * delta / e
    printf "%.6f %.3f %s", delta, pct, (adelta <= tol*e ? "PASS" : "FAIL")
  }')
EOF

{
  echo "# state_point_rahman: rho*=0.8177 T*=0.7880 (Rahman 1964, Phys. Rev. 136, A405)"
  echo "# state_point_verlet: rho*=0.8442 T*=0.7280 (Verlet 1967, Phys. Rev. 159, 98)"
  echo "# git_sha: $GIT_SHA"
  echo "# host: $HOST"
  echo "# units: reduced Lennard-Jones units (sigma = epsilon = m = 1)"
  echo "# cells: $CELLS (N=864)"
  echo "# cutoff: $CUTOFF"
  echo "# skin: $SKIN"
  echo "# dt: $DT"
  echo "# equilibrate_steps: $EQUILIBRATE (Langevin)"
  echo "# production_steps: $PRODUCTION (NVE, centre-of-mass momentum re-zeroed at the"
  echo "#       equilibration/production boundary so p_com stays zero through production)"
  echo "# n_replicas_per_state_point: $N_REPLICAS"
  echo "# error_bars: standard error of the mean across independent replicas (see script header)"
  echo "# tolerance_estimators: |mean(D_msd) - mean(D_vacf)| <= se(D_msd) + se(D_vacf), per state point"
  echo "# tolerance_experiment: |mean(D_msd) - D_exp_star| <= ${D_EXP_TOL_FRAC} * D_exp_star, Rahman point only"
  echo "# experimental_reference: D = 2.43e-5 cm^2/s (Rahman 1964, liquid argon, 94.4 K, 1.374 g/cm3)"
  echo "#       reduced by sigma*sqrt(eps/m) = 5.3767e-4 cm^2/s (sigma=3.405 A, eps/kB=119.8 K,"
  echo "#       M=39.948 g/mol) -> D_exp_star = $D_EXP_STAR"
  echo "# result_rahman: mean_D_msd=$R_MEAN_MSD se=$R_SE_MSD mean_D_vacf=$R_MEAN_VACF se=$R_SE_VACF"
  echo "#       D_vacf-D_msd=$R_DIFF combined_se=$R_COMBINED estimators=$R_STATUS"
  echo "# result_verlet: mean_D_msd=$V_MEAN_MSD se=$V_SE_MSD mean_D_vacf=$V_MEAN_VACF se=$V_SE_VACF"
  echo "#       D_vacf-D_msd=$V_DIFF combined_se=$V_COMBINED estimators=$V_STATUS"
  echo "# result_experiment: mean_D_msd=$R_MEAN_MSD vs D_exp_star=$D_EXP_STAR"
  echo "#       delta=$EXP_DELTA (${EXP_PCT}%) experiment=$EXP_STATUS"
  echo "# findings: docs/findings/2026-09-21-exp04-einstein-vs-green-kubo.md"
  echo "label,rho_star,T_star,seed,D_msd,D_vacf,D_exp_star,estimator_status,experiment_status"
  while IFS=, read -r label rho temp seed d_msd d_vacf; do
    if [ "$label" = "rahman" ]; then
      echo "$label,$rho,$temp,$seed,$d_msd,$d_vacf,$D_EXP_STAR,$R_STATUS,$EXP_STATUS"
    else
      echo "$label,$rho,$temp,$seed,$d_msd,$d_vacf,NA,$V_STATUS,NA"
    fi
  done < "$ROWS"
} > "$SUMMARY_CSV"

echo "exp04: wrote $SUMMARY_CSV"

FAIL=0

echo "exp04[rahman]: mean_D_msd=$R_MEAN_MSD se_D_msd=$R_SE_MSD mean_D_vacf=$R_MEAN_VACF se_D_vacf=$R_SE_VACF D_vacf-D_msd=$R_DIFF combined_se=$R_COMBINED $R_STATUS"
[ "$R_STATUS" = "PASS" ] || FAIL=1

echo "exp04[verlet]: mean_D_msd=$V_MEAN_MSD se_D_msd=$V_SE_MSD mean_D_vacf=$V_MEAN_VACF se_D_vacf=$V_SE_VACF D_vacf-D_msd=$V_DIFF combined_se=$V_COMBINED $V_STATUS"
[ "$V_STATUS" = "PASS" ] || FAIL=1

echo "exp04[rahman vs experiment]: mean_D_msd=$R_MEAN_MSD vs Rahman 1964 measured D*=$D_EXP_STAR -- delta=$EXP_DELTA (${EXP_PCT}%), tolerance ${D_EXP_TOL_FRAC} relative: $EXP_STATUS"
if [ "$EXP_STATUS" != "PASS" ]; then
  echo "exp04: FAIL -- the Einstein D does not reproduce Rahman's measured value within ${D_EXP_TOL_FRAC} relative" >&2
  FAIL=1
fi

exit $FAIL
