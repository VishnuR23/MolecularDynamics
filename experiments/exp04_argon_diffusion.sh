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
# result. No plotting here -- see plot_exp04.py.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT/build/apps/moldyn_run"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

GIT_SHA="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
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
  echo "# production_steps: $PRODUCTION (NVE)"
  echo "# n_replicas_per_state_point: $N_REPLICAS"
  echo "# error_bars: standard error of the mean across independent replicas (see script header)"
  echo "# tolerance: |mean(D_msd) - mean(D_vacf)| <= se(D_msd) + se(D_vacf), per state point"
  echo "label,rho_star,T_star,seed,D_msd,D_vacf"
  cat "$ROWS"
} > "$SUMMARY_CSV"

echo "exp04: wrote $SUMMARY_CSV"

check_label() {
  local label="$1"
  awk -F, -v label="$label" '
    $1==label { n++; sum_msd+=$5; sum_vacf+=$6; msd[n]=$5; vacf[n]=$6 }
    END {
      if (n<2) { print "insufficient replicas"; exit 1 }
      mean_msd=sum_msd/n; mean_vacf=sum_vacf/n
      for (i=1;i<=n;i++) { ss_msd+=(msd[i]-mean_msd)^2; ss_vacf+=(vacf[i]-mean_vacf)^2 }
      se_msd=sqrt(ss_msd/(n-1))/sqrt(n)
      se_vacf=sqrt(ss_vacf/(n-1))/sqrt(n)
      diff = (mean_msd-mean_vacf<0) ? mean_vacf-mean_msd : mean_msd-mean_vacf
      combined = se_msd+se_vacf
      printf "mean_D_msd=%.6f se_D_msd=%.6f mean_D_vacf=%.6f se_D_vacf=%.6f |diff|=%.6f combined_se=%.6f %s\n", \
             mean_msd, se_msd, mean_vacf, se_vacf, diff, combined, (diff<=combined ? "PASS" : "FAIL")
      exit (diff<=combined ? 0 : 1)
    }
  ' "$ROWS"
}

FAIL=0

RAHMAN_REPORT=$(check_label rahman) || FAIL=1
echo "exp04[rahman]: $RAHMAN_REPORT"

VERLET_REPORT=$(check_label verlet) || FAIL=1
echo "exp04[verlet]: $VERLET_REPORT"

exit $FAIL
