#!/usr/bin/env bash
# Experiment 02: velocity-Verlet energy-conservation order.
#
# Same setup as tests/test_integrator.cpp's "velocity-Verlet energy error
# scales as dt^2" case: fcc(3), rho*=0.80, T*=1.0, t=1.0 of physical time,
# linear-force-shift truncation (continuous force at rc, so the measured
# drift is the integrator's own discretization error, not the extra,
# non-dt^2 error a plain truncated potential adds when atoms cross rc --
# moldyn_run defaults to plain truncation, so --truncation is set
# explicitly here). Asserts the fitted log-log slope of relative energy
# error vs dt lands in [1.85, 2.15]. No plotting here -- see plot_exp02.py.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT/build/apps/moldyn_run"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

. "$(dirname "${BASH_SOURCE[0]}")/git_sha.sh"
GIT_SHA="$(moldyn_git_sha "$ROOT")"
HOST="$(uname -s) $(uname -m)"
SEED=12345
TOTAL_TIME=1.0

DTS=(0.001 0.002 0.004 0.008)
SUMMARY_CSV="$RESULTS/exp02_energy_conservation.csv"
COMMANDS_FILE="$(mktemp)"
DATA_FILE="$(mktemp)"
trap 'rm -f "$COMMANDS_FILE" "$DATA_FILE"' EXIT

for dt in "${DTS[@]}"; do
  steps=$(awk -v dt="$dt" -v t="$TOTAL_TIME" 'BEGIN { printf "%d", (t/dt)+0.5 }')
  thermo_csv="$RESULTS/exp02_thermo_dt${dt}.csv"
  cmd="$BIN --cells 3 --density 0.80 --temperature 1.0 --cutoff 2.5 --skin 0.4 \
--dt $dt --equilibrate 0 --production $steps --thermostat none --seed $SEED \
--truncation linear-force-shift --thermo-out $thermo_csv"
  echo "exp02: $cmd"
  eval "$cmd"
  echo "$cmd" >> "$COMMANDS_FILE"

  # max |E(t) - E(0)| / |E(0)| over the whole trajectory, from the
  # step,time,phase,potential_energy,kinetic_energy,total_energy,... CSV.
  awk -F, -v dt="$dt" '
    /^#/ { next }
    !seen_header { seen_header=1; next }
    {
      row++
      e = $6
      if (row == 1) e0 = e
      dev = (e - e0 < 0) ? e0 - e : e - e0
      if (dev > maxdev) maxdev = dev
    }
    END { printf "%s,%.10e,%.10e,%.10e\n", dt, e0, maxdev, maxdev / ((e0 < 0) ? -e0 : e0) }
  ' "$thermo_csv" >> "$DATA_FILE"
done

# Least-squares slope of log(relative error) vs log(dt), matching
# tests/test_integrator.cpp exactly.
SLOPE=$(awk -F, '
  { n++; ldt=log($1); lerr=log($4); sdt+=ldt; serr+=lerr; ldts[n]=ldt; lerrs[n]=lerr }
  END {
    mdt=sdt/n; merr=serr/n;
    for (i=1;i<=n;i++) { num += (ldts[i]-mdt)*(lerrs[i]-merr); den += (ldts[i]-mdt)^2 }
    printf "%.6f", num/den
  }
' "$DATA_FILE")

{
  n=0
  while IFS= read -r line; do
    n=$((n+1))
    echo "# command_$n: $line"
  done < "$COMMANDS_FILE"
  echo "# git_sha: $GIT_SHA"
  echo "# host: $HOST"
  echo "# units: reduced Lennard-Jones units (sigma = epsilon = m = 1)"
  echo "# seed: $SEED"
  echo "# total_time: $TOTAL_TIME"
  echo "# fitted_slope: $SLOPE"
  echo "# expected_slope: 2.0 (velocity-Verlet is a second-order symplectic integrator)"
  echo "# tolerance: slope in [1.85, 2.15]"
  echo "dt,e0,max_abs_energy_drift,relative_energy_error"
  cat "$DATA_FILE"
} > "$SUMMARY_CSV"

echo "exp02: fitted slope = $SLOPE"
echo "exp02: wrote $SUMMARY_CSV"

PASS=$(awk -v s="$SLOPE" 'BEGIN { print (s>1.85 && s<2.15) ? "yes" : "no" }')
if [ "$PASS" != "yes" ]; then
  echo "exp02: FAIL -- fitted slope $SLOPE is outside [1.85, 2.15]" >&2
  exit 1
fi

echo "exp02: PASS -- fitted slope $SLOPE is within [1.85, 2.15]"
