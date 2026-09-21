#!/usr/bin/env bash
# Experiment 01: the 12-case NIST Standard Reference LJ energy/virial table.
#
# moldyn_run's --nist-config mode already computes every value and exits
# non-zero if any of the 12 cases exceeds 5e-5 relative deviation, so that
# exit code IS the pass/fail threshold this script asserts. This script's
# own job is only to capture that table into a committed CSV. No plotting
# here -- see plot_exp01.py.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$ROOT/build/apps/moldyn_run"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

CONFIG_DIR="$ROOT/data/nist/lj"
CMD="$BIN --nist-config $CONFIG_DIR"
OUT_CSV="$RESULTS/exp01_nist_table.csv"
RAW_TXT="$(mktemp)"
trap 'rm -f "$RAW_TXT"' EXIT

echo "exp01: $CMD"
set +e
$CMD >"$RAW_TXT"
STATUS=$?
set -e
cat "$RAW_TXT"

GIT_SHA="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
HOST="$(uname -s) $(uname -m)"

{
  echo "# command: $CMD"
  echo "# git_sha: $GIT_SHA"
  echo "# host: $HOST"
  echo "# units: reduced Lennard-Jones units (sigma = epsilon = m = 1)"
  echo "# tolerance: 5e-5 relative deviation (moldyn_run's own pass/fail threshold)"
  echo "# source: NIST Standard Reference Simulation Website, Lennard-Jones fluid"
  echo "#         reference calculations, cuboid cell -- see data/nist/README.md"
  echo "cfg,trunc,rc,u_calc,u_ref,u_dev,w_calc,w_ref,w_dev,ulrc_calc,ulrc_ref,ulrc_dev,status"
  awk 'NF==13 && $1 ~ /^[0-9]+$/ {
        OFS=",";
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13
      }' "$RAW_TXT"
} > "$OUT_CSV"

echo "exp01: wrote $OUT_CSV"

if [ "$STATUS" -ne 0 ]; then
  echo "exp01: FAIL -- moldyn_run reported one or more cases exceeding 5e-5 relative deviation" >&2
  exit 1
fi

echo "exp01: PASS -- all 12 NIST cases within 5e-5 relative deviation"
