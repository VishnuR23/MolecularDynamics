#!/usr/bin/env bash
# Experiment 01: the 12-case NIST Standard Reference LJ energy/virial table.
#
# moldyn_run's --nist-config mode already computes every value and exits
# non-zero if any of the 12 cases falls outside its reference's own
# precision, so that exit code IS the pass/fail threshold this script
# asserts. This script's own job is only to capture that table into a
# committed CSV.
#
# The threshold is per value, not flat: NIST publishes every number to
# exactly 5 significant figures, so the true value lies within half a unit
# of the last printed digit, and that half-width is how precisely the
# reference is known. It ranges from 5.34e-6 relative (W = -935.78,
# config 3 / LFS / rc=3, the tightest) to 4.90e-5 (U = -1021.0). The flat
# 5e-5 this replaced was up to 9.4x looser than the reference's actual
# precision -- loose enough that a regression breaking the fifth published
# digit would still have passed. Every case now has to match every digit
# NIST prints. See apps/moldyn_run.cpp (roundingHalfWidth) and
# tests/test_nist_lj.cpp.
# No plotting here -- see plot_exp01.py.
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

# `git describe --always --dirty`, not `rev-parse --short HEAD`: a run
# from a modified working tree is stamped `<sha>-dirty`, so it cannot
# masquerade as the clean commit it was derived from.
GIT_SHA="$(git -C "$ROOT" describe --always --dirty 2>/dev/null || echo unknown)"
HOST="$(uname -s) $(uname -m)"

{
  echo "# command: $CMD"
  echo "# git_sha: $GIT_SHA"
  echo "# host: $HOST"
  echo "# units: reduced Lennard-Jones units (sigma = epsilon = m = 1)"
  echo "# tolerance: per-value, half a unit in NIST's own 5th published significant figure"
  echo "#            (the *_tol columns; *_dev and *_tol are both relative, so dev/tol is the"
  echo "#            fraction of the reference's own precision used). The ulrc columns for the"
  echo "#            LFS rows are absolute instead: NIST's U_LRC there is an exact structural"
  echo "#            zero, not a rounded measurement, so it is checked against 1e-9 absolute."
  echo "# source: NIST Standard Reference Simulation Website, Lennard-Jones fluid"
  echo "#         reference calculations, cuboid cell -- see data/nist/README.md"
  echo "cfg,trunc,rc,u_calc,u_ref,u_dev,u_tol,w_calc,w_ref,w_dev,w_tol,ulrc_calc,ulrc_ref,ulrc_dev,ulrc_tol,status"
  awk 'NF==16 && $1 ~ /^[0-9]+$/ {
        OFS=",";
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16
      }' "$RAW_TXT"
} > "$OUT_CSV"

echo "exp01: wrote $OUT_CSV"

if [ "$STATUS" -ne 0 ]; then
  echo "exp01: FAIL -- moldyn_run reported one or more cases outside NIST's own 5-significant-figure rounding half-width" >&2
  exit 1
fi

echo "exp01: PASS -- all 12 NIST cases match every digit NIST publishes"
