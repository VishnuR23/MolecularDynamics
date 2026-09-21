#!/usr/bin/env bash
# Experiment 03: the radial distribution function g(r) of liquid argon at
# two historic MD state points, both in Lennard-Jones reduced units
# (sigma=3.405 A, eps/kB=119.8 K, M=39.948 g/mol):
#
#   Rahman (1964), Phys. Rev. 136, A405: 94.4 K, 1.374 g/cm^3
#     -> rho* = 0.8177, T* = 0.7880   (primary state point)
#   Verlet (1967), Phys. Rev. 159, 98: the state point most widely
#     tabulated afterwards as "the" LJ argon liquid point
#     -> rho* = 0.8442, T* = 0.7280   (secondary state point)
#
# The task brief's plan mislabelled 0.8442/0.7280 as Rahman's; that is
# actually Verlet's. Both are run here and labelled correctly -- see
# data/nist/README.md is NOT the source for these (they are not from the
# NIST SRSW; see the reduction above and the task-15 report for how the
# reduced values were obtained from the physical (T, rho) each paper
# reports).
#
# N=864 (fcc, 6 cells/side), a standard choice in this literature. Asserts
# the g(r) first-peak position and height land in a band drawn from the
# standard LJ-argon-near-triple-point literature (e.g. Rahman 1964 fig. 2;
# Verlet 1967; Hansen & McDonald, "Theory of Simple Liquids"): first peak
# near r* ~ 1.0-1.25, height ~ 2.3-3.6. This is a qualitative literature
# band, not an exact tabulated NIST number -- unlike experiments 01 and 05,
# no source publishes g(r) to enough digits for a tight numeric check.
# No plotting here -- see plot_exp03.py.
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
PRODUCTION=20000

PEAK_LO=0.85
PEAK_HI=1.40
PEAK_POS_MIN=1.00
PEAK_POS_MAX=1.25
PEAK_HEIGHT_MIN=2.3
PEAK_HEIGHT_MAX=3.6

run_one() {
  local label="$1" rho="$2" temp="$3" seed="$4"
  local rdf_csv="$RESULTS/exp03_rdf_${label}.csv"
  local cmd="$BIN --cells $CELLS --density $rho --temperature $temp --cutoff $CUTOFF \
--skin $SKIN --dt $DT --equilibrate $EQUILIBRATE --production $PRODUCTION \
--thermostat none --seed $seed --rdf-out $rdf_csv"
  echo "exp03[$label]: $cmd"
  eval "$cmd"

  # Locate the first peak: the maximum g(r) with r in [PEAK_LO, PEAK_HI],
  # a window chosen to contain argon's first coordination shell and
  # exclude the second shell near r* ~ 2.
  read -r peak_r peak_g <<EOF
$(awk -F, -v lo="$PEAK_LO" -v hi="$PEAK_HI" '
    /^#/ { next }
    !seen_header { seen_header=1; next }
    { r=$1; g=$2; if (r>=lo && r<=hi && g>maxg) { maxg=g; maxr=r } }
    END { printf "%.6f %.6f", maxr, maxg }
  ' "$rdf_csv")
EOF

  echo "exp03[$label]: first peak at r*=$peak_r, g=$peak_g"
  echo "$label,$rho,$temp,$seed,$peak_r,$peak_g,$rdf_csv" >> "$SUMMARY_ROWS"
}

SUMMARY_CSV="$RESULTS/exp03_argon_rdf.csv"
SUMMARY_ROWS="$(mktemp)"
trap 'rm -f "$SUMMARY_ROWS"' EXIT

run_one rahman 0.8177 0.7880 3001
run_one verlet 0.8442 0.7280 3002

{
  echo "# state_point_rahman: rho*=0.8177 T*=0.7880 (Rahman 1964, Phys. Rev. 136, A405; 94.4 K, 1.374 g/cm3, sigma=3.405 A, eps/kB=119.8 K, M=39.948 g/mol)"
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
  echo "# peak_search_window: r* in [$PEAK_LO, $PEAK_HI]"
  echo "# tolerance: first-peak r* in [$PEAK_POS_MIN, $PEAK_POS_MAX], height in [$PEAK_HEIGHT_MIN, $PEAK_HEIGHT_MAX]"
  echo "# reference: qualitative literature band, not an exact tabulated value (Rahman 1964; Verlet 1967; Hansen & McDonald, Theory of Simple Liquids)"
  echo "label,rho_star,T_star,seed,peak_r_star,peak_g,rdf_csv"
  cat "$SUMMARY_ROWS"
} > "$SUMMARY_CSV"

echo "exp03: wrote $SUMMARY_CSV"

FAIL=0
while IFS=, read -r label rho temp seed peak_r peak_g rdf_csv; do
  ok=$(awk -v r="$peak_r" -v g="$peak_g" -v rmin="$PEAK_POS_MIN" -v rmax="$PEAK_POS_MAX" \
            -v gmin="$PEAK_HEIGHT_MIN" -v gmax="$PEAK_HEIGHT_MAX" \
            'BEGIN { print (r>=rmin && r<=rmax && g>=gmin && g<=gmax) ? "yes" : "no" }')
  if [ "$ok" != "yes" ]; then
    echo "exp03[$label]: FAIL -- peak r*=$peak_r g=$peak_g outside [$PEAK_POS_MIN,$PEAK_POS_MAX] x [$PEAK_HEIGHT_MIN,$PEAK_HEIGHT_MAX]" >&2
    FAIL=1
  else
    echo "exp03[$label]: PASS -- peak r*=$peak_r g=$peak_g within literature band"
  fi
done < "$SUMMARY_ROWS"

exit $FAIL
