#!/usr/bin/env bash
# Runs every experiment (moldyn_run + threshold assertions), then every
# plot script, in order.
#
# Experiments 01 and 02 check the *engine* (the NIST energy/virial table,
# the integrator's convergence order): a failure there means something is
# broken, so this script stops immediately, matching CI.
#
# Experiments 03-05 check *physics reproduction* against the literature.
# The standing rule for this repository is that a computed number which
# disagrees with the literature is a finding to investigate and report,
# not a tolerance to widen. A failing threshold there is a result, not a
# bug, so this script keeps going (every experiment still gets its CSV and
# PNG) but remembers that a disagreement was found and exits non-zero at
# the end, so `make reproduce`'s own exit status is still honest.
#
# Two experiments currently disagree, and both are written up:
#   docs/findings/2026-09-21-exp04-einstein-vs-green-kubo.md
#   docs/findings/2026-09-21-rho090-outlier.md
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PYTHON="$ROOT/.venv/bin/python3"

if [ ! -x "$PYTHON" ]; then
  cat >&2 <<MSG
run_all: $PYTHON not found.

The plotting virtualenv is gitignored, so a fresh clone does not have one.
'make reproduce' builds it for you; if you are running this script directly,
create it first with exactly:

  make venv

or, by hand:

  uv venv .venv && uv pip install --python .venv/bin/python3 -r experiments/requirements.txt

or, without uv:

  python3 -m venv .venv && .venv/bin/python3 -m pip install -r experiments/requirements.txt
MSG
  exit 1
fi

ENGINE_EXPERIMENTS="exp01_nist_table exp02_energy_conservation"
PHYSICS_EXPERIMENTS="exp03_argon_rdf exp04_argon_diffusion exp05_lj_eos exp05b_rho090_diagnostic"

START=$(date +%s)

for exp in $ENGINE_EXPERIMENTS; do
  tag="${exp%%_*}"
  echo "==== $exp: simulating ===="
  bash "$SCRIPT_DIR/${exp}.sh" || { echo "run_all: $exp FAILED -- this indicates a broken engine, stopping" >&2; exit 1; }
  echo "==== $exp: plotting ===="
  "$PYTHON" "$SCRIPT_DIR/plot_${tag}.py"
done

FOUND_DISAGREEMENT=0
for exp in $PHYSICS_EXPERIMENTS; do
  tag="${exp%%_*}"
  echo "==== $exp: simulating ===="
  if ! bash "$SCRIPT_DIR/${exp}.sh"; then
    echo "run_all: $exp reported a literature disagreement -- see its output above and docs/findings/" >&2
    FOUND_DISAGREEMENT=1
  fi
  echo "==== $exp: plotting ===="
  "$PYTHON" "$SCRIPT_DIR/plot_${tag}.py"
done

END=$(date +%s)
echo "==== run_all: all experiments complete in $((END-START))s ===="

if [ "$FOUND_DISAGREEMENT" -ne 0 ]; then
  echo "run_all: one or more physics experiments (03-05) disagreed with the literature outside tolerance -- see docs/findings/ for the investigations" >&2
  exit 1
fi
