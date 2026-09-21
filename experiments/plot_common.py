"""Shared plotting infrastructure for the moldyn reproduction experiments.

Every plot_expNN.py script uses apply_style(), load_csv() and save() from
here rather than styling itself, so the figures read as one system. See
.superpowers/sdd/2026-09-20-phase1-lennard-jones/task-15-brief.md and the
dataviz skill for the design rationale.

Only numpy and matplotlib are used anywhere in experiments/ (no seaborn, no
pandas) -- see experiments/requirements.txt.
"""

from __future__ import annotations

import pathlib
import subprocess

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

ROOT = pathlib.Path(__file__).resolve().parent.parent

# Colors: slot 1 (blue) is always "our simulation", slot 2 (orange) is
# always "the reference/literature value". This order is fixed across every
# figure in this repository, per the dataviz skill's rule that categorical
# hues are assigned by identity, never cycled -- and validated as a CVD-safe
# adjacent pair (worst adjacent CVD Delta E 9.1, normal-vision 19.6).
OURS = "#2a78d6"
REFERENCE = "#eb6834"
BAD = "#e34948"
GOOD = "#0ca30c"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#d8d7d2"
SURFACE = "#fcfcfb"


def apply_style() -> None:
    """Set the shared figure style once. Call before creating any figure."""
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.size": 11,
            "font.family": "sans-serif",
            "text.color": TEXT_PRIMARY,
            "axes.edgecolor": GRID,
            "axes.labelcolor": TEXT_PRIMARY,
            "axes.titlecolor": TEXT_PRIMARY,
            "axes.titleweight": "bold",
            "axes.titlesize": 13,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.7,
            "grid.alpha": 0.8,
            "xtick.color": TEXT_SECONDARY,
            "ytick.color": TEXT_SECONDARY,
            "legend.frameon": False,
            "legend.fontsize": 9.5,
            "lines.linewidth": 2.0,
            "lines.markersize": 6,
            "figure.dpi": 150,
            "savefig.dpi": 200,
        }
    )


def load_csv(path):
    """Load one of this project's provenance-stamped CSVs.

    Lines beginning with '#' are provenance/metadata: 'key: value' lines are
    collected into a header dict (later duplicate keys are suffixed _2, _3,
    ... so multiple '# command: ...' lines are all kept); the first
    non-comment line is the CSV column-name row. Returns (header, data),
    where data is a numpy structured array (mixed-type columns, e.g. a
    string 'phase' or 'status' column, come through fine) indexed by field
    name, e.g. data['r'], data['g'].
    """
    import io

    header: dict[str, str] = {}
    seen: dict[str, int] = {}
    body_lines: list[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                content = line[1:].strip()
                if ":" in content:
                    key, value = content.split(":", 1)
                    key, value = key.strip(), value.strip()
                    if key in seen:
                        seen[key] += 1
                        key = f"{key}_{seen[key]}"
                    else:
                        seen[key] = 1
                    header[key] = value
            else:
                body_lines.append(line)

    # Comment lines are stripped by hand above (rather than left to
    # genfromtxt's own comments= handling) because a provenance line like
    # "# source: ... see data/nist/README.md" can otherwise confuse
    # genfromtxt's column-count auto-detection pre-pass.
    data = np.genfromtxt(io.StringIO("".join(body_lines)), delimiter=",", names=True, dtype=None,
                          encoding="utf-8")
    return header, data


def _git_sha() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=ROOT,
            capture_output=True,
            check=True,
            text=True,
        )
        return out.stdout.strip()
    except Exception:
        return "unknown"


def save(fig, path) -> None:
    """Stamp the git SHA into the figure footer and write the PNG."""
    sha = _git_sha()
    fig.text(
        0.995,
        0.005,
        f"moldyn @ {sha}",
        ha="right",
        va="bottom",
        fontsize=7.5,
        color=TEXT_SECONDARY,
        alpha=0.75,
    )
    fig.savefig(path)
    plt.close(fig)
