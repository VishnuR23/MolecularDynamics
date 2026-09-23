"""Shared plotting infrastructure for the moldyn reproduction experiments.

Every plot_expNN.py script uses apply_style(), load_csv() and save() from
here rather than styling itself, so the figures read as one system: one
palette in which blue is always "ours" and orange always "the reference",
one axis per panel (never a dual-y chart), a legend on every multi-series
figure, and the git revision stamped into every footer.

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
    # A single-data-row CSV genfromtxt's to a 0-d structured array (a bare
    # record, not a length-1 array), which breaks any caller that does
    # data[0] or data['col'][i]. Always return a 1-d array so every caller
    # can index it uniformly regardless of row count.
    data = np.atleast_1d(data)
    return header, data


def _git_sha() -> str:
    """The revision this figure was generated from.

    The revision of the CODE, so a figure rendered from a modified working
    tree is stamped `<sha>-dirty` and cannot masquerade as the clean commit
    it came from. Dirtiness ignores results/: the experiments write into the
    repository, so regenerating one would otherwise mark every later figure
    dirty for a reason unrelated to the code. Matches experiments/git_sha.sh
    and apps/GitSha.cmake, so all three agree.
    """
    try:
        sha = subprocess.run(
            ["git", "describe", "--always"],
            cwd=ROOT, capture_output=True, check=True, text=True,
        ).stdout.strip()
        if not sha:
            return "unknown"
        dirt = subprocess.run(
            ["git", "status", "--porcelain", "--", ".", ":(exclude)results"],
            cwd=ROOT, capture_output=True, check=True, text=True,
        ).stdout.strip()
        return f"{sha}-dirty" if dirt else sha
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
