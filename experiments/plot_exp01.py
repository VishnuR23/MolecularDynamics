#!/usr/bin/env python3
"""Experiment 01 plot: how much of NIST's own published precision each of
our 12 LJ energy/virial/tail-correction values uses.

NIST publishes every reference value to 5 significant figures, so the true
value lies within half a unit of the last printed digit. That per-value
rounding half-width -- not a flat tolerance -- is the bound plotted here,
and each bar is the ratio deviation/bound: 1.0 means "just inside the last
digit NIST printed", and anything below 1.0 means our value agrees with
every digit NIST publishes.

Reads results/exp01_nist_table.csv, writes results/exp01_nist_table.png.
No simulation here -- see exp01_nist_table.sh.
"""
import numpy as np

from plot_common import BAD, GOOD, ROOT, TEXT_SECONDARY, apply_style, load_csv, save

OUT_CSV = ROOT / "results" / "exp01_nist_table.csv"
OUT_PNG = ROOT / "results" / "exp01_nist_table.png"

# Ratios below this are drawn at the floor so a log axis can show them.
FLOOR = 1e-3


def ratios(dev, tol):
    """deviation / per-value bound, floored so a log axis can render zeros."""
    out = np.array([d / t if t > 0 else 0.0 for d, t in zip(dev, tol)])
    return np.where(out > FLOOR, out, FLOOR)


def main() -> None:
    apply_style()
    header, data = load_csv(OUT_CSV)

    n = len(data)
    labels = [f"{int(c)}/{t}/{rc:g}" for c, t, rc in zip(data["cfg"], data["trunc"], data["rc"])]

    u = ratios(data["u_dev"].astype(float), data["u_tol"].astype(float))
    w = ratios(data["w_dev"].astype(float), data["w_tol"].astype(float))
    ulrc = ratios(data["ulrc_dev"].astype(float), data["ulrc_tol"].astype(float))

    x = np.arange(n)
    width = 0.26

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.set_yscale("log")
    ax.bar(x - width, u, width, label="U (pair energy)", color="#2a78d6")
    ax.bar(x, w, width, label="W (virial)", color="#1baf7a")
    ax.bar(x + width, ulrc, width, label="U_LRC (tail correction)", color="#eda100")
    ax.axhline(
        1.0,
        color=BAD,
        linestyle="--",
        linewidth=1.5,
        label="NIST's own rounding half-width\n(half a unit in the 5th published digit)",
    )

    worst = max(u.max(), w.max(), ulrc.max())
    worst_i = int(np.argmax(np.maximum(np.maximum(u, w), ulrc)))
    ax.annotate(
        f"worst case {worst:.3f}x of the bound",
        (x[worst_i], worst),
        textcoords="offset points",
        xytext=(0, 7),
        ha="left",
        fontsize=8.5,
        color=TEXT_SECONDARY,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("NIST case: config / truncation / $r_c^*$")
    ax.set_ylabel("deviation / NIST's own rounding half-width")
    ax.set_title("Experiment 01 — LJ energy, virial and tail correction vs. the 12 NIST cases",
                 pad=34)
    all_pass = all(s == "PASS" for s in data["status"])
    ax.text(
        0.01,
        1.015,
        (
            "ALL 12 CASES MATCH EVERY DIGIT NIST PUBLISHES"
            if all_pass
            else "ONE OR MORE CASES FALL OUTSIDE NIST'S PUBLISHED PRECISION"
        ),
        transform=ax.transAxes,
        fontsize=10,
        color=GOOD if all_pass else BAD,
        fontweight="bold",
    )
    ax.text(
        0.01,
        -0.30,
        f"Bars below 1.0 agree with NIST to every published digit. Values at {FLOOR:g} are exact "
        "to the precision plotted (the four LinearForceShift tail corrections are structurally zero).",
        transform=ax.transAxes,
        fontsize=8,
        color=TEXT_SECONDARY,
    )
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()

    save(fig, OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
