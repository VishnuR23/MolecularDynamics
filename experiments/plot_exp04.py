#!/usr/bin/env python3
"""Experiment 04 plot: D from MSD (Einstein) vs. D from VACF (Green-Kubo),
each with a standard-error bar across independent replicas, at both argon
state points -- and, at the Rahman point, against Rahman's own 1964
laboratory measurement of liquid argon reduced to LJ units.

That third mark is the point of the whole experiment: a self-diffusion
coefficient measured in a laboratory in 1964, reproduced from nothing but
a pair potential and Newton's equations. Reads
results/exp04_argon_diffusion.csv, writes
results/exp04_argon_diffusion.png. No simulation here -- see
exp04_argon_diffusion.sh.
"""
import numpy as np

from plot_common import (
    BAD,
    GOOD,
    OURS,
    REFERENCE,
    ROOT,
    TEXT_SECONDARY,
    apply_style,
    load_csv,
    save,
)

IN_CSV = ROOT / "results" / "exp04_argon_diffusion.csv"
OUT_PNG = ROOT / "results" / "exp04_argon_diffusion.png"

# A third identity beyond plot_common's "ours" (blue) and "reference"
# (orange): a laboratory measurement, not a tabulated simulation result.
EXPERIMENT = "#6b3fa0"

STATE_POINTS = [
    ("rahman", "Rahman (1964)\n$\\rho^*$=0.8177, $T^*$=0.7880"),
    ("verlet", "Verlet (1967)\n$\\rho^*$=0.8442, $T^*$=0.7280"),
]


def mean_se(x: np.ndarray):
    n = len(x)
    mean = x.mean()
    se = x.std(ddof=1) / np.sqrt(n)
    return mean, se


def main() -> None:
    apply_style()
    header, data = load_csv(IN_CSV)
    n_replicas = header.get("n_replicas_per_state_point", "?")

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    x = np.arange(len(STATE_POINTS))
    width = 0.18

    all_top = []
    for i, (label, _) in enumerate(STATE_POINTS):
        rows = data[data["label"] == label]
        d_msd, d_vacf = rows["D_msd"].astype(float), rows["D_vacf"].astype(float)
        m_msd, se_msd = mean_se(d_msd)
        m_vacf, se_vacf = mean_se(d_vacf)
        agree = abs(m_msd - m_vacf) <= se_msd + se_vacf

        # Rahman 1964 measured D for real liquid argon at his own state
        # point; Verlet 1967 is a simulation paper and has no measured
        # counterpart, so the CSV carries "NA" there.
        d_exp_raw = str(rows["D_exp_star"][0])
        d_exp = float(d_exp_raw) if d_exp_raw not in ("NA", "nan", "") else None
        if d_exp is not None:
            ax.hlines(d_exp, x[i] - 0.34, x[i] + 0.34, color=EXPERIMENT, linestyle="-",
                      linewidth=2.5,
                      label="measured (Rahman 1964, $2.43\\times10^{-5}$ cm$^2$/s)"
                      if i == 0 else None)
            pct = 100.0 * (m_msd - d_exp) / d_exp
            ax.annotate(
                f"$D^*_{{exp}}$ = {d_exp:.4f}\nEinstein $D^*$ is {pct:+.1f}% from it",
                (x[i] + 0.34, d_exp),
                textcoords="offset points",
                xytext=(6, -4),
                fontsize=8.5,
                color=TEXT_SECONDARY,
                va="top",
            )

        ax.errorbar(x[i] - width / 2, m_msd, yerr=se_msd, fmt="o", color=OURS,
                    capsize=5, markersize=8,
                    label="$D$ from MSD (Einstein)" if i == 0 else None)
        ax.errorbar(x[i] + width / 2, m_vacf, yerr=se_vacf, fmt="s", color=REFERENCE,
                    capsize=5, markersize=8,
                    label="$D$ from VACF (Green-Kubo)" if i == 0 else None)

        ymax = max(m_msd + se_msd, m_vacf + se_vacf, d_exp or 0.0)
        all_top.append(ymax)
        ax.text(x[i], ymax + 0.003, "agree" if agree else "DISAGREE",
                ha="center", fontsize=9, fontweight="bold", color=(GOOD if agree else BAD))

    ax.set_xlim(-0.5, len(STATE_POINTS) - 0.5)
    ax.set_ylim(0, max(all_top) * 1.35)
    ax.set_xticks(x)
    ax.set_xticklabels([lbl for _, lbl in STATE_POINTS])
    ax.set_ylabel(r"self-diffusion coefficient $D^*$")
    ax.set_title(f"Experiment 04 — argon self-diffusion ({n_replicas} replicas/point)", fontsize=12)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()

    save(fig, OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
