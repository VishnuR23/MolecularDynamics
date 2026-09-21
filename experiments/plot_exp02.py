#!/usr/bin/env python3
"""Experiment 02 plot: log-log energy-error convergence of velocity-Verlet.
Reads results/exp02_energy_conservation.csv, writes
results/exp02_energy_conservation.png. No simulation here -- see
exp02_energy_conservation.sh.
"""
import numpy as np

from plot_common import GOOD, OURS, REFERENCE, ROOT, apply_style, load_csv, save

IN_CSV = ROOT / "results" / "exp02_energy_conservation.csv"
OUT_PNG = ROOT / "results" / "exp02_energy_conservation.png"


def main() -> None:
    apply_style()
    header, data = load_csv(IN_CSV)

    dt = data["dt"].astype(float)
    err = data["relative_energy_error"].astype(float)
    slope = float(header["fitted_slope"])

    order = np.argsort(dt)
    dt, err = dt[order], err[order]

    # dt^2 reference line, offset above the data (by a constant factor, so
    # slope 2 stays visible as a distinct guide rather than sitting exactly
    # on top of the near-perfect fit).
    ref = 1.6 * err[0] * (dt / dt[0]) ** 2

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.loglog(dt, err, "o-", color=OURS, label="our simulation (fitted slope %.3f)" % slope,
              zorder=3)
    ax.loglog(dt, ref, "--", color=REFERENCE, label=r"exact $dt^2$ reference (slope 2)", zorder=2)

    for x, y in zip(dt, err):
        ax.annotate(f"{x:g}", (x, y), textcoords="offset points", xytext=(6, -3), fontsize=8,
                    color="#52514e")

    ax.set_xlabel(r"timestep $dt^*$")
    ax.set_ylabel(r"relative energy error  $\max_t|E(t)-E(0)|\,/\,|E(0)|$")
    ax.set_title("Experiment 02 — velocity-Verlet energy conservation order")
    pass_ = 1.85 < slope < 2.15
    ax.text(
        0.02, 0.96,
        f"slope = {slope:.3f}  (tolerance [1.85, 2.15]: {'PASS' if pass_ else 'FAIL'})",
        transform=ax.transAxes, fontsize=10, color=(GOOD if pass_ else "#e34948"),
        fontweight="bold", va="top",
    )
    ax.legend(loc="lower right")
    fig.tight_layout()

    save(fig, OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
