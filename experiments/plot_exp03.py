#!/usr/bin/env python3
"""Experiment 03 plot: g(r) for both argon state points, with the
literature first-peak band marked. Reads results/exp03_rdf_{rahman,
verlet}.csv and results/exp03_argon_rdf.csv, writes
results/exp03_argon_rdf.png. No simulation here -- see exp03_argon_rdf.sh.
"""
import matplotlib.pyplot as plt

from plot_common import GOOD, OURS, REFERENCE, ROOT, apply_style, load_csv, save

SUMMARY_CSV = ROOT / "results" / "exp03_argon_rdf.csv"
OUT_PNG = ROOT / "results" / "exp03_argon_rdf.png"

STATE_POINTS = [
    ("rahman", "Rahman (1964): $\\rho^*$=0.8177, $T^*$=0.7880", OURS),
    ("verlet", "Verlet (1967): $\\rho^*$=0.8442, $T^*$=0.7280", REFERENCE),
]


def main() -> None:
    apply_style()
    header, summary = load_csv(SUMMARY_CSV)
    pos_lo, pos_hi = (float(x) for x in header["tolerance"].split("[")[1].split("]")[0].split(","))
    height_lo, height_hi = (
        float(x) for x in header["tolerance"].split("[")[2].split("]")[0].split(",")
    )

    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.axvspan(pos_lo, pos_hi, color="#d8d7d2", alpha=0.5, zorder=0,
               label=f"literature 1st-peak band ($r^*$: [{pos_lo:g},{pos_hi:g}], "
                     f"g: [{height_lo:g},{height_hi:g}])")
    ax.axhspan(height_lo, height_hi, color="#d8d7d2", alpha=0.5, zorder=0)

    for label, legend_label, color in STATE_POINTS:
        rdf_csv = ROOT / "results" / f"exp03_rdf_{label}.csv"
        _, rdf = load_csv(rdf_csv)
        ax.plot(rdf["r"], rdf["g"], color=color, label=legend_label, linewidth=1.8)

        row = summary[summary["label"] == label][0]
        peak_r, peak_g = float(row["peak_r_star"]), float(row["peak_g"])
        in_band = pos_lo <= peak_r <= pos_hi and height_lo <= peak_g <= height_hi
        ax.plot([peak_r], [peak_g], marker="o", color=(GOOD if in_band else "#e34948"),
                markersize=8, zorder=5, markeredgecolor="white", markeredgewidth=1)
        ax.annotate(f"peak: $r^*$={peak_r:.3f}, g={peak_g:.2f}", (peak_r, peak_g),
                    textcoords="offset points", xytext=(10, 8), fontsize=9, color=color)

    ax.axhline(1.0, color="#b8b7b2", linewidth=1, linestyle=":", zorder=0)
    ax.set_xlim(0, 3.5)
    ax.set_xlabel(r"$r^*$")
    ax.set_ylabel(r"$g(r^*)$")
    ax.set_title("Experiment 03 — argon radial distribution function")
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()

    save(fig, OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
