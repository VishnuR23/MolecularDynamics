#!/usr/bin/env python3
"""Oxygen-oxygen radial distribution function of SPC/E water at ambient."""

import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import plot_common as pc  # noqa: E402

CSV = pc.ROOT / "results" / "water_rdf.csv"
PNG = pc.ROOT / "results" / "water_rdf.png"


def main() -> int:
    pc.apply_style()
    header, data = pc.load_csv(CSV)
    r = data["r_angstrom"]
    g = data["g_oo"]

    peak_r = float(header["first_peak_r_angstrom"])
    peak_g = float(header["first_peak_g"])
    ref_r = float(header["reference_first_peak_r_angstrom"])
    meanT = float(header["mean_production_temperature_K"])
    targetT = float(header["target_temperature_K"])
    tail = float(header["tail_mean_g_beyond_8.5A"])
    molecules = int(header["molecules"])

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(r, g, color=pc.OURS, lw=2.0, label="this engine (SPC/E, Ewald, rigid)")
    ax.axhline(1.0, color="0.6", lw=1.0, ls=":", zorder=1)
    ax.axvline(ref_r, color=pc.REFERENCE, lw=1.6, ls="--",
               label=f"published first peak, {ref_r:.2f} $\\AA$")
    ax.plot([peak_r], [peak_g], "o", color=pc.GOOD, ms=9, zorder=5)

    ax.annotate(
        f"first peak $r$ = {peak_r:.4f} $\\AA$, $g$ = {peak_g:.3f}\n"
        f"published $\\approx$ {ref_r:.2f} $\\AA$, $g \\approx$ 3.0  "
        f"({abs(peak_r - ref_r) / ref_r * 100:.2f}% in position)",
        xy=(peak_r, peak_g), xytext=(peak_r + 1.4, peak_g - 0.25),
        color=pc.OURS, fontsize=10,
        arrowprops=dict(arrowstyle="->", color=pc.OURS, lw=1.2),
    )
    ax.annotate(
        f"$g(r) \\to$ {tail:.4f} beyond 8.5 $\\AA$\n(normalisation check)",
        xy=(9.2, 1.0), xytext=(7.0, 1.55), fontsize=9, color="0.35",
        arrowprops=dict(arrowstyle="->", color="0.55", lw=1.0),
    )

    ax.set_xlabel(r"$r$  ($\AA$)")
    ax.set_ylabel(r"$g_{\mathrm{OO}}(r)$")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, max(3.4, peak_g * 1.12))
    ax.legend(loc="upper right", frameon=False)
    ax.set_title("SPC/E water — oxygen-oxygen structure at ambient conditions")
    # Sits under the axes, so the figure needs room reserved for it -- a
    # negative y here renders outside the canvas and gets clipped.
    fig.subplots_adjust(bottom=0.20)
    fig.text(
        0.5, 0.02,
        f"{molecules} rigid molecules, 0.997 g/cm$^3$, 30 ps equilibration + 20 ps production; "
        f"mean production temperature {meanT:.1f} K against a {targetT:.0f} K target.",
        ha="center", fontsize=9, color="0.35",
    )
    pc.save(fig, PNG)
    print(f"plot_water_rdf: wrote {PNG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
