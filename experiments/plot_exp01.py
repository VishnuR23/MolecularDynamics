#!/usr/bin/env python3
"""Experiment 01 plot: relative deviation of our LJ energy/virial/tail
correction from the 12 NIST Standard Reference values. Reads
results/exp01_nist_table.csv, writes results/exp01_nist_table.png.
No simulation here -- see exp01_nist_table.sh.
"""
import pathlib

import numpy as np

from plot_common import BAD, GOOD, ROOT, apply_style, load_csv, save

OUT_CSV = ROOT / "results" / "exp01_nist_table.csv"
OUT_PNG = ROOT / "results" / "exp01_nist_table.png"


def main() -> None:
    apply_style()
    header, data = load_csv(OUT_CSV)

    n = len(data)
    labels = [f"{int(c)}/{t}/{rc:g}" for c, t, rc in zip(data["cfg"], data["trunc"], data["rc"])]
    u_dev = data["u_dev"].astype(float)
    w_dev = data["w_dev"].astype(float)
    ulrc_dev = np.array(
        [d if np.isfinite(d) and d > 0 else 1e-9 for d in data["ulrc_dev"].astype(float)]
    )

    tol = float(header.get("tolerance", "5e-5").split()[0])

    x = np.arange(n)
    width = 0.26

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.set_yscale("log")
    ax.bar(x - width, u_dev, width, label="U (pair energy)", color="#2a78d6")
    ax.bar(x, w_dev, width, label="W (virial)", color="#1baf7a")
    ax.bar(x + width, ulrc_dev, width, label="U_LRC (tail correction)", color="#eda100")
    ax.axhline(tol, color=BAD, linestyle="--", linewidth=1.5, label=f"tolerance ({tol:.0e})")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("NIST case: config / truncation / $r_c^*$")
    ax.set_ylabel("relative deviation from NIST reference")
    ax.set_title("Experiment 01 — LJ energy, virial and tail correction vs. the 12 NIST cases")
    all_pass = all(s == "PASS" for s in data["status"])
    subtitle_color = GOOD if all_pass else BAD
    ax.text(
        0.01,
        1.06,
        "ALL 12 CASES PASS" if all_pass else "ONE OR MORE CASES FAIL",
        transform=ax.transAxes,
        fontsize=10,
        color=subtitle_color,
        fontweight="bold",
    )
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()

    save(fig, OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
