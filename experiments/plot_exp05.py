#!/usr/bin/env python3
"""Experiment 05 plot: our simulated U* and p* vs. NIST's tabulated MD
isotherm (mmlapps.nist.gov/srs/LJ_PURE/md.htm), as a function of density.
Two stacked panels (one axis each -- never a dual-axis chart) rather than
one panel with two y-scales. Reads results/exp05_lj_eos.csv, writes
results/exp05_lj_eos.png. No simulation here -- see exp05_lj_eos.sh.
"""
import numpy as np

from plot_common import BAD, OURS, REFERENCE, ROOT, apply_style, load_csv, save

IN_CSV = ROOT / "results" / "exp05_lj_eos.csv"
OUT_PNG = ROOT / "results" / "exp05_lj_eos.png"


def main() -> None:
    apply_style()
    header, data = load_csv(IN_CSV)

    order = np.argsort(data["rho_star"])
    d = data[order]
    rho = d["rho_star"].astype(float)

    import matplotlib.pyplot as plt

    fig, (ax_u, ax_p) = plt.subplots(2, 1, figsize=(7.5, 8), sharex=True)

    ax_u.plot(rho, d["our_U_star"], "o-", color=OURS, label="our simulation")
    ax_u.plot(rho, d["nist_U_star"], "s--", color=REFERENCE, label="NIST (md.htm)")
    ax_u.set_ylabel(r"$U^*$ per atom")
    ax_u.set_title("Experiment 05 — LJ equation of state vs. NIST, $T^*\\approx0.85$")
    ax_u.legend(loc="lower left")

    ax_p.plot(rho, d["our_p_star"], "o-", color=OURS, label="our simulation")
    ax_p.plot(rho, d["nist_p_star"], "s--", color=REFERENCE, label="NIST (md.htm)")
    ax_p.set_ylabel(r"$p^*$")
    ax_p.set_xlabel(r"$\rho^*$")

    for i, row in enumerate(d):
        du = abs(float(row["our_U_star"]) - float(row["nist_U_star"]))
        dp = abs(float(row["our_p_star"]) - float(row["nist_p_star"]))
        y_off = 8 if i % 2 == 0 else 18
        ax_u.annotate(f"$\\Delta$={du:.3f}", (float(row["rho_star"]), float(row["our_U_star"])),
                      textcoords="offset points", xytext=(0, y_off), fontsize=7.5, ha="center",
                      color="#52514e")
        ax_p.annotate(f"$\\Delta$={dp:.3f}", (float(row["rho_star"]), float(row["our_p_star"])),
                      textcoords="offset points", xytext=(0, y_off), fontsize=7.5, ha="center",
                      color="#52514e")

    # Flag rho*=0.9: far outside tolerance and investigated (see
    # data/nist/README.md / task-15 report) -- likely residual order from
    # the fcc starting lattice that 50 t* of equilibration did not fully
    # randomize this close to the freezing line, not a units/convention bug
    # (the other five points, and the independent exp01 NIST energy/virial
    # check, all agree well).
    du_last = abs(float(d["our_U_star"][-1]) - float(d["nist_U_star"][-1]))
    if du_last > 0.1:
        ax_u.annotate("investigated outlier\n(see report)", (rho[-1], d["our_U_star"][-1]),
                      textcoords="offset points", xytext=(-70, -18), fontsize=8, color=BAD,
                      fontweight="bold", ha="center")

    fig.tight_layout()
    save(fig, OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
