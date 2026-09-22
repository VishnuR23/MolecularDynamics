#!/usr/bin/env python3
"""Experiment 05 plot: our simulated U*, p* and achieved T* vs. NIST's
tabulated MD isotherm (mmlapps.nist.gov/srs/LJ_PURE/md.htm), as a function
of density. Three stacked panels (one axis each -- never a dual-axis
chart) rather than one panel with two y-scales. The T* panel exists
because this experiment's production now runs NVT (Nose-Hoover), pinned at
each row's own NIST-quoted T*, specifically so a reader can see the
comparison is on-target rather than confounded by a temperature offset
(see exp05_lj_eos.sh's header for why). Reads results/exp05_lj_eos.csv,
writes results/exp05_lj_eos.png. No simulation here -- see
exp05_lj_eos.sh.
"""
import numpy as np

from plot_common import BAD, GOOD, OURS, REFERENCE, ROOT, apply_style, load_csv, save

IN_CSV = ROOT / "results" / "exp05_lj_eos.csv"
OUT_PNG = ROOT / "results" / "exp05_lj_eos.png"


def main() -> None:
    apply_style()
    header, data = load_csv(IN_CSV)
    t_tol = float(header.get("tolerance", "").split("delta T*| <= ")[-1].split(",")[0].split()[0]) \
        if "delta T*" in header.get("tolerance", "") else 0.02

    order = np.argsort(data["rho_star"])
    d = data[order]
    rho = d["rho_star"].astype(float)

    import matplotlib.pyplot as plt

    fig, (ax_t, ax_u, ax_p) = plt.subplots(3, 1, figsize=(7.5, 10.5), sharex=True)

    ax_t.plot(rho, d["our_T_star"], "o-", color=OURS, label="our simulation (NVT, Nose-Hoover)")
    ax_t.plot(rho, d["nist_T_star"], "s--", color=REFERENCE, label="NIST target")
    ax_t.set_ylabel(r"$T^*$ (achieved)")
    ax_t.set_title("Experiment 05 — LJ EOS vs. NIST (NVT, pinned at NIST's own $T^*$ per row)",
                    fontsize=11.5)
    ax_t.legend(loc="lower left", fontsize=8.5)
    t_xy_offsets = [(-22, 14), (22, 22), (0, 12), (0, 24), (0, -16), (0, 14)]
    for i, row in enumerate(d):
        dt = abs(float(row["our_T_star"]) - float(row["nist_T_star"]))
        dx, dy = t_xy_offsets[i % len(t_xy_offsets)]
        ax_t.annotate(f"$\\Delta T^*$={dt:.4f}", (float(row["rho_star"]), float(row["our_T_star"])),
                      textcoords="offset points", xytext=(dx, dy), fontsize=7,
                      ha="center", color=(GOOD if dt <= t_tol else BAD))

    ax_u.plot(rho, d["our_U_star"], "o-", color=OURS, label="our simulation")
    ax_u.plot(rho, d["nist_U_star"], "s--", color=REFERENCE, label="NIST (md.htm)")
    ax_u.set_ylabel(r"$U^*$ per atom")
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

    # Flag rho*=0.9: far outside tolerance, and diagnosed. The
    # simulation never melts there -- it stays in the fcc lattice it was
    # started from, with D_msd ~ 1e-5 against NIST's tabulated 0.027 in
    # all three diagnostic seeds (experiment 05b). An ordered lattice is
    # both more cohesive (U* too negative) and less compressed at fixed
    # density (p* far too low), which is exactly the pair of signs seen
    # here. Not a units/convention bug -- the other five densities use
    # the identical code path and agree, as does the independent exp01
    # NIST energy/virial check -- and not a temperature-offset artifact,
    # since T* is pinned here too. Full write-up:
    # docs/findings/2026-09-21-rho090-outlier.md
    du_last = abs(float(d["our_U_star"][-1]) - float(d["nist_U_star"][-1]))
    if du_last > 0.1:
        ax_u.annotate("never melts — stays fcc\n(docs/findings/, exp05b)",
                      (rho[-1], d["our_U_star"][-1]),
                      textcoords="offset points", xytext=(-72, -20), fontsize=8, color=BAD,
                      fontweight="bold", ha="center")

    fig.tight_layout()
    save(fig, OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
