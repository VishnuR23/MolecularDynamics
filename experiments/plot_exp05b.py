#!/usr/bin/env python3
"""Experiment 05b plot: the diagnostic evidence behind experiment 05's
rho*=0.9 finding. Leads with three direct structural/dynamical signatures,
all read from the committed trajectory data:

  - VACF: a genuine backscattering minimum (~-0.30 at t*~0.2) -- the
    velocity decorrelates and partially reverses, as in a real liquid,
    not a rigid lattice.
  - MSD: a ballistic-to-caging-to-escape shape, not a clean straight
    diffusive line -- there is a real caging plateau (t*~0.5-2.5) before
    the atoms escape into the slow upward curve the D_msd fit uses.
  - g(r): does NOT cleanly decay to 1 within the accessible range -- after
    nearly reaching baseline near r*~2.4-2.8 it re-strengthens to
    g~1.3 near r*~3.05, an oscillation regrowing rather than damping out.
    Contrasted directly against experiment 03's Rahman g(r) (same code,
    lower density), which decays cleanly to ~1 by r*~4.5-5.

Together these say: genuinely mobile (not a rigid, zero-diffusion
crystal), but retaining real structural order (not a normal liquid
either) -- a partially-ordered, still-diffusing state near the freezing
line. The NIST D* comparison (our D_msd ~25% below NIST's D*=0.027) is
kept as supporting, not headline, evidence -- it is a single trajectory
(n=1), unlike experiment 04's replica-averaged, error-barred D.

Also states, directly on the figure, that neither diffusion estimate is
converged at this density: the running Green-Kubo integral is still
rising at the end of the accessible window (maxLag=2000 steps = t*=10,
capped by kMaxLagCap in apps/moldyn_run.cpp -- see the task-15 report for
why this isn't just raised), and the Einstein fit window straddles the
plateau-to-escape crossover, so its local slope is not yet linear.

Reads results/exp05b_rho090_diagnostic.csv, results/exp05b_msd_rho0.900.csv,
results/exp05b_rdf_rho0.900.csv, results/exp03_rdf_rahman.csv, writes
results/exp05b_rho090_diagnostic.png. No simulation here -- all figures
and numbers come from already-committed CSVs -- see
exp05b_rho090_diagnostic.sh.
"""
import numpy as np

from plot_common import BAD, GOOD, OURS, REFERENCE, ROOT, apply_style, load_csv, save

SUMMARY_CSV = ROOT / "results" / "exp05b_rho090_diagnostic.csv"
MSD_CSV = ROOT / "results" / "exp05b_msd_rho0.900.csv"
RDF_CSV = ROOT / "results" / "exp05b_rdf_rho0.900.csv"
RAHMAN_RDF_CSV = ROOT / "results" / "exp03_rdf_rahman.csv"
OUT_PNG = ROOT / "results" / "exp05b_rho090_diagnostic.png"


def running_gk_integral(times: np.ndarray, vacf: np.ndarray, d_vacf_full: float) -> np.ndarray:
    """Cumulative-trapezoid Green-Kubo D(t), calibrated so the value at the
    last committed lag reproduces moldyn_run's own reported D_vacf exactly
    (the CSV's vacf column is normalised to 1 at lag 0; moldyn_run
    integrates the unnormalised VACF, so this rescales rather than
    re-deriving the physical velocity scale)."""
    full = np.trapezoid(vacf, times)
    scale = 3.0 * d_vacf_full / full
    running = np.zeros_like(times)
    # cumulative trapezoid
    seg = 0.5 * (vacf[1:] + vacf[:-1]) * np.diff(times)
    running[1:] = np.cumsum(seg)
    return running * scale / 3.0


def main() -> None:
    apply_style()
    _, summary = load_csv(SUMMARY_CSV)
    row = summary[0]
    d_msd = float(row["D_msd"])
    nist_d = float(row["nist_D_star"])
    d_floor = 0.005

    _, msd = load_csv(MSD_CSV)
    _, rdf = load_csv(RDF_CSV)
    _, rdf_rahman = load_csv(RAHMAN_RDF_CSV)

    times = msd["time"].astype(float)
    vacf = msd["vacf"].astype(float)
    d_vacf_line = next(line for line in open(MSD_CSV) if line.startswith("# D_vacf"))
    d_vacf_full = float(d_vacf_line.split(":")[1].strip())
    running_d = running_gk_integral(times, vacf, d_vacf_full)
    d_at_5 = running_d[np.argmin(np.abs(times - 5.0))]
    d_at_10 = running_d[-1]

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    ax_vacf, ax_msd, ax_rdf = axes

    # --- Panel 1: VACF backscattering ---
    ax_vacf.plot(times, vacf, color=OURS, linewidth=1.6)
    ax_vacf.axhline(0.0, color="#b8b7b2", linewidth=1)
    imin = int(np.argmin(vacf))
    ax_vacf.plot([times[imin]], [vacf[imin]], marker="o", color=BAD, markersize=7, zorder=5)
    ax_vacf.annotate(f"backscattering\nmin={vacf[imin]:.2f} at $t^*$={times[imin]:.2f}",
                      (times[imin], vacf[imin]), textcoords="offset points", xytext=(35, 10),
                      fontsize=8.5, color=BAD)
    ax_vacf.set_xlim(0, 3)
    ax_vacf.set_xlabel(r"$t^*$")
    ax_vacf.set_ylabel(r"VACF (normalised)")
    ax_vacf.set_title("1. Velocity backscatters — genuinely dynamic", fontsize=11.5)
    ax_vacf.text(
        0.97, 0.95,
        f"running $D_{{GK}}$: {d_at_5:.4f} ($t^*$=5) $\\to$ {d_at_10:.4f} ($t^*$=10)\n"
        "still rising -- NOT converged\n"
        "(maxLag caps window at $t^*$=10)",
        transform=ax_vacf.transAxes, fontsize=8, va="top", ha="right", color=BAD,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, edgecolor="none"),
    )

    # --- Panel 2: MSD regimes ---
    ax_msd.plot(times, msd["msd"], color=OURS, linewidth=1.8, zorder=3)
    ax_msd.axvspan(0, 0.5, color="#eda100", alpha=0.12, zorder=0)
    ax_msd.axvspan(0.5, 2.5, color=BAD, alpha=0.10, zorder=0)
    ax_msd.axvspan(2.5, 7.5, color="#1baf7a", alpha=0.10, zorder=0)
    ax_msd.text(0.25, 0.03, "ballistic", transform=ax_msd.get_xaxis_transform(),
                fontsize=7.5, ha="center", color="#8a6d00")
    ax_msd.text(1.5, 0.03, "caging\nplateau", transform=ax_msd.get_xaxis_transform(),
                fontsize=7.5, ha="center", color=BAD)
    ax_msd.text(5.0, 0.03, "fit window\n(non-linear)", transform=ax_msd.get_xaxis_transform(),
                fontsize=7.5, ha="center", color="#0a7a4a")
    ax_msd.set_xlabel(r"$t^*$")
    ax_msd.set_ylabel(r"MSD $\langle |\Delta r(t)|^2 \rangle$")
    ax_msd.set_title("2. Caging plateau then escape — fit window not linear", fontsize=11.5)
    ax_msd.text(
        0.03, 0.95,
        f"local $D$ on [2.5,5): 0.0153\nlocal $D$ on [5,7.5): 0.0254\n"
        f"reported $D_{{msd}}$={d_msd:.4f} (likely an underestimate)",
        transform=ax_msd.transAxes, fontsize=8, va="top", color=BAD,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, edgecolor="none"),
    )

    # --- Panel 3: RDF re-strengthening ---
    ax_rdf.plot(rdf["r"], rdf["g"], color=OURS, linewidth=1.8,
                label=r"$\rho^*$=0.900 (this diagnostic)")
    ax_rdf.plot(rdf_rahman["r"], rdf_rahman["g"], color=REFERENCE, linewidth=1.3, alpha=0.8,
                linestyle="--", label=r"$\rho^*$=0.8177 (exp03 Rahman, same code)")
    ax_rdf.axhline(1.0, color="#b8b7b2", linewidth=1, linestyle=":")
    r_local = rdf["r"][(rdf["r"] > 2.8) & (rdf["r"] < 3.3)]
    g_local = rdf["g"][(rdf["r"] > 2.8) & (rdf["r"] < 3.3)]
    j = int(np.argmax(g_local))
    ax_rdf.annotate(f"re-strengthens\nto g={g_local[j]:.2f}", (r_local[j], g_local[j]),
                     textcoords="offset points", xytext=(15, 12), fontsize=8, color=BAD)
    ax_rdf.set_xlim(0, 5)
    ax_rdf.set_xlabel(r"$r^*$")
    ax_rdf.set_ylabel(r"$g(r^*)$")
    ax_rdf.set_title("3. g(r) does not decay — tail has not settled", fontsize=11.5)
    ax_rdf.legend(loc="upper right", fontsize=7.5)

    fig.suptitle(
        r"Experiment 05b — $\rho^*$=0.900, $T^*$=0.851 (the exp05 outlier): "
        "not frozen, but not a normal liquid either",
        fontsize=13,
    )
    fig.text(
        0.5, 0.015,
        f"Supporting (n=1, not error-barred like experiment 04): our $D_{{msd}}$={d_msd:.4f} vs. "
        f"NIST's tabulated $D^*$={nist_d:.4f} at this exact state point (floor for 'not frozen': "
        f"{d_floor:.3f}) — suppressed mobility, consistent with partial order.",
        ha="center", fontsize=8.5, color="#52514e",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))

    save(fig, OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
