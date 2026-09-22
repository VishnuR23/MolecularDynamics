#!/usr/bin/env python3
"""Experiment 05b plot: what state the system is actually in at rho*=0.900,
T*=0.851 -- the one density where experiment 05 disagrees with NIST.

Three panels, all read from the committed trajectory CSVs of the primary
seed, plus the per-seed verdicts from the summary CSV:

  1. MSD: flat. It rises ballistically and then stops. A liquid's MSD
     grows linearly forever; a bounded one is atoms vibrating about fixed
     lattice sites. The panel also reports the Lindemann ratio (per-axis
     rms displacement over nearest-neighbour distance), which lands on the
     textbook ~0.15 of a crystal at its melting point.
  2. VACF: backscatters at t*~0.2 and decays to zero. On its own that
     does not separate liquid from solid -- an atom rattling in a cage and
     an atom rattling in a lattice both reverse -- but with a flat MSD it
     says the reversals never add up to transport.
  3. g(r) against experiment 03's Rahman g(r) (identical code and
     normalisation, lower density, a genuine liquid). The principal peaks
     DO decay monotonically; an earlier version of this figure claimed an
     "oscillation regrowing rather than damping out", which the committed
     data never supported. The true statement is stronger: the envelope
     damps several times more slowly than the reference liquid's, and
     there is a split second-shell shoulder the reference does not have
     at all.

Every number drawn on the figure is computed here from the CSVs; none is
hardcoded. See docs/findings/2026-09-21-rho090-outlier.md.

Reads results/exp05b_rho090_diagnostic.csv,
results/exp05b_{msd,rdf}_rho0.900_seed*.csv and
results/exp03_rdf_rahman.csv; writes
results/exp05b_rho090_diagnostic.png. No simulation here -- see
exp05b_rho090_diagnostic.sh.
"""
import numpy as np

from plot_common import BAD, OURS, REFERENCE, ROOT, TEXT_SECONDARY, apply_style, load_csv, save

SUMMARY_CSV = ROOT / "results" / "exp05b_rho090_diagnostic.csv"
RAHMAN_RDF_CSV = ROOT / "results" / "exp03_rdf_rahman.csv"
OUT_PNG = ROOT / "results" / "exp05b_rho090_diagnostic.png"


def local_extrema(r, g):
    """Interior local maxima and minima of g(r), in increasing r."""
    maxima, minima = [], []
    for i in range(1, len(g) - 1):
        if g[i] >= g[i - 1] and g[i] > g[i + 1]:
            maxima.append((r[i], g[i]))
        if g[i] <= g[i - 1] and g[i] < g[i + 1]:
            minima.append((r[i], g[i]))
    return maxima, minima


def principal_peaks(maxima):
    """Coordination-shell peaks: the local maxima that rise above the bulk
    value g = 1. Sub-unity local maxima are shoulders inside a trough, not
    shells, and counting them as peaks is what made an earlier version of
    this figure report a decaying envelope as a growing one."""
    return [(r, g) for r, g in maxima if g > 1.0]


def main() -> None:
    header, summary = load_csv(SUMMARY_CSV)
    apply_style()

    # Provenance header values are free text ("5006 (the trajectory
    # plotted in the figure)"), so take the leading token.
    primary = int(header.get("primary_seed", str(summary["seed"][0])).split()[0])
    nist_d = float(summary["nist_D_star"][0])
    d_floor = float(header.get("d_floor", "0.005").split()[0])
    row = summary[summary["seed"] == primary][0]
    d_msd = float(row["D_msd"])

    _, msd_csv = load_csv(ROOT / "results" / f"exp05b_msd_rho0.900_seed{primary}.csv")
    _, rdf = load_csv(ROOT / "results" / f"exp05b_rdf_rho0.900_seed{primary}.csv")
    _, rdf_ref = load_csv(RAHMAN_RDF_CSV)

    times = msd_csv["time"].astype(float)
    msd = msd_csv["msd"].astype(float)
    vacf = msd_csv["vacf"].astype(float)

    import matplotlib.pyplot as plt

    fig, (ax_msd, ax_vacf, ax_rdf) = plt.subplots(1, 3, figsize=(16, 5.2))

    # --- Panel 1: the MSD does not grow ------------------------------
    ax_msd.plot(times, msd, color=OURS, linewidth=1.8, zorder=3)

    plateau = msd[times >= 1.0]
    lo, hi = plateau.min(), plateau.max()
    ax_msd.axhspan(lo, hi, color=BAD, alpha=0.12, zorder=0)

    # For scale: what a genuinely diffusive liquid at NIST's own D* would
    # have done over the same window, 6*D*t.
    ax_msd.plot(times, 6.0 * nist_d * times, color=REFERENCE, linestyle="--", linewidth=1.5,
                label=rf"a liquid at NIST's $D^*$={nist_d:g}  ($6D^*t$)")
    ax_msd.plot(times, msd, color=OURS, linewidth=1.8, zorder=3,
                label=r"measured, $\rho^*$=0.900")

    rms_axis = np.sqrt(msd[-1] / 3.0)
    # Nearest-neighbour distance: the position of the first coordination
    # shell, i.e. the first principal peak of this run's own g(r).
    _first_maxima, _ = local_extrema(rdf["r"].astype(float), rdf["g"].astype(float))
    nn = principal_peaks(_first_maxima)[0][0]
    lindemann = rms_axis / nn

    ax_msd.set_ylim(0, max(hi * 3.0, 6.0 * nist_d * times[-1] * 0.12))
    ax_msd.set_xlabel(r"$t^*$")
    ax_msd.set_ylabel(r"MSD $\langle |\Delta r(t)|^2 \rangle$")
    ax_msd.set_title("1. The MSD is flat — atoms vibrate, they do not travel", fontsize=11.5)
    ax_msd.text(
        0.04, 0.95,
        f"bounded at {lo:.4f}–{hi:.4f} for $t^*$≥1\n"
        f"per-axis rms {rms_axis:.3f}$\\sigma$ / nn {nn:.3f}$\\sigma$\n"
        f"Lindemann ratio {lindemann:.3f} (a melting crystal is ≈0.15)\n"
        f"$D_{{msd}}$={d_msd:.2e}, floor for 'not frozen' {d_floor:g}",
        transform=ax_msd.transAxes, fontsize=8, va="top", color=BAD,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="none"),
    )
    ax_msd.legend(loc="lower right", fontsize=8)

    # --- Panel 2: VACF ------------------------------------------------
    ax_vacf.plot(times, vacf, color=OURS, linewidth=1.6)
    ax_vacf.axhline(0.0, color="#b8b7b2", linewidth=1)
    imin = int(np.argmin(vacf))
    ax_vacf.plot([times[imin]], [vacf[imin]], marker="o", color=BAD, markersize=7, zorder=5)
    ax_vacf.annotate(
        f"backscattering min {vacf[imin]:.2f}\nat $t^*$={times[imin]:.2f}",
        (times[imin], vacf[imin]), textcoords="offset points", xytext=(30, 6),
        fontsize=8.5, color=BAD,
    )
    ax_vacf.set_xlim(0, 3)
    ax_vacf.set_xlabel(r"$t^*$")
    ax_vacf.set_ylabel("VACF (normalised)")
    ax_vacf.set_title("2. Velocity reverses, but the reversals go nowhere", fontsize=11.5)
    ax_vacf.text(
        0.96, 0.95,
        "an atom in a cage and an atom on a\nlattice site both backscatter —\n"
        "panel 1 is what separates them",
        transform=ax_vacf.transAxes, fontsize=8, va="top", ha="right", color=TEXT_SECONDARY,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="none"),
    )

    # --- Panel 3: g(r) envelope vs. the reference liquid --------------
    r_ours, g_ours = rdf["r"].astype(float), rdf["g"].astype(float)
    r_ref, g_ref = rdf_ref["r"].astype(float), rdf_ref["g"].astype(float)

    ax_rdf.plot(r_ours, g_ours, color=OURS, linewidth=1.8,
                label=r"$\rho^*$=0.900 (this diagnostic)")
    ax_rdf.plot(r_ref, g_ref, color=REFERENCE, linewidth=1.3, alpha=0.85, linestyle="--",
                label=r"$\rho^*$=0.8177 (exp03 Rahman liquid, same code)")
    ax_rdf.axhline(1.0, color="#b8b7b2", linewidth=1, linestyle=":")

    max_ours, min_ours = local_extrema(r_ours, g_ours)
    max_ref, _ = local_extrema(r_ref, g_ref)
    peaks_ours, peaks_ref = principal_peaks(max_ours), principal_peaks(max_ref)

    third_r, third_g = peaks_ours[2]
    ref_r, ref_g = peaks_ref[2]
    env_ours, env_ref = abs(third_g - 1.0), abs(ref_g - 1.0)
    ratio = env_ours / env_ref

    ax_rdf.plot([third_r], [third_g], marker="o", color=BAD, markersize=6, zorder=5)
    ax_rdf.annotate(f"3rd peak $|g-1|$={env_ours:.2f}", (third_r, third_g),
                     textcoords="offset points", xytext=(10, 12), fontsize=8, color=BAD)
    ax_rdf.plot([ref_r], [ref_g], marker="s", color=REFERENCE, markersize=5, zorder=5)
    ax_rdf.annotate(f"reference $|g-1|$={env_ref:.3f}", (ref_r, ref_g),
                     textcoords="offset points", xytext=(6, -24), fontsize=8, color=REFERENCE)

    # The split second-shell shoulder: two minima with a sub-unity local
    # maximum between them, where the reference has one smooth trough.
    shoulder = [(r, g) for r, g in max_ours if 2.2 < r < 3.0 and g < 1.0]
    troughs = [(r, g) for r, g in min_ours if 2.2 < r < 3.0]
    if shoulder and len(troughs) >= 2:
        ax_rdf.axvspan(troughs[0][0], troughs[-1][0], color=BAD, alpha=0.10, zorder=0)
        ax_rdf.text(
            0.5 * (troughs[0][0] + troughs[-1][0]), 0.05,
            f"split shoulder\n(minima {troughs[0][0]:.2f}, {troughs[-1][0]:.2f};\n"
            f"no counterpart in the reference)",
            transform=ax_rdf.get_xaxis_transform(), fontsize=7.5, ha="center", color=BAD,
        )

    peak_text = " → ".join(f"{g:.3f}" for _, g in peaks_ours[:3])
    ax_rdf.set_xlim(0, 5)
    ax_rdf.set_xlabel(r"$r^*$")
    ax_rdf.set_ylabel(r"$g(r^*)$")
    ax_rdf.set_title(f"3. Peaks decay ({peak_text}), but the\n"
                      f"envelope sits {ratio:.1f}x further from 1 than the liquid's",
                      fontsize=10.5)
    ax_rdf.legend(loc="upper right", fontsize=7.5)

    n_pass = int((summary["status"] == "PASS").sum())
    n_seeds = len(summary)
    seed_text = ",  ".join(
        f"seed {int(s)}: $D_{{msd}}$={float(d):.1e}" for s, d in zip(summary["seed"],
                                                                     summary["D_msd"])
    )

    fig.suptitle(
        r"Experiment 05b — $\rho^*$=0.900, $T^*$=0.851 (experiment 05's one failing row): "
        "the simulation never melts",
        fontsize=13,
    )
    fig.text(
        0.5, 0.055,
        f"{n_pass} of {n_seeds} independent seeds cleared the a-priori 'not frozen' floor "
        f"$D_{{msd}}>{d_floor:g}$.   {seed_text}.",
        ha="center", fontsize=8.5, color=TEXT_SECONDARY,
    )
    fig.text(
        0.5, 0.022,
        f"NIST tabulates $D^*$={nist_d:g} for this same state point — about "
        f"{nist_d / float(summary['D_msd'].astype(float).mean()):.0f}x larger than the mean "
        "measured here.   See docs/findings/2026-09-21-rho090-outlier.md.",
        ha="center", fontsize=8.5, color=TEXT_SECONDARY,
    )
    fig.tight_layout(rect=(0, 0.085, 1, 0.93))

    save(fig, OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
