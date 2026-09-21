#!/usr/bin/env python3
"""Experiment 05b plot: the diagnostic evidence behind experiment 05's
rho*=0.9 finding -- MSD growing over time (evidence the state is genuinely
diffusing, not caged) and g(r) decaying toward 1 (evidence it is
structurally liquid-like, not crystalline), with our D_msd compared
directly against NIST's own tabulated D*=0.027 at this exact state point.
Reads results/exp05b_rho090_diagnostic.csv, results/exp05b_msd_rho0.900.csv,
results/exp05b_rdf_rho0.900.csv, writes results/exp05b_rho090_diagnostic.png.
No simulation here -- see exp05b_rho090_diagnostic.sh.
"""
from plot_common import BAD, GOOD, OURS, REFERENCE, ROOT, apply_style, load_csv, save

SUMMARY_CSV = ROOT / "results" / "exp05b_rho090_diagnostic.csv"
MSD_CSV = ROOT / "results" / "exp05b_msd_rho0.900.csv"
RDF_CSV = ROOT / "results" / "exp05b_rdf_rho0.900.csv"
OUT_PNG = ROOT / "results" / "exp05b_rho090_diagnostic.png"


def main() -> None:
    apply_style()
    _, summary = load_csv(SUMMARY_CSV)
    row = summary[0]
    d_msd = float(row["D_msd"])
    nist_d = float(row["nist_D_star"])
    d_floor = 0.005

    _, msd = load_csv(MSD_CSV)
    _, rdf = load_csv(RDF_CSV)

    import matplotlib.pyplot as plt

    fig, (ax_msd, ax_rdf) = plt.subplots(1, 2, figsize=(12, 5))

    ax_msd.plot(msd["time"], msd["msd"], color=OURS, linewidth=1.8)
    ax_msd.set_xlabel(r"$t^*$")
    ax_msd.set_ylabel(r"MSD $\langle |\Delta r(t)|^2 \rangle$")
    ax_msd.set_title("MSD grows -- genuinely diffusing, not caged")
    ax_msd.text(
        0.03, 0.95,
        f"$D_{{msd}}$={d_msd:.4f}\nNIST $D^*$={nist_d:.4f}\nfloor={d_floor:.4f}",
        transform=ax_msd.transAxes, fontsize=9.5, va="top",
        color=(GOOD if d_msd > d_floor else BAD),
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7, edgecolor="none"),
    )

    ax_rdf.plot(rdf["r"], rdf["g"], color=OURS, linewidth=1.8)
    ax_rdf.axhline(1.0, color="#b8b7b2", linewidth=1, linestyle=":")
    ax_rdf.set_xlabel(r"$r^*$")
    ax_rdf.set_ylabel(r"$g(r^*)$")
    ax_rdf.set_title("g(r) decays to 1 -- liquid-like, not crystalline")

    fig.suptitle(
        r"Experiment 05b — diagnostic at $\rho^*$=0.900, $T^*$=0.851 "
        "(the exp05 outlier): not frozen, but suppressed mobility",
        fontsize=12.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    save(fig, OUT_PNG)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
