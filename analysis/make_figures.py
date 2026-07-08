"""S5: publication data figures (vector PDF) for the manuscript.

Produces, into the manuscript's figures/ directory:
    morphometry_panel.pdf   along-tube morphometry (replaces placeholder)
    roof_panel.pdf          roof thickness + stability (replaces placeholder)
    planetary_panel.pdf     Mars/Moon comparison (replaces placeholder)

Style: one palette across all panels, colour-blind-safe, 8 pt text, no chart
junk. Every panel regenerates from analysis_out/ by re-running this script.

Run:  .venv/bin/python analysis/make_figures.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
FIGS = (ROOT.parent / "_Nature__Autonomous_aerial_reconnaissance_of_a_basaltic_"
        "lava_tube_reveals_interior_morphology_and_roof_thickness_distribution_"
        "for_planetary_subsurface" / "figures")

C = dict(  # shared palette (Okabe-Ito based)
    tube="#0072B2", surface="#E69F00", skylight="#D55E00",
    intact="#009E73", flagged="#999999", envelope="#000000",
    mars="#D55E00", moon="#56B4E9", earth="#009E73",
)
plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42,
})


def load_roof():
    rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    f = lambda k: np.array([float(r[k]) for r in rows])
    return rows, f("s"), f("x"), f("y"), f("tau"), f("span"), \
        np.array([r["class"] for r in rows])


def skylight_bands(s, klass):
    """Contiguous s-intervals of skylight stations."""
    bands, start = [], None
    for i in range(len(s)):
        if klass[i] == "skylight" and start is None:
            start = s[i]
        elif klass[i] != "skylight" and start is not None:
            bands.append((start, s[i - 1])); start = None
    if start is not None:
        bands.append((start, s[-1]))
    return bands


def fig_morphometry():
    rows = list(csv.DictReader(open(OUT / "morphometry.csv")))
    ok = [r for r in rows if r["area"]]
    s = np.array([float(r["s"]) for r in ok])
    A = np.array([float(r["area"]) for r in ok])
    W = np.array([float(r["width"]) for r in ok])
    H = np.array([float(r["height"]) for r in ok])
    eta = np.array([float(r["eta"]) for r in ok])
    _, s_all, x, y, tau, span, klass = load_roof()
    bands = skylight_bands(s_all, klass)

    fig = plt.figure(figsize=(7.1, 7.6))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.5, 1, 1, 1], hspace=0.5,
                          wspace=0.28)

    ax = fig.add_subplot(gs[0, :])   # (a) plan view with centreline
    ax.plot(x, y, color=C["tube"], lw=1.8)
    for i, (lo, hi) in enumerate(bands):
        m = (s_all >= lo) & (s_all <= hi)
        ax.plot(x[m], y[m], color=C["skylight"], lw=3.5,
                label="skylight" if i == 0 else None)
    for tick in range(0, int(s.max()) + 1, 50):
        j = np.argmin(np.abs(s_all - tick))
        ax.annotate(f"{tick}", (x[j], y[j]), fontsize=6, color="0.35",
                    xytext=(3, 3), textcoords="offset points")
    ax.set_aspect("equal")
    ax.set_xlabel("easting (m, local frame)")
    ax.set_ylabel("northing (m)")
    ax.set_title("a  Centreline of the surveyed conduit (1 m stations, "
                 "ticks every 50 m)", loc="left")
    ax.legend(frameon=False)

    panels = [("b  Cross-section area", A, r"$A$ (m$^2$)"),
              ("c  Width and height", None, "extent (m)"),
              ("d  Aspect ratio", eta, r"$\eta$")]
    for i, (title, val, ylab) in enumerate(panels):
        ax = fig.add_subplot(gs[1 + i, 0])
        for lo, hi in bands:
            ax.axvspan(lo, hi, color=C["skylight"], alpha=0.15, lw=0)
        if title.startswith("c"):
            ax.plot(s, W, lw=0.9, color=C["tube"], label="width")
            ax.plot(s, H, lw=0.9, color=C["intact"], label="height")
            ax.legend(frameon=False, ncols=2)
        else:
            ax.plot(s, val, lw=0.9, color=C["tube"])
        ax.set_ylabel(ylab)
        ax.set_title(title, loc="left")
        if i == 2:
            ax.set_xlabel("distance along tube $s$ (m)")

    ax = fig.add_subplot(gs[1, 1])   # (e) area histogram
    ax.hist(A, bins=24, color=C["tube"], alpha=0.85)
    med, q25, q75 = np.median(A), *np.percentile(A, [25, 75])
    ax.axvline(med, color="k", lw=1)
    ax.axvspan(q25, q75, color="k", alpha=0.08, lw=0)
    ax.set_xlabel(r"$A$ (m$^2$)")
    ax.set_ylabel("stations")
    ax.set_title("e  Area distribution (median, IQR)", loc="left")

    ax = fig.add_subplot(gs[2, 1])   # (f) sinuosity: recompute per window
    cl = list(csv.DictReader(open(OUT / "centreline.csv")))
    P = np.array([[float(r["x"]), float(r["y"]), float(r["z"])] for r in cl])
    Wn = 50
    sv = np.full(len(P), np.nan)
    for k in range(len(P)):
        lo, hi = max(0, k - Wn // 2), min(len(P) - 1, k + Wn // 2)
        if hi - lo >= Wn // 2:
            sv[k] = (hi - lo) / max(np.linalg.norm(P[hi] - P[lo]), 1e-6)
    ax.plot(np.arange(len(P)), sv, lw=0.9, color=C["tube"])
    for lo, hi in bands:
        ax.axvspan(lo, hi, color=C["skylight"], alpha=0.15, lw=0)
    ax.set_ylabel("sinuosity (50 m)")
    ax.set_xlabel("$s$ (m)")
    ax.set_title("f  Sinuosity", loc="left")

    ax = fig.add_subplot(gs[3, 1])   # (g) elevation profile
    ax.plot(np.arange(len(P)), P[:, 2], lw=1.1, color=C["tube"])
    for lo, hi in bands:
        ax.axvspan(lo, hi, color=C["skylight"], alpha=0.15, lw=0)
    ax.set_ylabel("centreline z (m)")
    ax.set_xlabel("$s$ (m)")
    ax.set_title("g  Elevation profile (65 m descent)", loc="left")

    fig.savefig(FIGS / "morphometry_panel.pdf", bbox_inches="tight")
    plt.close(fig)
    print("morphometry_panel.pdf")


def fig_roof():
    rows, s, x, y, tau, span, klass = load_roof()
    sig = float(rows[0]["sigma_tau"])
    it = klass == "intact"
    mp = klass == "multipass"
    bands = skylight_bands(s, klass)

    fig = plt.figure(figsize=(7.1, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.42,
                          wspace=0.3)

    ax = fig.add_subplot(gs[0, :])   # (a) tau profile
    for lo, hi in bands:
        ax.axvspan(lo, hi, color=C["skylight"], alpha=0.15, lw=0)
    ax.errorbar(s[it], tau[it], yerr=sig, fmt="o", ms=2.5, lw=0,
                elinewidth=0.5, color=C["intact"], label="intact roof")
    ax.plot(s[mp], np.zeros(mp.sum()) - 0.6, "|", ms=6, color=C["flagged"],
            label="excluded (multipass artefact)")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlim(-5, 180)
    ax.set_xlabel("distance along tube $s$ (m)")
    ax.set_ylabel(r"roof thickness $\tau$ (m)")
    ax.set_title("a  Roof thickness where the surface model overlaps the "
                 "survey (orange bands: skylights)", loc="left")
    ax.legend(frameon=False, loc="upper left")

    ax = fig.add_subplot(gs[1, 0])   # (b) histogram
    ax.hist(tau[it], bins=16, color=C["intact"], alpha=0.9)
    med = np.median(tau[it])
    ax.axvline(med, color="k", lw=1)
    ax.set_xlabel(r"$\tau$ (m)")
    ax.set_ylabel("stations")
    ax.set_title(f"b  Distribution (median {med:.1f} m)", loc="left")
    ax2 = ax.twiny()
    ax2.set_xlim(np.array(ax.get_xlim()) * 300.0 / 1000.0)
    ax2.set_xlabel(r"areal shielding mass (10$^3$ g cm$^{-2}$)", fontsize=7)

    ax = fig.add_subplot(gs[1, 1])   # (c) tau/L envelope
    r_it = tau[it] / span[it]
    kappa = float(json.loads((OUT / "roof_summary.json").read_text())["kappa_env"])
    ax.scatter(span[it], tau[it], s=10, color=C["intact"], label="intact roof")
    sky_mask = klass == "skylight"
    ax.scatter(span[sky_mask], np.zeros(sky_mask.sum()), s=14, marker="v",
               color=C["skylight"], label=r"skylights ($\tau\to0$)")
    Ls = np.linspace(4, 27, 50)
    ax.plot(Ls, kappa * Ls, "k--", lw=1,
            label=rf"$\tau/L={kappa:.3f}$ (observed min.)")
    ax.set_xlabel("local span $L$ (m)")
    ax.set_ylabel(r"$\tau$ (m)")
    ax.set_title("c  Thickness against span", loc="left")
    ax.legend(frameon=False, loc="upper left")

    fig.savefig(FIGS / "roof_panel.pdf", bbox_inches="tight")
    plt.close(fig)
    print("roof_panel.pdf")


def fig_planetary():
    cat = list(csv.DictReader(open(OUT / "planetary_catalogue.csv")))
    mars = np.array([float(r["aperture_long_axis_m"]) for r in cat
                     if r["body"] == "mars" and r["aperture_long_axis_m"]])
    moon = np.array([float(r["aperture_long_axis_m"]) for r in cat
                     if r["body"] == "moon" and r["aperture_long_axis_m"]])
    sky = json.loads((OUT / "aerial_holes.json").read_text())
    ours = [2 * h["semi_major"] for h in sky["holes"][:3]]
    pla = json.loads((OUT / "planetary_summary.json").read_text())

    fig = plt.figure(figsize=(7.1, 3.1))
    gs = fig.add_gridspec(1, 2, wspace=0.3)

    ax = fig.add_subplot(gs[0, 0])   # (a) aperture ECDFs
    for arr, key, name in ((mars, "mars", f"Mars APCs (n={len(mars)})"),
                           (moon, "moon", f"lunar pits (n={len(moon)})")):
        xs = np.sort(arr)
        ax.step(xs, np.arange(1, len(xs) + 1) / len(xs), where="post",
                color=C[key], label=name)
    for i, d in enumerate(sorted(ours)):
        ax.axvline(d, color=C["earth"], lw=1.2, ls=":",
                   label="Raufarhólshellir skylights" if i == 0 else None)
    ax.set_xscale("log")
    ax.set_xlabel("aperture long axis (m)")
    ax.set_ylabel("cumulative fraction")
    ax.set_title("a  Orbital aperture catalogues against the surveyed "
                 "skylights", loc="left")
    ax.legend(frameon=False, loc="upper left")

    ax = fig.add_subplot(gs[0, 1])   # (b) gravity-scaled stable span
    kappa = pla["kappa_env_earth"]
    taus = np.linspace(0, 14, 100)
    for body, g_scale, key in (("Earth", 1.0, "earth"),
                               ("Mars", pla["span_scale"]["mars"], "mars"),
                               ("Moon", pla["span_scale"]["moon"], "moon")):
        ax.plot(taus, g_scale * taus / kappa, color=C[key], lw=1.4,
                label=f"{body} (x{g_scale:.2f})" if body != "Earth"
                else "Earth (measured)")
    ax.fill_between(taus, 0, taus / kappa, color=C["earth"], alpha=0.07, lw=0)
    ax.set_xlabel(r"roof thickness $\tau$ (m)")
    ax.set_ylabel("max. span at the observed envelope (m)")
    ax.set_title("b  Stable span propagated by $L\\propto g^{-1/2}$",
                 loc="left")
    ax.legend(frameon=False, loc="upper left")

    fig.savefig(FIGS / "planetary_panel.pdf", bbox_inches="tight")
    plt.close(fig)
    print("planetary_panel.pdf")


if __name__ == "__main__":
    FIGS.mkdir(exist_ok=True)
    fig_morphometry()
    fig_roof()
    fig_planetary()
