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

import os

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out")))
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
    # The conduit is strongly elongated (~4:1), so an equal-aspect plan in the
    # E/N frame collapses to a thin vertical strip. Rotate to the centreline's
    # principal axis so the plan view fills the wide panel.
    XY = np.column_stack([x, y])
    ctr = XY.mean(0)
    _, _, Vt = np.linalg.svd(XY - ctr, full_matrices=False)
    uv = (XY - ctr) @ Vt.T
    u, v = uv[:, 0], uv[:, 1]
    if u[-1] < u[0]:                  # orient travel direction to +u
        u = -u
    ax.plot(u, v, color=C["tube"], lw=1.8)
    for i, (lo, hi) in enumerate(bands):
        m = (s_all >= lo) & (s_all <= hi)
        ax.plot(u[m], v[m], color=C["skylight"], lw=3.5,
                label="skylight" if i == 0 else None)
    for tick in range(0, int(s.max()) + 1, 50):
        j = np.argmin(np.abs(s_all - tick))
        ax.annotate(f"{tick} m", (u[j], v[j]), fontsize=6, color="0.35",
                    xytext=(0, 6), textcoords="offset points", ha="center")
    ax.set_aspect("equal")
    ax.margins(x=0.02, y=0.25)
    ax.set_xlabel("along principal axis (m, local frame)")
    ax.set_ylabel("across (m)")
    ax.set_title("a  Plan view of the surveyed centreline (stations every 50 m)",
                 loc="left")
    ax.legend(frameon=False, loc="lower right")

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
            ax.legend(frameon=False, ncols=2, loc="upper left")
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
    desc = np.percentile(P[:, 2], 99) - np.percentile(P[:, 2], 1)
    ax.set_title(f"g  Elevation profile ({desc:.0f} m relief, no net descent)",
                 loc="left")

    fig.savefig(FIGS / "morphometry_panel.pdf", bbox_inches="tight")
    plt.close(fig)
    print("morphometry_panel.pdf")


def fig_roof():
    rows, s, x, y, tau, span, klass = load_roof()
    sig = float(rows[0]["sigma_tau"])
    it = klass == "intact"
    mp = np.isin(klass, ["multipass", "inconsistent", "low_coverage", "no_dem"])
    bands = skylight_bands(s, klass)

    fig = plt.figure(figsize=(7.1, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.42,
                          wspace=0.3)

    ax = fig.add_subplot(gs[0, :])   # (a) tau profile
    for lo, hi in bands:
        ax.axvspan(lo, hi, color=C["skylight"], alpha=0.15, lw=0)
    # sigma_tau is dominated by the common registration datum, a systematic
    # band shared by every station, not independent per-station noise.
    o = np.argsort(s[it])
    ax.fill_between(s[it][o], tau[it][o] - sig, tau[it][o] + sig,
                    color=C["intact"], alpha=0.18, lw=0,
                    label=rf"systematic datum band ($\pm{sig:.1f}$ m)")
    ax.plot(s[it][o], tau[it][o], "o", ms=2.2, color=C["intact"],
            label="intact roof")
    if mp.any():
        ax.plot(s[mp], np.zeros(mp.sum()) - 0.6, "|", ms=6, color=C["flagged"],
                label="excluded by consistency screen")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlim(-5, s.max() + 5)
    ax.set_xlabel("distance along tube $s$ (m)")
    ax.set_ylabel(r"roof thickness $\tau$ (m)")
    ax.set_title("a  Roof thickness along the surveyed length "
                 "(orange bands: skylights)", loc="left")
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

    ax = fig.add_subplot(gs[1, 1])   # (c) tau-L with plate-model iso-strength
    ax.scatter(span[it], tau[it], s=10, color=C["intact"], label="intact roof")
    sky_mask = klass == "skylight"
    ax.scatter(span[sky_mask], np.zeros(sky_mask.sum()), s=14, marker="v",
               color=C["skylight"], label=r"skylights (failed, $\tau\to0$)")
    # minimum thickness for stability under the clamped-strip plate model,
    # tau_min = beta rho g L^2 / sigma_t, at rock-mass tensile strengths
    beta, rho = 0.5, 2600.0
    Ls = np.linspace(4, 28, 80)
    for st_mpa, ls in ((1.0, "--"), (5.0, "-."), (10.0, ":")):
        ax.plot(Ls, beta * rho * 9.81 * Ls ** 2 / (st_mpa * 1e6), "k",
                ls=ls, lw=0.9,
                label=rf"$\sigma_t={st_mpa:.0f}$ MPa")
    ax.set_ylim(-0.8, 16.5)
    ax.set_xlabel("local span $L$ (m)")
    ax.set_ylabel(r"$\tau$ (m)")
    ax.set_title("c  Thickness against span, with the plate-model\n"
                 "stability limit", loc="left")
    ax.legend(frameon=False, loc="upper left", fontsize=6.2)

    fig.savefig(FIGS / "roof_panel.pdf", bbox_inches="tight")
    plt.close(fig)
    print("roof_panel.pdf")


def fig_planetary():
    """a: aperture ECDFs split by feature type, inner apertures for lunar
    pits (the opening into the void, not the funnel rim). b: plate-model
    stable span, L_max = sqrt(sigma_t tau / (beta rho g)), as a sensitivity
    band over rock-mass tensile strength, scaled to Mars and the Moon."""
    cat = list(csv.DictReader(open(OUT / "planetary_catalogue.csv")))

    def vals(rows, key):
        out = []
        for r in rows:
            v = (r.get(key) or "").strip()
            if v:
                out.append(float(v))
        return np.array(out)

    mars_apc = vals([r for r in cat if r["body"] == "mars"
                     and r["type"] == "APC"], "aperture_long_axis_m")
    moon_melt = vals([r for r in cat if r["body"] == "moon"
                      and r["type"] == "pit (impact melt)"], "inner_long_axis_m")
    moon_mare = vals([r for r in cat if r["body"] == "moon"
                      and r["type"] in ("pit (mare)", "pit (highland)")],
                     "inner_long_axis_m")
    sky = json.loads((OUT / "aerial_holes.json").read_text())
    ours = [2 * h["semi_major"] for h in sky["holes"][:3]]

    rev = json.loads((OUT / "review_stats.json").read_text())
    env = rev["envelope"]
    beta, rho = env["beta"], env["rho_mid"]

    fig = plt.figure(figsize=(7.1, 3.2))
    gs = fig.add_gridspec(1, 2, wspace=0.34)

    ax = fig.add_subplot(gs[0, 0])   # (a) aperture ECDFs by type
    for arr, col, ls, name in (
            (mars_apc, C["mars"], "-",
             f"Mars APCs, outer (n={len(mars_apc)})"),
            (moon_melt, C["moon"], "--",
             f"lunar impact-melt pits, inner (n={len(moon_melt)})"),
            (moon_mare, C["moon"], "-",
             f"lunar mare+highland pits, inner (n={len(moon_mare)})")):
        xs = np.sort(arr)
        ax.step(xs, np.arange(1, len(xs) + 1) / len(xs), where="post",
                color=col, ls=ls, label=name)
    for i, d in enumerate(sorted(ours)):
        ax.axvline(d, color=C["earth"], lw=1.2, ls=":",
                   label="Raufarhólshellir skylights" if i == 0 else None)
    ax.set_xscale("log")
    ax.set_xlabel("aperture long axis (m)")
    ax.set_ylabel("cumulative fraction")
    ax.set_title("a  Catalogued apertures by feature type", loc="left")
    ax.legend(frameon=False, loc="upper left", fontsize=6.2)

    ax = fig.add_subplot(gs[0, 1])   # (b) plate-model stable span band
    taus = np.linspace(0.5, 15, 120)
    g_body = dict(Earth=9.81, Mars=3.71, Moon=1.62)
    for (body, g), key in zip(g_body.items(), ("earth", "mars", "moon")):
        lo = np.sqrt(1e6 * taus / (beta * rho * g))     # sigma_t = 1 MPa
        hi = np.sqrt(10e6 * taus / (beta * rho * g))    # sigma_t = 10 MPa
        ax.fill_between(taus, lo, hi, color=C[key], alpha=0.18, lw=0)
        ax.plot(taus, np.sqrt(5e6 * taus / (beta * rho * g)), color=C[key],
                lw=1.3, label=f"{body}")
    it_tau = rev["stats"]["median"]
    ax.plot([it_tau], [np.sqrt(5e6 * it_tau / (beta * rho * 9.81))], "o",
            ms=4, color="k")
    ax.annotate("median roof,\nthis survey", (it_tau, 45), fontsize=6.5,
                ha="left", xytext=(it_tau + 0.6, 20))
    ax.set_xlabel(r"roof thickness $\tau$ (m)")
    ax.set_ylabel(r"$L_{\max}=\sqrt{\sigma_t\,\tau/(\beta\rho g)}$ (m)")
    ax.set_title(r"b  Plate-model stable span "
                 r"($\sigma_t$ = 1--10 MPa)", loc="left")
    ax.legend(frameon=False, loc="upper left")

    fig.savefig(FIGS / "planetary_panel.pdf", bbox_inches="tight")
    plt.close(fig)
    print("planetary_panel.pdf")


if __name__ == "__main__":
    FIGS.mkdir(exist_ok=True)
    fig_morphometry()
    fig_roof()
    fig_planetary()
