"""Extended Data figure panels (vector PDF) for the manuscript.

Produces, into the manuscript's figures/ directory:
    ed_registration_panel.pdf  ED Fig: skylight constellation match, vertical
                               residuals, and the DLIO-vs-FAST-LIO registration
                               sensitivity that motivated the gravity-preserving
                               4-DOF solution.
    ed_consistency_panel.pdf   ED Fig: map consistency screen (longitudinal
                               profile against the DEM, above-surface return
                               fraction, station classification footprint).
    ed_centreline_panel.pdf    ED Fig: L1-medial centreline over the shell
                               (plan and side views) and a gallery of
                               representative cross-sections.

Same style contract as make_figures.py: one Okabe-Ito palette, 8 pt text,
no chart junk, every panel regenerates from analysis_out/ + maps/.

Run:  env -u PYTHONPATH .venv/bin/python analysis/make_ed_figures.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

from lava_pcd.io.pcd_reader import BinaryPcdReader

import os

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out")))
OUT_BASE = ROOT / "analysis_out"
MAPS = ROOT / "maps"
# Tube working copies for the shell/section panels (Ed-frame by default when
# ANALYSIS_OUT points at the Ed run).
SHELL_PCD = Path(os.environ.get("EDFIG_SHELL_PCD", str(MAPS / "flf_30cm_aerial.pcd")))
SECT_PCD = Path(os.environ.get("EDFIG_SECT_PCD", str(MAPS / "flf_10cm_aerial.pcd")))
FIGS = (ROOT.parent / "_Nature__Autonomous_aerial_reconnaissance_of_a_basaltic_"
        "lava_tube_reveals_interior_morphology_and_roof_thickness_distribution_"
        "for_planetary_subsurface" / "figures")

C = dict(  # shared palette (Okabe-Ito based), identical to make_figures.py
    tube="#0072B2", surface="#E69F00", skylight="#D55E00",
    intact="#009E73", flagged="#999999", envelope="#000000",
    dlio="#CC79A7", shell="#BBBBBB",
)
CLASS_COLOUR = {"intact": C["intact"], "skylight": C["skylight"],
                "multipass": C["flagged"], "inconsistent": "#CC79A7"}
plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42,
})


def load_xyz(path: Path, keep: int | None = None) -> np.ndarray:
    parts = []
    with BinaryPcdReader(path) as r:
        for chunk in r.chunks():
            parts.append(chunk[:, :3].astype(np.float64))
    X = np.vstack(parts)
    X = X[np.isfinite(X).all(axis=1)]
    if keep is not None and len(X) > keep:
        rng = np.random.default_rng(0)
        X = X[rng.choice(len(X), keep, replace=False)]
    return X


def load_roof():
    rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    f = lambda k: np.array([float(r[k]) for r in rows])
    return (f("s"), f("x"), f("y"), f("z_dem"), f("z_ceil"), f("tau"),
            f("ghost_frac"), np.array([r["class"] for r in rows]))


def load_centreline():
    rows = list(csv.DictReader(open(OUT / "centreline.csv")))
    f = lambda k: np.array([float(r[k]) for r in rows])
    return f("s"), np.column_stack([f("x"), f("y"), f("z")]), \
        np.column_stack([f("tx"), f("ty"), f("tz")])


def principal_frame(xy: np.ndarray):
    """Rotation of the horizontal frame onto the centreline principal axis."""
    ctr = xy.mean(0)
    _, _, Vt = np.linalg.svd(xy - ctr, full_matrices=False)
    def to_uv(p):
        uv = (np.asarray(p)[:, :2] - ctr) @ Vt.T
        return uv[:, 0], uv[:, 1]
    u0, _ = to_uv(xy[[0, -1]])
    flip = u0[1] < u0[0]
    def to_uv_oriented(p):
        u, v = to_uv(p)
        return (-u if flip else u), v
    return to_uv_oriented


# ----------------------------------------------------------------------
def fig_registration():
    """Slice-based registration: match, anchor validation, datum correction,
    and the two-SLAM drift bound."""
    rep = json.load(open(OUT / "registration_report.json"))
    T = np.array(json.load(open(OUT_BASE / "transform_ed_slice.json"))["matrix"])
    val = json.load(open(OUT / "registration_validation.json"))
    chain = json.load(open(OUT_BASE / "transform_flf_ed_chain.json"))

    fig = plt.figure(figsize=(7.1, 7.0))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.35, 1, 1], hspace=0.55,
                          wspace=0.30, width_ratios=[1.6, 1])

    # (a) skylight openings in the surface frame with the registered interior
    # rim centroids. Constellation is elongated ~south-north, so northing runs
    # along the panel x-axis; the consensus-rejected fourth detection is inset.
    ax = fig.add_subplot(gs[0, 0])
    matched_aer = {p[1] for p in rep["landmark"]["inliers"]}
    cluster = []
    for h in rep["aerial_holes"]:
        cx, cy, _ = h["centroid"]
        if h["id"] not in matched_aer:
            continue
        cluster.append((cy, cx))
        ang = 90.0 - np.degrees(h.get("orientation", 0.0))
        e = Ellipse((cy, cx), 2 * max(h["semi_major"], 0.8),
                    2 * max(h["semi_minor"], 0.8), angle=ang,
                    fill=False, lw=1.4, edgecolor=C["surface"])
        ax.add_patch(e)
        ax.annotate(f"S{sorted(matched_aer).index(h['id']) + 1}",
                    (cy, cx), xytext=(0, 12), textcoords="offset points",
                    ha="center", fontsize=8)
    th = np.array([h["centroid"] for h in rep["tube_holes"]])
    th_a = (T @ np.column_stack([th, np.ones(len(th))]).T).T[:, :3]
    ax.scatter(th_a[:, 1], th_a[:, 0], marker="+", s=55, lw=1.4,
               color=C["tube"], label="interior rim centroid (registered)",
               zorder=5)
    ax.scatter([], [], marker="o", facecolor="none", edgecolor=C["surface"],
               s=55, label="surface-model opening")
    cl = np.array(cluster)
    ax.set_xlim(cl[:, 0].min() - 8, cl[:, 0].max() + 8)
    ax.set_ylim(cl[:, 1].min() - 10, cl[:, 1].max() + 12)
    ax.set_aspect("equal")
    ax.set_xlabel("northing (m, local frame)")
    ax.set_ylabel("easting (m)")
    ax.set_title("a  Skylight openings and registered rim centroids",
                 loc="left")
    ax.legend(frameon=False, loc="upper left", fontsize=6.5, ncols=1,
              handletextpad=0.4, borderaxespad=0.3)
    # location inset, boxed, in the empty lower-right corner: the three
    # accepted openings plus the consensus-rejected fourth detection
    axi = ax.inset_axes([0.70, 0.05, 0.28, 0.30])
    rej = [h for h in rep["aerial_holes"] if h["id"] not in matched_aer]
    axi.scatter(cl[:, 0], cl[:, 1], s=8, color=C["surface"], lw=0)
    for h in rej:
        axi.plot([h["centroid"][1]], [h["centroid"][0]], "x", ms=5,
                 color=C["flagged"])
    axi.set_aspect("equal")
    axi.margins(0.35)
    axi.set_xticks([]); axi.set_yticks([])
    for sp in axi.spines.values():
        sp.set_visible(True); sp.set_color("0.6"); sp.set_linewidth(0.6)
    axi.text(0.5, 0.86, "rejected 4th detection ($\\times$)",
             transform=axi.transAxes, fontsize=5.6, color="0.35",
             ha="center")

    # (b) anchor validation: signed offset between the through-skylight floor
    # seen by the photogrammetry at S3 and the lidar floor, recomputed here
    # for the histogram (stats cached in registration_validation.json).
    ax = fig.add_subplot(gs[0, 1])
    from scipy.spatial import cKDTree
    aer_holes = json.load(open(OUT_BASE / "aerial_holes.json"))
    s3 = np.array(aer_holes["holes"][2]["centroid"])

    def near_s3(q):
        return np.hypot(q[:, 0] - s3[0], q[:, 1] - s3[1]) < 5.0

    def load_near(p):
        parts = []
        with BinaryPcdReader(p) as r:
            for c in r.chunks():
                q = c[:, :3].astype(np.float64)
                parts.append(q[near_s3(q)])
        return np.vstack(parts)

    aer = load_near(MAPS / "aerial_10cm.pcd")
    tub = load_near(SECT_PCD)
    surf = np.percentile(aer[:, 2], 90)
    a_in = aer[aer[:, 2] < surf - 3.0]
    txy = cKDTree(tub[:, :2])
    dz = []
    for p in a_in:
        idx = txy.query_ball_point(p[:2], 0.5)
        if len(idx) < 5:
            continue
        z = tub[idx, 2]
        dz.append(p[2] - np.median(z[z < np.percentile(z, 30)]))
    dz = np.array(dz)
    ax.hist(dz, bins=25, color=C["tube"], alpha=0.85)
    st = val["floor_through_skylight_s3"]
    ax.axvline(0, color="k", lw=0.8)
    ax.axvline(st["median_m"], color=C["skylight"], lw=1.1, ls="--")
    ax.annotate(f"median {st['median_m']:.2f} m\nRMS {st['rms_m']:.2f} m",
                (0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
                fontsize=7, color="0.25")
    ax.set_xlabel("floor offset at S3 (m)")
    ax.set_ylabel("points")
    ax.set_title("b  Floor seen through S3\nvs interior lidar floor",
                 loc="left")

    # (c) the datum correction: interior ceiling under the rim-centroid
    # solution and under the slice registration, matched horizontally.
    base_rows = list(csv.DictReader(open(OUT_BASE / "roof_thickness.csv")))
    ed_rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    g = lambda rows, k: np.array([float(r[k]) if r[k] else np.nan
                                  for r in rows])
    bx, by, bz = g(base_rows, "x"), g(base_rows, "y"), g(base_rows, "z_ceil")
    es, ex, ey = g(ed_rows, "s"), g(ed_rows, "x"), g(ed_rows, "y")
    ez, ed_dem = g(ed_rows, "z_ceil"), g(ed_rows, "z_dem")
    t = cKDTree(np.c_[bx, by])
    dist, idx = t.query(np.c_[ex, ey])
    bz_m = np.where(dist < 2.0, bz[idx], np.nan)

    ax = fig.add_subplot(gs[1, :])
    ax.plot(es, ed_dem, lw=1.1, color=C["surface"], label="surface DEM")
    ax.plot(es, bz_m, lw=1.1, color=C["dlio"],
            label="ceiling, rim-centroid solve (4-DOF)")
    ax.plot(es, ez, lw=1.1, color=C["tube"],
            label="ceiling, slice registration")
    ax.set_ylabel("elevation (m)")
    ymax = np.nanmax(ed_dem)
    ax.set_ylim(None, ymax + 5.5)
    ax.set_title("c  Interior ceiling against the surface under the two "
                 "registrations", loc="left")
    ax.legend(frameon=False, loc="upper left", ncols=3, fontsize=6.5,
              handlelength=1.6, columnspacing=1.0, borderaxespad=0.2)

    # (d) two-SLAM consistency in the slice-registered frame (drift bound).
    zones = chain["ceiling_offset_validation"]
    ax = fig.add_subplot(gs[2, :])
    mids, means, stds, labels = [], [], [], []
    for k, v in zones.items():
        lo, hi = k[1:].split("-")
        mids.append((float(lo) + float(hi)) / 2)
        means.append(v["mean"]); stds.append(v["std"])
        labels.append(f"{lo}-{hi} m")
    ax.errorbar(mids, means, yerr=stds, fmt="o", ms=4, lw=1.1, capsize=3,
                color=C["envelope"])
    ax.axhline(0, color="0.5", lw=0.8, ls="--")
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel("distance along tube $s$ (m)")
    ax.set_ylabel(r"ceiling offset (m)")
    ax.set_title("d  Ceiling agreement between the two independent SLAM "
                 "solutions in the registered frame", loc="left")

    fig.savefig(FIGS / "ed_registration_panel.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  ed_registration_panel.pdf")


# ----------------------------------------------------------------------
def fig_consistency():
    s, x, y, z_dem, z_ceil, tau, ghost, klass = load_roof()

    fig = plt.figure(figsize=(7.1, 6.4))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.4, 1, 1.2], hspace=0.55)

    # (a) longitudinal profile: DEM vs measured ceiling, coloured by class
    ax = fig.add_subplot(gs[0])
    ax.plot(s, z_dem, lw=1.2, color=C["surface"], label="surface DEM")
    for cls, col in CLASS_COLOUR.items():
        m = klass == cls
        if m.any():
            ax.scatter(s[m], z_ceil[m], s=4, color=col, label=cls, lw=0)
    ax.set_ylabel("elevation (m)")
    ax.set_title("a  Surface elevation and interior ceiling along the tube",
                 loc="left")
    ax.legend(frameon=False, ncols=5, loc="lower left", fontsize=6.5)

    # (b) above-surface return fraction (the multipass screen)
    ax = fig.add_subplot(gs[1])
    for cls, col in CLASS_COLOUR.items():
        m = klass == cls
        if m.any():
            ax.scatter(s[m], np.clip(ghost[m], 0, 1), s=4, color=col, lw=0)
    ax.axhline(0.05, color="k", lw=0.9, ls="--")
    ax.annotate("5% screen threshold", (0.99, 0.06), xycoords=("axes fraction", "data"),
                ha="right", va="bottom", fontsize=7, color="0.25")
    ax.set_ylabel("interior returns above\nlocal surface (fraction)")
    ax.set_xlabel("distance along tube $s$ (m)")
    ax.set_title("b  Column consistency statistic per station", loc="left")

    # (c) station classification footprint (plan view, principal frame)
    ax = fig.add_subplot(gs[2])
    to_uv = principal_frame(np.column_stack([x, y]))
    u, v = to_uv(np.column_stack([x, y]))
    for cls, col in CLASS_COLOUR.items():
        m = klass == cls
        if m.any():
            ax.scatter(u[m], v[m], s=6, color=col, lw=0, label=cls)
    for tick in range(0, int(s.max()) + 1, 50):
        j = int(np.argmin(np.abs(s - tick)))
        ax.annotate(f"{tick} m", (u[j], v[j]), fontsize=6, color="0.35",
                    xytext=(2 if tick == 0 else 0, 7),
                    textcoords="offset points",
                    ha="left" if tick == 0 else "center")
    ax.set_aspect("equal")
    ax.margins(x=0.02, y=0.3)
    ax.set_xlabel("along principal axis (m)")
    ax.set_ylabel("across (m)")
    ax.set_title("c  Station classification for the roof statistics "
                 f"({int((klass=='intact').sum())} of {len(s)} retained)",
                 loc="left")

    fig.savefig(FIGS / "ed_consistency_panel.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  ed_consistency_panel.pdf")


# ----------------------------------------------------------------------
def fig_centreline():
    s, P, Tn = load_centreline()
    _, _, _, _, _, _, _, klass = load_roof()
    shell = load_xyz(SHELL_PCD, keep=350_000)

    fig = plt.figure(figsize=(7.1, 8.2))
    gs = fig.add_gridspec(4, 3, height_ratios=[1.55, 0.95, 1, 1],
                          hspace=0.5, wspace=0.3)

    to_uv = principal_frame(P[:, :2])
    su, sv = to_uv(shell)
    cu, cv = to_uv(P)

    # (a) plan view: shell + skeleton centreline
    ax = fig.add_subplot(gs[0, :])
    ax.scatter(su, sv, s=0.3, color=C["shell"], lw=0, rasterized=True)
    ax.plot(cu, cv, lw=1.6, color=C["tube"], label="L1-medial centreline")
    sky = klass == "skylight"
    ax.plot(cu[sky], cv[sky], lw=3.2, color=C["skylight"], label="skylights")
    for tick in range(0, int(s.max()) + 1, 50):
        j = int(np.argmin(np.abs(s - tick)))
        ax.annotate(f"{tick} m", (cu[j], cv[j]), fontsize=6, color="0.35",
                    xytext=(2 if tick == 0 else 0, 8),
                    textcoords="offset points",
                    ha="left" if tick == 0 else "center")
    ax.set_aspect("equal")
    ax.set_xlim(cu.min() - 25, cu.max() + 25)
    ax.set_ylim(cv.min() - 30, cv.max() + 30)
    ax.set_ylabel("across (m)")
    ax.set_title("a  Interior shell and extracted centreline, plan view",
                 loc="left")
    ax.legend(frameon=False, loc="lower right", fontsize=6.5)

    # (b) side view: shell + centreline elevation
    ax = fig.add_subplot(gs[1, :])
    ax.scatter(su, shell[:, 2], s=0.3, color=C["shell"], lw=0, rasterized=True)
    ax.plot(cu, P[:, 2], lw=1.6, color=C["tube"])
    ax.set_xlim(cu.min() - 25, cu.max() + 25)
    ax.set_ylim(P[:, 2].min() - 15, P[:, 2].max() + 15)
    ax.set_xlabel("along principal axis (m)")
    ax.set_ylabel("elevation (m)")
    rel = np.percentile(P[:, 2], 99) - np.percentile(P[:, 2], 1)
    net = P[-1, 2] - P[0, 2]
    trend = ("no net descent" if abs(net) < 3
             else ("net descent" if net < 0 else "net rise"))
    ax.set_title(f"b  Side view: close to horizontal, {rel:.0f} m relief, "
                 f"{trend}", loc="left")

    # (c-h) cross-section gallery from the 10 cm map
    cloud = load_xyz(SECT_PCD)
    from scipy.spatial import cKDTree
    tree = cKDTree(cloud)
    stations = [25, 75, 125, 175, 225, 275]
    letters = "cdefgh"
    morpho = {float(r["s"]): r for r in
              csv.DictReader(open(OUT / "morphometry.csv"))}
    for i, st in enumerate(stations):
        ax = fig.add_subplot(gs[2 + i // 3, i % 3])
        j = int(np.argmin(np.abs(s - st)))
        stn, T = P[j], Tn[j]
        idx = tree.query_ball_point(stn, 20.0)
        d = cloud[idx] - stn
        m = np.abs(d @ T) <= 0.375
        sec = d[m]
        uvec = np.cross(T, [0, 0, 1.0])
        uvec = uvec / max(np.linalg.norm(uvec), 1e-6)
        vvec = np.cross(T, uvec)
        if vvec[2] < 0:
            vvec = -vvec
        ax.scatter(sec @ uvec, sec @ vvec, s=0.5, color=C["tube"], lw=0,
                   rasterized=True)
        row = morpho.get(round(float(s[j]), 4)) or \
            morpho.get(min(morpho, key=lambda k: abs(k - s[j])))
        area = f"{float(row['area']):.0f}" if row and row["area"] else "-"
        ax.set_aspect("equal")
        ax.set_title(f"{letters[i]}  $s$={s[j]:.0f} m,  $A$={area} m$^2$",
                     loc="left")
        if i % 3 == 0:
            ax.set_ylabel("height (m)")
        if i // 3 == 1:
            ax.set_xlabel("across (m)")

    fig.savefig(FIGS / "ed_centreline_panel.pdf", bbox_inches="tight",
                dpi=300)
    plt.close(fig)
    print("  ed_centreline_panel.pdf")


if __name__ == "__main__":
    FIGS.mkdir(exist_ok=True)
    print("Extended Data panels ->", FIGS)
    fig_registration()
    fig_consistency()
    fig_centreline()
