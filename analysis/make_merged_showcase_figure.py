"""Session 2026-07-10: composite showcase figure of the co-registered
surface-subsurface model (upgrades fig:merged in the manuscript).

Panels:
  a  oblique cutaway render of the full-extent merged model
     (analysis_out/renders/showcase_oblique.png, from render_merged_showcase.py)
  b  plan-view roof-thickness map: DEM hillshade + per-cell tau raster
  c  geological long-section along the centreline: surface, ceiling, floor,
     roof fill coloured by tau, station classes honest (excluded = grey)
  d  transverse slice through the conduit at a skylight: both point clouds
     in section, with the roof thickness dimensioned

Colour logic is shared with the render: sequential inferno segment, BRIGHT =
THIN roof (the hazard end), DARK = THICK (the shielded end).

Output: <manuscript>/figures/merged_showcase_panel.pdf (+ .png preview)

Run:  env -u PYTHONPATH ./.venv/bin/python analysis/make_merged_showcase_figure.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LightSource, ListedColormap, Normalize
from matplotlib.patches import Ellipse
from scipy.spatial import cKDTree

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
import os
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out")))
TUBE_PCD = Path(os.environ.get("SHOWCASE_TUBE_PCD", str(ROOT / "maps/flf_10cm_aerial.pcd")))
FIGS = (ROOT.parent / "_Nature__Autonomous_aerial_reconnaissance_of_a_basaltic_"
        "lava_tube_reveals_interior_morphology_and_roof_thickness_distribution_"
        "for_planetary_subsurface" / "figures")

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42,
})

TAU_MAX = 16.5
# bright orange (thin) -> deep purple-black (thick); matches the 3D render
TAU_CMAP = ListedColormap(plt.get_cmap("inferno")(np.linspace(0.88, 0.14, 256)))
TAU_NORM = Normalize(0, TAU_MAX)
SKY_C = "#D55E00"
CONDUIT_C = "#c9d7e4"


# ---------------------------------------------------------------- data loading
def cols(path, *names):
    rows = list(csv.DictReader(open(path)))
    return [np.array([float(r[n]) for r in rows]) for n in names], rows


def load_xyz_rgb(path: Path):
    xyz_parts, rgb_parts = [], []
    with BinaryPcdReader(path) as r:
        has_rgb = "rgb" in r.fields
        rgb_i = r.fields.index("rgb") if has_rgb else -1
        for chunk in r.chunks():
            xyz_parts.append(chunk[:, :3].astype(np.float64))
            if has_rgb:
                packed = chunk[:, rgb_i].astype(np.float32).view(np.uint32)
                rgb_parts.append(np.stack([(packed >> 16) & 255, (packed >> 8) & 255,
                                           packed & 255], 1).astype(np.float64) / 255.0)
    return np.vstack(xyz_parts), (np.vstack(rgb_parts) if rgb_parts else None)


(rs, rx, ry, rzdem, rzceil, rtau, rspan), rrows = cols(
    OUT / "roof_thickness.csv", "s", "x", "y", "z_dem", "z_ceil", "tau", "span")
rclass = np.array([r["class"] for r in rrows])
(ms, mzfloor, mzceil), _ = cols(OUT / "morphometry.csv", "s", "z_floor", "z_ceil")
(cs, cx, cy, cz, ctx, cty), _ = cols(OUT / "centreline.csv", "s", "x", "y", "z", "tx", "ty")
holes = json.load(open(OUT / "aerial_holes.json"))["holes"][:3]

dem = np.load(OUT / "dem_grid_full_surface_and_subsurface_merged.npz")
DEMZ, XMIN, YMIN, RES = dem["z"], float(dem["xmin"]), float(dem["ymin"]), float(dem["res"])
NX, NY = DEMZ.shape   # DEM cache is x-major: DEMZ[ix, iy]

# principal axis of the centreline -> rotate the map so the tube runs left-right
d0 = np.stack([cx, cy], 1) - [cx.mean(), cy.mean()]
u, sv, vt = np.linalg.svd(d0, full_matrices=False)
ax0 = vt[0] / np.linalg.norm(vt[0])
if (d0[-1] - d0[0]) @ ax0 < 0:
    ax0 = -ax0  # s increases to the right
ROT = np.array([[ax0[0], ax0[1]], [-ax0[1], ax0[0]]])


def rot(x, y):
    p = ROT @ np.stack([x - cx.mean(), y - cy.mean()])
    return p[0], p[1]


def dem_at(x, y):
    ix = np.clip(((x - XMIN) / RES).astype(np.int64), 0, NX - 1)
    iy = np.clip(((y - YMIN) / RES).astype(np.int64), 0, NY - 1)
    return DEMZ[ix, iy]


# ------------------------------------------------------------------- panel b
def panel_b(ax):
    """Plan-view roof-thickness map on DEM hillshade, rotated frame."""
    crx, cry = rot(cx, cy)
    pad = 32.0
    gx = np.arange(crx.min() - pad, crx.max() + pad, 1.0)
    gy = np.arange(cry.min() - pad, cry.max() + pad, 1.0)
    GX, GY = np.meshgrid(gx, gy)
    # inverse-rotate the grid into the aerial frame and sample the DEM
    inv = ROT.T @ np.stack([GX.ravel(), GY.ravel()])
    QX, QY = inv[0] + cx.mean(), inv[1] + cy.mean()
    Z = dem_at(QX, QY).reshape(GX.shape)
    ls = LightSource(azdeg=315, altdeg=50)
    hs = ls.hillshade(Z, vert_exag=1.5, dx=1.0, dy=1.0)
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1.35, origin="lower",
              extent=[gx[0], gx[-1], gy[0], gy[-1]], interpolation="bilinear")

    # per-cell roof thickness from the FAST-LIO tube ceiling under the DEM
    txyz, _ = load_xyz_rgb(TUBE_PCD)
    tzdem = dem_at(txyz[:, 0], txyz[:, 1])
    below = txyz[:, 2] < tzdem - 0.3
    txyz = txyz[below]
    trx, try_ = rot(txyz[:, 0], txyz[:, 1])
    ix = ((trx - gx[0]) / 1.0).astype(np.int64)
    iy = ((try_ - gy[0]) / 1.0).astype(np.int64)
    ok = (ix >= 0) & (ix < len(gx)) & (iy >= 0) & (iy < len(gy))
    ceil = np.full(GX.shape, -np.inf)
    np.maximum.at(ceil, (iy[ok], ix[ok]), txyz[ok, 2])
    cnt = np.zeros(GX.shape)
    np.add.at(cnt, (iy[ok], ix[ok]), 1.0)
    tau_cell = Z - ceil
    tau_cell[(cnt < 4) | ~np.isfinite(tau_cell)] = np.nan

    # honesty mask: distinct flat tint where the nearest station is excluded
    # by the consistency screen (mapped, but not part of the tau statistics)
    srx, sry = rot(rx, ry)
    _, near = cKDTree(np.stack([srx, sry], 1)).query(
        np.stack([GX.ravel(), GY.ravel()], 1), workers=-1)
    excl = np.isin(rclass[near].reshape(GX.shape), ("multipass", "inconsistent"))
    have = np.isfinite(tau_cell)
    excl_img = np.zeros(GX.shape + (4,))
    excl_img[excl & have] = matplotlib.colors.to_rgba("#8fb4c9", 0.95)
    ax.imshow(excl_img, origin="lower",
              extent=[gx[0], gx[-1], gy[0], gy[-1]], interpolation="nearest")
    shown = np.ma.masked_invalid(np.where(excl, np.nan, tau_cell))
    im = ax.imshow(np.clip(shown, 0, TAU_MAX), cmap=TAU_CMAP, norm=TAU_NORM,
                   origin="lower", extent=[gx[0], gx[-1], gy[0], gy[-1]],
                   interpolation="nearest")
    ax.plot([], [], marker="s", ms=6, ls="none", mfc="#8fb4c9", mec="none",
            label="mapped, excluded by screen")
    ax.legend(loc="lower right", frameon=False, fontsize=6.5,
              handletextpad=0.4, borderaxespad=0.2)

    ax.plot(crx, cry, color="w", lw=0.6, alpha=0.85)
    for j, h in enumerate(holes):
        hxr, hyr = rot(np.array([h["centroid"][0]]), np.array([h["centroid"][1]]))
        ang = np.degrees(h["orientation"]) - np.degrees(np.arctan2(ax0[1], ax0[0]))
        ax.add_patch(Ellipse((hxr[0], hyr[0]), 2 * h["semi_major"] + 4,
                             2 * h["semi_minor"] + 4, angle=ang, fill=False,
                             ec=SKY_C, lw=1.1))
        ax.annotate(f"S{j+1}", (hxr[0], hyr[0]), textcoords="offset points",
                    xytext=(8, 9), color=SKY_C, fontsize=8, fontweight="bold")
    # arc-length ticks every 50 m
    for sk in range(0, int(cs.max()) + 1, 50):
        i = np.argmin(np.abs(cs - sk))
        ax.annotate(f"{sk}", (crx[i], cry[i]), textcoords="offset points",
                    xytext=(0, -11), fontsize=6, color="0.25", ha="center")
        ax.plot(crx[i], cry[i], marker="o", ms=2, color="w", mec="0.2", mew=0.4)
    # scale bar
    x0, y0 = crx.min() + 5, cry.min() - pad + 12
    ax.plot([x0, x0 + 50], [y0, y0], color="k", lw=1.6)
    ax.text(x0 + 25, y0 + 4, "50 m", ha="center", fontsize=7)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.033, pad=0.01, aspect=14,
                      ticks=[0, 4, 8, 12, 16])
    cb.set_label("overburden $\\tau$ (m)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    return gx, gy


# ------------------------------------------------------------------- panel c
def panel_c(ax):
    """Geological long-section along the centreline."""
    intact = rclass == "intact"
    zfloor = np.interp(rs, ms, mzfloor)
    # conduit void
    ax.fill_between(rs, zfloor, rzceil, color=CONDUIT_C, lw=0, label="conduit void")
    # roof fill coloured by tau where intact, grey where excluded
    for i in range(len(rs) - 1):
        sseg = [rs[i], rs[i + 1]]
        c = (TAU_CMAP(TAU_NORM(min(rtau[i], TAU_MAX)))
             if intact[i] else (0.78, 0.78, 0.78, 1.0))
        ax.fill_between(sseg, [rzceil[i], rzceil[i + 1]],
                        [rzdem[i], rzdem[i + 1]], color=c, lw=0)
    ax.plot(rs, rzdem, color="k", lw=0.9)
    ax.plot(rs, rzceil, color="0.15", lw=0.5)
    ax.plot(rs, zfloor, color="0.4", lw=0.5)
    # skylights
    for j, h in enumerate(holes):
        i = np.argmin(np.hypot(rx - h["centroid"][0], ry - h["centroid"][1]))
        ax.annotate(f"S{j+1}", (rs[i], rzdem[i] + 2.2), ha="center",
                    color=SKY_C, fontsize=8, fontweight="bold")
        ax.plot([rs[i]], [rzdem[i] + 0.8], marker="v", ms=4, color=SKY_C)
    ax.text(152, rzdem[np.argmin(np.abs(rs - 152))] + 3.0, "surface (DEM)",
            fontsize=7)
    ax.text(178, 229.5, "conduit void", fontsize=7, color="0.3")
    ax.annotate("roof, coloured by $\\tau$", (300, 242.5), (248, 248.5),
                fontsize=7, arrowprops=dict(arrowstyle="-", color="0.3", lw=0.6))
    ax.annotate("mapped, excluded by screen", (55, 238.5), (8, 246.5),
                fontsize=6.5, color="0.35",
                arrowprops=dict(arrowstyle="-", color="0.5", lw=0.6))
    ax.set_xlim(0, rs[-1])
    ax.set_ylim(zfloor.min() - 1.5, rzdem.max() + 5.5)
    ax.set_xlabel("along-tube distance $s$ (m)")
    ax.set_ylabel("elevation (m, local)")
    ax.text(0.995, 0.03, "vertical exaggeration 3$\\times$",
            transform=ax.transAxes, ha="right", fontsize=6.5, color="0.35")
    ax.set_aspect(3.0)


# ------------------------------------------------------------------- panel d
def panel_d(ax):
    """Transverse slice at the thickest intact roof: both clouds in section,
    with the roof slab dimensioned. The rock between the photogrammetric
    surface and the lidar ceiling is exactly tau."""
    # thickest intact roof away from the survey end (end stations have a
    # partially observed ring)
    intact = (rclass == "intact") & (rs >= 240) & (rs <= 295)
    i_r = int(np.argmax(np.where(intact, rtau, -np.inf)))
    i = np.argmin(np.abs(cs - rs[i_r]))
    p0 = np.array([cx[i], cy[i]])
    t = np.array([ctx[i], cty[i]]); t /= np.linalg.norm(t)
    n = np.array([-t[1], t[0]])
    LAT = 20.0

    def slab(xyz, half, rgb=None):
        d_al = (xyz[:, :2] - p0) @ t
        d_lat = (xyz[:, :2] - p0) @ n
        keep = (np.abs(d_al) < half) & (np.abs(d_lat) < LAT)
        return d_lat[keep], xyz[keep, 2], (rgb[keep] if rgb is not None else None)

    sxyz, srgb = load_xyz_rgb(ROOT / "maps/merged_surface_10cm.pcd")
    # keep only the surface skin (the split leaves stray deep points near DEM gaps)
    skin = sxyz[:, 2] > dem_at(sxyz[:, 0], sxyz[:, 1]) - 2.0
    su, sz, sc = slab(sxyz[skin], 2.5, srgb[skin])
    txyz, _ = load_xyz_rgb(TUBE_PCD)
    tu, tz, _ = slab(txyz, 1.5)
    ax.scatter(su, sz, s=1.6, c=np.clip(sc * 1.05, 0, 1), lw=0, rasterized=True)
    ax.scatter(tu, tz, s=1.6, color="#0072B2", lw=0, rasterized=True)
    # dimension tau at the crown
    zc, zd = rzceil[i_r], rzdem[i_r]
    uc = 0.0
    ax.annotate("", (uc, zd), (uc, zc),
                arrowprops=dict(arrowstyle="<->", color="k", lw=1.0))
    ax.text(uc + 1.6, (zc + zd) / 2,
            f"$\\tau$ = {rtau[i_r]:.1f} m", fontsize=8, ha="left")
    ax.text(0.03, 0.06, f"$s$ = {rs[i_r]:.0f} m", transform=ax.transAxes,
            fontsize=7, va="bottom", color="0.3")
    ax.text(-LAT + 1.5, sz.max() - 0.6, "surface\n(photogrammetry)", fontsize=6.5,
            color="0.25", va="top")
    ax.text(LAT - 1.5, tz.min() + 2.0, "conduit\n(interior lidar)", fontsize=6.5,
            color="#0072B2", ha="right")
    ax.set_xlabel("across-tube distance (m)")
    ax.set_ylabel("elevation (m, local)")
    ax.set_aspect("equal")
    ax.set_xlim(-LAT, LAT)
    ax.set_ylim(np.percentile(tz, 0.5) - 1.0, sz.max() + 2.0)


# ------------------------------------------------------------------- assembly
def main() -> None:
    fig = plt.figure(figsize=(7.2, 8.6))
    gs = fig.add_gridspec(3, 5, height_ratios=[1.45, 0.95, 1.0],
                          hspace=0.28, wspace=0.9)
    axa = fig.add_subplot(gs[0, :])
    axb = fig.add_subplot(gs[1, :])
    axc = fig.add_subplot(gs[2, :3])
    axd = fig.add_subplot(gs[2, 3:])

    hero = plt.imread(OUT / "renders/showcase_oblique.png")
    # trim white margins
    nonwhite = np.where((hero[..., :3].min(-1) < 0.97).any(1))[0]
    nonwhite_c = np.where((hero[..., :3].min(-1) < 0.97).any(0))[0]
    hero = hero[nonwhite[0]:nonwhite[-1] + 1, nonwhite_c[0]:nonwhite_c[-1] + 1]
    axa.imshow(hero, interpolation="bilinear")
    axa.set_xticks([]); axa.set_yticks([])
    for sp in axa.spines.values():
        sp.set_visible(False)

    panel_b(axb)
    panel_c(axc)
    panel_d(axd)

    for ax, lab in ((axa, "a"), (axb, "b"), (axc, "c"), (axd, "d")):
        ax.text(-0.015, 1.02, lab, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom", ha="right")

    FIGS.mkdir(exist_ok=True)
    for ext, dpi in (("pdf", 300), ("png", 220)):
        fig.savefig(FIGS / f"merged_showcase_panel.{ext}", dpi=dpi,
                    bbox_inches="tight")
    print("wrote", FIGS / "merged_showcase_panel.pdf")


if __name__ == "__main__":
    main()
