"""Nature initial submission: the five main figures (plan/13 Step 4).

Reuses the panel code of make_figures.py (roof/morphometry helpers, the
planetary figure) and make_merged_showcase_figure.py (plan view, long section,
transverse slice) and the showcase render, reading every number from
ANALYSIS_OUT (default analysis_out_v9). Nothing is hand-typed: values that
appear on the figures come from paper_numbers.json, review_stats.json or the
CSVs, and anything computed here is printed so it can be checked.

Nature rules applied: 180 mm wide, at most 170 mm tall; Nimbus Sans
(Helvetica clone) embedded; 7 pt text, 6 pt ticks, bold 8 pt panel letters;
axis labels lower-case with a first capital and a space before the unit;
thousands with commas; keys on the figure only where they aid reading.

Outputs (PDF + 300 dpi PNG) in analysis_out_v10/figures_nature/:
    fig1_anatomy_site   tube anatomy schematic (vector) + site map on the ortho
    fig2_survey         photographs, map growth, global graph (rasters)
    fig3_merged_model   co-registered model (render, plan, long section, slice)
    fig4_geometry_cover along-tube geometry and cover (six panels)
    fig5_planetary      planetary comparison (fig_planetary of make_figures)

Run from the repo root:
    env -u PYTHONPATH PYTHONPATH=src .venv/bin/python analysis/make_nature_figures.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
OUT10 = ROOT / "analysis_out_v10"
FIGS = OUT10 / "figures_nature"
SCRATCH = OUT10 / "figures_nature_scratch"   # where the reused scripts drop their own copies
MS_FIGS = (ROOT.parent / "_Nature__Autonomous_aerial_reconnaissance_of_a_basaltic_"
           "lava_tube_reveals_interior_morphology_and_roof_thickness_distribution_"
           "for_planetary_subsurface" / "figures")   # read-only source of rasters
ORTHO = Path("/home/aakapatel/Nature Journal Article Autonomous lava tube mapping/"
             "Birgir_data_and_papers_08_sep/20250828_Raufarholshellir_ORTHO_ISN16.tif")
ORIGIN = (479158.0, 7089826.0)          # local frame = EPSG:32627 minus this
MAPS = ROOT / "maps"

# the reused modules read their inputs from the environment at import time
os.environ.setdefault("ANALYSIS_OUT", str(OUT))
os.environ["FIGS_DIR"] = str(SCRATCH)
os.environ.setdefault("ROOF_DEM_PCD", str(MAPS / "aerial_isn16_ortho_dem.pcd"))
os.environ.setdefault("SHOWCASE_TUBE_PCD", str(MAPS / "flf_10cm_ed_ortho.pcd"))
os.environ.setdefault("SHOWCASE_SURFACE_PCD", str(MAPS / "aerial_isn16_ortho_surface_10cm.pcd"))
os.environ.setdefault("EDFIG_SHELL_PCD", str(MAPS / "flf_30cm_ed_ortho.pcd"))
os.environ.setdefault("EDFIG_SECT_PCD", str(MAPS / "flf_10cm_ed_ortho.pcd"))
os.environ.setdefault("EDFIG_AERIAL_PCD", str(MAPS / "aerial_isn16_ortho_10cm.pcd"))
SCRATCH.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "analysis"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Polygon, Circle, Rectangle, FancyArrowPatch
from matplotlib.ticker import FuncFormatter
import matplotlib.patheffects as pe

import make_figures as mf                      # noqa: E402  (helpers + planetary)
import make_merged_showcase_figure as msf      # noqa: E402  (plan, long section, slice)
import make_ed_figures as edf                  # noqa: E402  (point-cloud loaders)
# the two panels promoted from Extended Data (the skylight constellation of
# Fig. 3a and the onboard-against-offline zones of Fig. 2d) are drawn by the
# Extended Data script, so there is one implementation of each
import make_nature_ed_figures as ned           # noqa: E402

MM = 1.0 / 25.4
# ned and the modules it imports set their own rcParams at import time; the
# Nature style below is applied afterwards and is the one that takes effect
plt.rcParams.update({
    "font.family": "Nimbus Sans",
    "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
    "legend.fontsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False,
    # Type 42 (CID, as in the Extended Data script). Type 3 cannot be written
    # once make_nature_ed_figures has been imported: its TextToPath instance
    # leaves the NimbusSans glyph map in a state matplotlib's Type 3 subsetter
    # rejects ("bytes must be in range(0, 256)").
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Nimbus Sans",
    "mathtext.it": "Nimbus Sans:italic", "mathtext.bf": "Nimbus Sans:bold",
    "savefig.dpi": 300,
})
C = mf.C
ELEV_LABEL = "Elevation (m a.s.l., ISH2004)"
_datum = json.loads((OUT / "datum.json").read_text())
assert _datum["elevation_axis_label"].lower() == ELEV_LABEL.lower(), _datum

PN = json.loads((OUT / "paper_numbers.json").read_text())
REV = json.loads((OUT / "review_stats.json").read_text())
HOLES = json.loads((OUT / "aerial_holes.json").read_text())["holes"][:3]


# ------------------------------------------------------------------ helpers
def mm_axes(fig, x, y, w, h, **kw):
    """Axes placed in millimetres from the lower-left corner of the figure."""
    W, H = fig.get_size_inches() / MM
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def comma_fmt(ax, which="x"):
    axis = ax.xaxis if which == "x" else ax.yaxis
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))


def capitalise_labels(fig):
    """Nature: lower-case labels with the first letter capitalised."""
    for ax in fig.axes:
        for get, set_ in ((ax.get_xlabel, ax.set_xlabel), (ax.get_ylabel, ax.set_ylabel)):
            t = get()
            if t and t[0].islower():
                set_(t[0].upper() + t[1:])


def clamp_fonts(fig, max_pt=7.0):
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_fontsize() > max_pt and t.get_fontweight() != "bold":
            t.set_fontsize(max_pt)


def place_letters(fig, panels):
    """panels: list of (axes or [axes], letter, group or None, (dx, dy) mm).
    The letter sits above the upper-left corner of the panel's tight bbox
    (ticks and labels included); panels in the same group share the
    left-most x so a column of letters lines up."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.get_size_inches() * fig.dpi
    boxes = []
    for axes, letter, group, off in panels:
        axes = axes if isinstance(axes, (list, tuple)) else [axes]
        bb = matplotlib.transforms.Bbox.union([a.get_tightbbox(r) for a in axes])
        boxes.append([bb.x0 / W, bb.y1 / H, letter, group, off])
    for g in {b[3] for b in boxes if b[3] is not None}:
        x0 = min(b[0] for b in boxes if b[3] == g)
        for b in boxes:
            if b[3] == g:
                b[0] = x0
    Wmm, Hmm = fig.get_size_inches() / MM
    for x, y, letter, _, (dx, dy) in boxes:
        fig.text(x + dx / Wmm, y + (0.6 + dy) / Hmm, letter, fontsize=8,
                 fontweight="bold", ha="left", va="bottom")


def save(fig, name):
    fig.savefig(FIGS / f"{name}.pdf")
    fig.savefig(FIGS / f"{name}.png", dpi=300)
    w, h = fig.get_size_inches() / MM
    print(f"{name}: {w:.0f} x {h:.0f} mm")
    plt.close(fig)


def halo(lw=1.6, color="w"):
    return [pe.withStroke(linewidth=lw, foreground=color)]


# ================================================================== figure 1
def anatomy_panel(ax):
    """Lava-tube anatomy, long section (own vector art after lt_fix.jpg)."""
    ax.set_axis_off()
    SOIL, ROOF, LINE, BED = "#d8c8a4", "#7a675a", "#5b4b41", "#6c5b50"
    VOID, SNOW, BLOCK, SKY = "#eaf0f4", "#ffffff", "#a9a198", "0.25"
    ax.set_xlim(-2, 103)
    ax.set_ylim(-1, 40)
    xs = np.linspace(0, 100, 401)
    z_s = 25 + 0.7 * np.sin(xs / 9.0) + 0.35 * np.sin(xs / 3.7 + 1)   # surface
    z_c = 15 + 0.9 * np.sin(xs / 7.0 + 1)                             # ceiling
    z_f = 6 + 0.5 * np.sin(xs / 5.0)                                  # floor
    xl, xr = 30.0, 38.0                     # skylight aperture at the ceiling
    xt = 87.0                               # collapse trench begins

    def seg(mask):
        return xs[mask], z_s[mask], z_c[mask], z_f[mask]

    # bedrock below the floor, everywhere
    ax.fill_between(xs, -1, z_f, color=BED, lw=0)
    # roof left of the skylight and right of it up to the trench; the rim is
    # flared (wider at the surface than at the ceiling)
    zs_at = lambda q: float(np.interp(q, xs, z_s))
    zc_at = lambda q: float(np.interp(q, xs, z_c))
    for lo, hi, lo_c, hi_c in ((0.0, xl - 1.2, 0.0, xl), (xr + 1.2, xt, xr, xt)):
        ms_ = (xs >= lo) & (xs <= hi)
        mc_ = (xs >= lo_c) & (xs <= hi_c)
        top = np.column_stack([xs[ms_], z_s[ms_]])
        bot = np.column_stack([xs[mc_][::-1], z_c[mc_][::-1]])
        poly = np.vstack([top, bot])
        ax.add_patch(Polygon(poly, closed=True, fc=ROOF, ec=LINE, lw=0.5))
        ax.fill_between(xs[ms_], z_s[ms_] - 1.6, z_s[ms_], color=SOIL, lw=0)
        ax.plot(xs[ms_], z_s[ms_], color=LINE, lw=0.6)
        for k, d in enumerate((3.2, 5.6, 7.6)):
            ax.plot(xs[ms_], z_s[ms_] - d + 0.25 * np.sin(xs[ms_] / 4 + k), color=LINE,
                    lw=0.35, alpha=0.8)
    # trench: surface drops to a rubble floor; the conduit opens into it
    ztr = 9.0
    ax.add_patch(Polygon([[xt, zs_at(xt)], [xt, zc_at(xt)], [xt, z_f[np.argmin(np.abs(xs - xt))]],
                          [xt + 4, ztr], [100, ztr + 0.8], [100, -1], [xt, -1]],
                         closed=True, fc=BED, ec="none"))
    ax.fill_between(xs[xs <= xt], z_f[xs <= xt], z_c[xs <= xt], color=VOID, lw=0)
    ax.plot([xt, xt], [zs_at(xt), zc_at(xt)], color=LINE, lw=0.6)
    ax.plot([xt, xt + 4, 100], [z_f[np.argmin(np.abs(xs - xt))], ztr, ztr + 0.8],
            color=LINE, lw=0.6)
    # snow cone under the skylight, breakdown blocks on the floor
    ax.add_patch(Polygon([[26, 6.3], [34, 13.5], [42.5, 6.2]], fc=SNOW, ec="0.4", lw=0.5))
    rng = np.random.default_rng(3)
    for bx, by, r in ((12, 6.6, 1.1), (17, 6.2, 0.8), (49, 6.4, 1.2), (55, 6.0, 0.7),
                      (63, 6.3, 0.9), (74, 6.1, 1.3), (80, 6.4, 0.8), (90, 9.6, 1.0),
                      (94, 9.9, 1.4), (98, 9.8, 0.8), (45, 6.1, 0.6), (24, 6.0, 0.7)):
        ang = rng.uniform(0, 2 * np.pi, 5)
        ang.sort()
        pts = np.column_stack([bx + r * np.cos(ang), by + 0.75 * r * np.sin(ang) + 0.5 * r])
        ax.add_patch(Polygon(pts, fc=BLOCK, ec="0.3", lw=0.4))
    # orbiter: body, panels, view cone to the aperture
    ox, oz = 34.0, 36.5
    ax.add_patch(Rectangle((ox - 0.9, oz - 0.6), 1.8, 1.2, fc=SKY, ec="none"))
    ax.plot([ox - 4.2, ox - 0.9], [oz, oz], color=SKY, lw=2.2)
    ax.plot([ox + 0.9, ox + 4.2], [oz, oz], color=SKY, lw=2.2)
    ax.plot([ox, xl - 1.2], [oz - 0.6, 25.4], color=SKY, lw=0.5, ls=(0, (2, 2)))
    ax.plot([ox, xr + 1.2], [oz - 0.6, 25.4], color=SKY, lw=0.5, ls=(0, (2, 2)))
    # robot in the conduit with lidar rays
    rx, rz = 62.0, 10.5
    ax.plot([rx - 1.6, rx + 1.6], [rz, rz], color=SKY, lw=1.6)
    for dx in (-1.6, 1.6):
        ax.add_patch(Ellipse((rx + dx, rz + 0.35), 1.6, 0.5, fc="none", ec=SKY, lw=0.6))
    for ang in np.linspace(0.35, np.pi - 0.35, 6):
        L = 30.0
        px, pz = rx + L * np.cos(ang), rz + L * np.sin(ang)
        tt = np.linspace(0, 1, 200)
        qx, qz = rx + tt * (px - rx), rz + tt * (pz - rz)
        inside = (qz <= np.interp(qx, xs, z_c)) & (qx >= 0) & (qx <= xt)
        k = np.argmax(~inside) if (~inside).any() else len(tt) - 1
        ax.plot(qx[:k], qz[:k], color=SKY, lw=0.4, ls=(0, (1.5, 2)), alpha=0.9)
    # overburden arrow at x = 72
    i = np.argmin(np.abs(xs - 72))
    ax.annotate("", (72, z_s[i]), (72, z_c[i]),
                arrowprops=dict(arrowstyle="<->", color="k", lw=0.7, shrinkA=0, shrinkB=0))
    ax.text(73.2, 0.5 * (z_s[i] + z_c[i]), "Overburden $\\tau$", fontsize=7, va="center")
    # labels: dark text over the white sky and the pale void, white text on rock
    kw = dict(fontsize=6.5, color="0.15")
    ax.text(2, 27.3, "Surface", **kw)
    ax.annotate("Soil, tephra", (14, z_s[np.argmin(np.abs(xs - 14))] - 0.8), (12, 31.5),
                fontsize=6.5, color="0.15", ha="center",
                arrowprops=dict(arrowstyle="-", color="0.3", lw=0.5))
    ax.text(4, 19.8, "Basalt roof\n(stacked flow units)", fontsize=6.5, color="w", va="center")
    ax.text(14, 11.5, "Drained\nconduit", fontsize=6.5, color="0.2", ha="center", va="center")
    ax.text(43.5, 9.3, "Snow cone", fontsize=6.5, color="0.2", ha="left", va="center")
    ax.annotate("Breakdown blocks", (24, 6.4), (30, 2.2), fontsize=6.5, color="w",
                ha="center", arrowprops=dict(arrowstyle="-", color="w", lw=0.5))
    ax.text(41.5, 29.0, "Skylight: the one feature seen\nfrom above and from within",
            fontsize=6.5, ha="left", va="center", color="0.15")
    ax.text(ox + 5.0, oz, "What orbit sees: the aperture", fontsize=6.5, va="center")
    ax.annotate("What the robot sees: the void from inside", (rx, rz - 0.6), (73, 2.2),
                fontsize=6.5, color="w", ha="center",
                arrowprops=dict(arrowstyle="-", color="w", lw=0.5))
    ax.text(93.5, 27.5, "Collapse trench\n(entrance)", fontsize=6.5, ha="center", va="bottom",
            color="0.15")
    ax.text(1.0, 2.2, "Basalt", fontsize=6.5, color="w")


def cross_section_panel(ax):
    """Small transverse section defining span L and the ceiling."""
    ax.set_axis_off()
    SOIL, ROOF, LINE, BED, VOID = "#d8c8a4", "#7a675a", "#5b4b41", "#6c5b50", "#eaf0f4"
    ax.set_xlim(-16, 16)
    ax.set_ylim(-1, 27)
    xs = np.linspace(-16, 16, 200)
    z_s = 22 + 0.4 * np.sin(xs / 3.0)
    # conduit outline: arched ceiling, flatter floor
    th = np.linspace(0, np.pi, 100)
    L = 18.0
    cx_, cz_ = (L / 2) * np.cos(th), 6.5 + 7.0 * np.sin(th) ** 0.9
    fx, fz = np.linspace(L / 2, -L / 2, 60), 6.5 - 1.2 * np.cos(np.linspace(-1.2, 1.2, 60))
    void = np.vstack([np.column_stack([cx_, cz_]), np.column_stack([fx, fz])])
    ax.fill_between(xs, -1, z_s, color=BED, lw=0)
    ax.fill_between(xs, z_s - 6.8, z_s, color=ROOF, lw=0)
    ax.fill_between(xs, z_s - 1.4, z_s, color=SOIL, lw=0)
    ax.plot(xs, z_s, color=LINE, lw=0.6)
    for d in (3.0, 5.0):
        ax.plot(xs, z_s - d, color=LINE, lw=0.35, alpha=0.8)
    ax.add_patch(Polygon(void, closed=True, fc=VOID, ec=LINE, lw=0.6))
    # span L at the widest, overburden at the crown
    yL = 6.9
    ax.annotate("", (-L / 2, yL), (L / 2, yL),
                arrowprops=dict(arrowstyle="<->", color="k", lw=0.7, shrinkA=0, shrinkB=0))
    ax.text(0, yL + 0.7, "Span $L$", fontsize=7, ha="center", va="bottom")
    zc = 13.5
    ax.annotate("", (0, float(np.interp(0, xs, z_s))), (0, zc),
                arrowprops=dict(arrowstyle="<->", color="k", lw=0.7, shrinkA=0, shrinkB=0))
    ax.text(0.8, 0.5 * (zc + 22), "$\\tau$", fontsize=7, va="center")
    ax.text(0, 25.2, "Cross-section", fontsize=6.5, ha="center", color="0.15")
    ax.text(0, 11.0, "Highest ceiling", fontsize=6, ha="center", color="0.3")


def iceland_coast():
    """Iceland coastline: Natural Earth 50 m, cached as a small JSON extract."""
    small = OUT10 / "iceland_coastline_ne50m.json"
    if small.exists():
        return json.loads(small.read_text())
    big = OUT10 / "ne_50m_coastline.geojson"
    if not big.exists():
        return None
    feats = json.loads(big.read_text())["features"]
    lines = []
    for f in feats:
        b = f.get("bbox") or [0, 0, 0, 0]
        if -25 < b[0] and b[2] < -13 and 63 < b[1] and b[3] < 67:
            g = f["geometry"]
            parts = g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]
            lines.extend(parts)
    small.write_text(json.dumps(lines))
    return lines


def site_map(fig, x0, y0, w, h):
    """Orthomosaic crop with the centreline, skylights, entrance, road."""
    import rasterio
    from pyproj import Transformer

    xlim, ylim = (1195.0, 1455.0), (545.0, 870.0)
    step = 0.25
    tr = Transformer.from_crs(32627, 8088, always_xy=True)
    xs = np.arange(xlim[0], xlim[1], step)
    ys = np.arange(ylim[0], ylim[1], step)
    XX, YY = np.meshgrid(xs, ys)
    Xp, Yp = tr.transform(XX.ravel() + ORIGIN[0], YY.ravel() + ORIGIN[1])
    with rasterio.open(ORTHO) as r:
        assert r.crs.to_epsg() == 8088, r.crs
        rr, cc = rasterio.transform.rowcol(r.transform, Xp, Yp)
        rr, cc = np.asarray(rr), np.asarray(cc)
        r0, c0, dec = rr.min(), cc.min(), 4      # read at ~12 cm, sample at 25 cm
        hh, ww = rr.max() - r0 + 1, cc.max() - c0 + 1
        tile = r.read([1, 2, 3], window=rasterio.windows.Window(c0, r0, ww, hh),
                      out_shape=(3, hh // dec, ww // dec))
    img = tile[:, np.clip((rr - r0) // dec, 0, tile.shape[1] - 1),
               np.clip((cc - c0) // dec, 0, tile.shape[2] - 1)]
    img = img.reshape(3, len(ys), len(xs)).transpose(1, 2, 0)

    ax = mm_axes(fig, x0, y0, w, h)
    ax.imshow(img, extent=[xlim[0], xlim[1], ylim[0], ylim[1]], origin="lower",
              interpolation="bilinear")
    # road (OpenStreetMap way 620673065, nodes lon/lat -> local frame)
    d = json.loads((OUT / "osm_way_620673065.json").read_text())
    nodes = {e["id"]: (e["lon"], e["lat"]) for e in d["elements"] if e["type"] == "node"}
    way = [e for e in d["elements"] if e["type"] == "way"][0]
    tr2 = Transformer.from_crs(4326, 32627, always_xy=True)
    road = np.array([tr2.transform(*nodes[n]) for n in way["nodes"]]) - np.array(ORIGIN)
    ax.plot(road[:, 0], road[:, 1], color="k", lw=0.6, ls=(0, (4, 2)))
    # 1970 map-sheet medial line beyond the survey (tube continuation)
    g = np.load(OUT / "survey1970_georef.npz")
    med = g["map_medial_local"]
    cl = list(csv.DictReader(open(OUT / "centreline.csv")))
    cx = np.array([float(r["x"]) for r in cl]); cy = np.array([float(r["y"]) for r in cl])
    cs = np.array([float(r["s"]) for r in cl])
    arc = g["map_medial_arc"]
    s70 = json.loads((OUT / "survey1970.json").read_text())["map_sheet_1970"]["medial"]
    beyond = arc >= s70["our_last_station_nearest_1970_arc_m"]
    ax.plot(med[beyond, 0], med[beyond, 1], color="w", lw=0.8, ls=(0, (2, 2)), alpha=0.9)
    # lidar centreline with 50 m ticks
    ax.plot(cx, cy, color="k", lw=2.0, alpha=0.5)
    ax.plot(cx, cy, color="w", lw=1.1)
    for sk in range(0, int(cs.max()) + 1, 50):
        i = np.argmin(np.abs(cs - sk))
        ax.plot(cx[i], cy[i], "o", ms=2.6, mfc="w", mec="k", mew=0.4)
        if sk > 0:
            ax.annotate(f"{sk}", (cx[i], cy[i]), xytext=(5, -1), textcoords="offset points",
                        fontsize=6, color="k", path_effects=halo(1.4), va="center")
    # skylights
    for j, hh_ in enumerate(HOLES):
        ax.add_patch(Ellipse(hh_["centroid"][:2], 2 * hh_["semi_major"] + 6,
                             2 * hh_["semi_minor"] + 6, angle=np.degrees(hh_["orientation"]),
                             fill=False, ec=C["skylight"], lw=1.0))
        ax.annotate(f"S{j+1}", hh_["centroid"][:2], xytext=(-7, 0), textcoords="offset points",
                    fontsize=7, fontweight="bold", color=C["skylight"], ha="right",
                    va="center", path_effects=halo(1.6))
    # entrance = station s = 0 (the 1970 tie point "entrance line midpoint")
    ax.plot(cx[0], cy[0], marker="v", ms=5, mfc="w", mec="k", mew=0.6)
    ax.annotate("Entrance", (cx[0], cy[0]), xytext=(6, -5), textcoords="offset points",
                fontsize=6.5, path_effects=halo(1.6))
    ax.annotate("Survey end,\n$s$ = %d m" % round(cs.max()), (cx[-1], cy[-1]),
                xytext=(9, 7), textcoords="offset points", fontsize=6, path_effects=halo(1.6))
    ax.text(1244, 720, "Route 39", fontsize=6.5, rotation=-66, ha="center", va="center",
            path_effects=halo(1.6))
    ax.text(1207, 852, "1970 survey", fontsize=6, color="k", rotation=0, ha="center",
            va="center", path_effects=halo(1.6))
    # scale bar and north arrow
    bx, by = xlim[0] + 14, ylim[0] + 14
    ax.plot([bx, bx + 50], [by, by], color="k", lw=2.4)
    ax.plot([bx, bx + 50], [by, by], color="w", lw=1.2)
    ax.text(bx + 25, by + 5, "50 m", ha="center", fontsize=6.5, path_effects=halo(1.6))
    nx_, ny_ = xlim[1] - 18, ylim[0] + 18
    ax.annotate("", (nx_, ny_ + 26), (nx_, ny_),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.0, mutation_scale=8))
    ax.text(nx_, ny_ + 30, "N", ha="center", va="bottom", fontsize=7, fontweight="bold",
            path_effects=halo(1.6))
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.5)

    # Iceland inset (upper-right corner of the map)
    tr3 = Transformer.from_crs(32627, 4326, always_xy=True)
    lon, lat = tr3.transform(cx[0] + ORIGIN[0], cy[0] + ORIGIN[1])
    coast = iceland_coast()
    ins = mm_axes(fig, x0 + w - 27.5, y0 + h - 21.5, 27, 21)
    ins.set_facecolor("w")
    kx = np.cos(np.radians(65.0))
    if coast:
        for line in coast:
            a = np.array(line)
            ins.plot(a[:, 0] * kx, a[:, 1], color="0.35", lw=0.4)
        ins.set_xlim(-24.8 * kx, -13.2 * kx); ins.set_ylim(63.2, 66.7)
        ins.text(-18.6 * kx, 65.0, "Iceland", fontsize=6, ha="center", color="0.25")
    else:
        ins.set_xlim(-24.8 * kx, -13.2 * kx); ins.set_ylim(63.2, 66.7)
        ins.text(-19 * kx, 65.3, "Iceland\n(no offline coastline)", fontsize=5.5, ha="center")
    ins.plot(lon * kx, lat, "o", ms=3.5, mfc=C["skylight"], mec="k", mew=0.4)
    # coordinates are given in the Methods; the inset carries only the marker
    # (a coordinate label overlapped the coastline at page scale)
    ins.set_xticks([]); ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.5)
    print(f"site: entrance at {lat:.4f} N, {abs(lon):.4f} W; survey end s = {cs.max():.0f} m; "
          f"coastline {'NE 50m' if coast else 'none'}")
    return ax, ins


def photo_quadrants():
    """The four frames of mission_tpv_collage.png (a 2 x 2 contact sheet):
    top-left  chamber beneath a skylight with the robot in flight,
    top-right the managed entrance section with its staircase,
    bottom-left  a dark constriction, bottom-right  a breakdown chamber."""
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    col = np.asarray(Image.open(MS_FIGS / "mission_tpv_collage.png").convert("RGB"))
    h, w = col.shape[0] // 2, col.shape[1] // 2
    return dict(robot_chamber=col[:h, :w], entrance=col[:h, w:],
                constriction=col[h:, :w], breakdown=col[h:, w:])


def raster_ax(fig, img, x, y, w, h=None, min_dpi=300.0, max_dpi=450.0):
    """Place a raster at an exact size in millimetres, report its resolution
    at that size and resample anything finer than max_dpi (never upsample)."""
    from PIL import Image
    ph, pw = img.shape[:2]
    if h is None:
        h = w * ph / pw
    dpi = pw / (w * MM)
    if dpi < min_dpi:
        raise ValueError(f"raster {pw}x{ph} is {dpi:.0f} dpi at {w:.0f} mm")
    if dpi > max_dpi:
        f = max_dpi / dpi
        img = np.asarray(Image.fromarray(img).resize(
            (int(round(pw * f)), int(round(ph * f))), Image.LANCZOS))
        dpi = img.shape[1] / (w * MM)
    ax = mm_axes(fig, x, y, w, h)
    ax.imshow(img, interpolation="none", aspect="auto")
    ax.set_axis_off()
    print(f"  raster {w:.1f} x {h:.1f} mm, {dpi:.0f} dpi")
    return ax, h


def fig1():
    """a schematic long section, b the site, c the robot in the tube,
    d two views of the closed section."""
    q = photo_quadrants()
    photos = [("c", q["robot_chamber"]), ("d1", q["constriction"]),
              ("d2", q["breakdown"])]
    wph = 56.6
    hph = wph * photos[0][1].shape[0] / photos[0][1].shape[1]
    fig = plt.figure(figsize=(180 * MM, 137 * MM))
    axp = {}
    for i, (key, img) in enumerate(photos):
        ax, _ = raster_ax(fig, img, 2 + i * (wph + 3.0), 5, wph, hph)
        axp[key] = ax
    axa = mm_axes(fig, 2, 60, 104, 62)
    anatomy_panel(axa)
    axb, _ = site_map(fig, 110, 47.5, 68, 85)
    place_letters(fig, [(axa, "a", None, (0, 0)), (axb, "b", None, (0, 0)),
                        (axp["c"], "c", None, (0, 0)),
                        (axp["d1"], "d", None, (0, 0))])
    save(fig, "fig1_anatomy_site")


# ================================================================== figure 2
def stage_frames():
    """Two stages of the onboard map from exploration_progression.png, a 2 x 2
    contact sheet of the same mission at four times (A top-left to D
    bottom-right); the baked capital letter in each frame is painted out."""
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    prog = np.asarray(Image.open(MS_FIGS / "exploration_progression.png").convert("RGB"))
    h, w = prog.shape[0] // 2, prog.shape[1] // 2
    early, later = prog[:h, :w].copy(), prog[h:, :w].copy()
    for st in (early, later):
        st[:120, :110] = 0
    b = 8                       # the rules of the contact sheet
    early, later = early[b:h - b, b:w - b], later[b:h - b, b:w - b]
    # both frames share one crop, taken from the content of the later one, so
    # the two stages stay in the same frame and can be compared directly
    lit = later.mean(2) > 40
    hh, ww = lit.shape
    rows = np.where(lit.sum(1) > 0.004 * ww)[0]
    cols = np.where(lit.sum(0) > 0.004 * hh)[0]
    pad = 30
    r0, r1 = max(int(rows[0]) - pad, 0), min(int(rows[-1]) + pad, hh - 1)
    c0, c1 = max(int(cols[0]) - pad, 0), min(int(cols[-1]) + pad, ww - 1)
    return (np.ascontiguousarray(early[r0:r1 + 1, c0:c1 + 1]),
            np.ascontiguousarray(later[r0:r1 + 1, c0:c1 + 1]))


def gain_top_view():
    """The top view of vertices_information_gain_instance.png: the executed
    trajectory, the sampled viewpoints and the highest-gain subset. The frame
    below it (the same instant in side view) and the leader that crosses into
    it are dropped, and the frame is cropped to its content."""
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    a = np.asarray(Image.open(MS_FIGS / "vertices_information_gain_instance.png")
                   .convert("RGB")).copy()
    split = 3483                      # the rule between the two views
    a = a[:split]
    a[2900:, 8600:] = 0               # leader pointing down into the side view
    lit = a.mean(2) > 12
    rr, cc = np.where(lit.any(1))[0], np.where(lit.any(0))[0]
    pad = 20
    r0, r1 = max(rr[0] - pad, 0), min(rr[-1] + pad, a.shape[0] - 1)
    c0, c1 = max(cc[0] - pad, 0), min(cc[-1] + pad, a.shape[1] - 1)
    return a[r0:r1 + 1, c0:c1 + 1]


def reconstruction_data():
    """Points of the final onboard-flight reconstruction in the centreline
    principal frame. The 30 cm working copy of the principal flight in the
    registered frame; a few hundred stray returns tens of metres off the
    conduit are dropped for display (no number on the figure comes from them)."""
    from scipy.spatial import cKDTree
    shell = edf.load_xyz(Path(os.environ["EDFIG_SHELL_PCD"]))
    d, _ = cKDTree(np.column_stack([msf.cx, msf.cy])).query(shell[:, :2],
                                                            workers=-1)
    keep = d < 30.0
    print(f"fig2 c: {int((~keep).sum())} of {len(shell):,} points lie more "
          f"than 30 m from the centreline and are not drawn")
    shell = shell[keep]
    u, v = msf.rot(shell[:, 0], shell[:, 1])
    return u, v, shell[:, 2]


def reconstruction_plan(fig, x, y, w, h, data, pad=8.0):
    """Plan view of the final reconstruction, coloured by elevation, with a
    50 m scale bar. Frame: the centreline principal axis, as in the other
    along-tube panels."""
    u, v, z = data
    cu, cv = msf.rot(msf.cx, msf.cy)
    ax = mm_axes(fig, x, y, w, h)
    lo, hi = np.percentile(z, [1.0, 99.0])
    sc = ax.scatter(u, v, c=z, s=0.12, lw=0, cmap="viridis", vmin=lo, vmax=hi,
                    rasterized=True)
    ax.set_aspect("equal")
    half_u = 0.5 * (u.max() - u.min()) + pad
    mid_u = 0.5 * (u.max() + u.min())
    half_v = half_u * h / w
    mid_v = 0.5 * (v.max() + v.min())
    if v.max() - v.min() > 2 * half_v:
        raise ValueError(f"fig2 c: the map spans {v.max() - v.min():.0f} m "
                         f"across, the panel only {2 * half_v:.0f} m")
    ax.set_xlim(mid_u - half_u, mid_u + half_u)
    ax.set_ylim(mid_v - half_v, mid_v + half_v)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    smax = float(msf.cs.max())
    ax.annotate("Entrance, $s$ = 0 m", (cu[0], cv[0]), xytext=(-2, 7),
                textcoords="offset points", fontsize=6.5, ha="left",
                va="bottom", color="0.15", path_effects=halo(1.5))
    ax.annotate(f"End of the mapped conduit, $s$ = {smax:.0f} m",
                (cu[-1], cv[-1]), xytext=(-2, 7), textcoords="offset points",
                fontsize=6.5, ha="right", va="bottom", color="0.15",
                path_effects=halo(1.5))
    x0 = mid_u - half_u + 8.0
    y0 = mid_v - half_v + 6.0
    ax.plot([x0, x0 + 50], [y0, y0], color="k", lw=1.6, solid_capstyle="butt")
    ax.text(x0 + 25, y0 + 2.5, "50 m", ha="center", va="bottom", fontsize=6.5,
            path_effects=halo(1.5))
    cax = mm_axes(fig, x + 4, y + h - 8.0, 32, 2.6)
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label(ELEV_LABEL, fontsize=6, labelpad=2)
    cb.ax.tick_params(labelsize=6, pad=1.5)
    cb.outline.set_linewidth(0.4)
    cb.ax.xaxis.set_ticks_position("top")
    cb.ax.xaxis.set_label_position("top")
    print(f"fig2 c: {len(u):,} points drawn, elevation {z.min():.1f} to "
          f"{z.max():.1f} m (colour scale clipped to the 1st and 99th "
          f"percentiles, {lo:.1f} to {hi:.1f} m), mapped length {smax:.0f} m")
    return ax, cax


def fig2():
    """a map growth, b viewpoint selection, c the final reconstruction,
    d onboard against offline ceiling by 80 m zone."""
    early, later = stage_frames()
    gain = gain_top_view()
    chain = json.loads((ROOT / "analysis_out" /
                        "transform_flf_ed_chain.json").read_text())
    fig = plt.figure(figsize=(180 * MM, 170 * MM))
    # a: two stages of the onboard map
    wst = 87.0
    hst = wst * early.shape[0] / early.shape[1]
    y_a = 170 - 5 - hst
    ax_a1, _ = raster_ax(fig, early, 2, y_a, wst, hst)
    ax_a2, _ = raster_ax(fig, later, 91, y_a, wst, hst)
    for ax, txt in ((ax_a1, "Early in the mission"), (ax_a2, "Later in the mission")):
        ax.text(0, -0.03, txt, transform=ax.transAxes, fontsize=6.5, va="top",
                ha="left", color="0.15")
    # c keeps the true plan-view proportions of the map; d has a fixed
    # height at the foot of the page; b takes whatever is left
    data = reconstruction_data()
    w_c, pad_c = 132.0, 12.0
    h_c = w_c * (data[1].max() - data[1].min()) / \
        (data[0].max() - data[0].min() + 2 * pad_c)
    h_d, y_d = 22.0, 9.0
    y_c = y_d + h_d + 12.0
    y_b_min = y_c + h_c + 5.0
    ratio = gain.shape[0] / gain.shape[1]
    wgn = min(176.0, (y_a - 7.0 - y_b_min) / ratio)
    hgn = wgn * ratio
    y_b = y_a - 7.0 - hgn
    ax_b, _ = raster_ax(fig, gain, 2.0, y_b, wgn, hgn)
    ax_c, cax = reconstruction_plan(fig, 2.0, y_c, w_c, h_c, data, pad=pad_c)
    ax_d = mm_axes(fig, 22, y_d, 140, h_d)
    ned.reg_two_slam(ax_d, chain)
    ax_d.set_xlim(-5, 310)
    ax_d.set_ylabel("Ceiling, onboard minus\noffline map (m)")
    capitalise_labels(fig)
    clamp_fonts(fig)
    place_letters(fig, [(ax_a1, "a", None, (0, 0)), (ax_b, "b", None, (0, 0)),
                        (ax_c, "c", None, (0, 0)), (ax_d, "d", None, (0, 0))])
    save(fig, "fig2_survey")


# ================================================================== figure 3
def showcase_overlay(fig, x, y, w, view="oblique", ticks=(0, 100, 200, 300)):
    """The archived cutaway render with the three skylights ringed and
    along-tube stations marked. Marker positions come from the camera
    recovered by analysis/showcase_camera.py, which is validated against the
    render itself (see analysis_out_v10/showcase_camera.json)."""
    cam_file = OUT10 / "showcase_camera.json"
    if not cam_file.exists():
        raise FileNotFoundError(
            f"{cam_file} missing; run analysis/showcase_camera.py first")
    cams = json.loads(cam_file.read_text())
    cam = cams["views"][view]
    if cam.get("check_cover_fraction", 0.0) < 0.9:
        raise ValueError(f"camera for {view} not validated: {cam}")
    img = plt.imread(OUT / f"renders/showcase_{view}.png")
    nz = (img[..., :3].min(-1) < 0.97)
    rws, cls = np.where(nz.any(1))[0], np.where(nz.any(0))[0]
    r0, c0 = int(rws[0]), int(cls[0])
    img = img[r0:rws[-1] + 1, c0:cls[-1] + 1]

    from importlib import import_module
    sc = import_module("showcase_camera")

    def px(P):
        a, b, d = sc.project(cam, P)
        return a - c0, b - r0, d

    cl = list(csv.DictReader(open(OUT / "centreline.csv")))
    P = np.array([[float(r_["x"]), float(r_["y"]), float(r_["z"])] for r_ in cl])
    cs = np.array([float(r_["s"]) for r_ in cl])
    tx, ty, _ = px(P)
    hx, hy, _ = px(np.array([h_["centroid"] for h_ in HOLES]))
    # trim the terrain beyond the model: keep the tube and a margin of a
    # tenth of the frame on each side
    m = 0.10 * img.shape[1]
    x_lo = max(int(min(tx.min(), hx.min()) - m), 0)
    x_hi = min(int(max(tx.max(), hx.max()) + m), img.shape[1] - 1)
    img = img[:, x_lo:x_hi + 1]
    tx, hx = tx - x_lo, hx - x_lo
    h = w * img.shape[0] / img.shape[1]
    ax = mm_axes(fig, x, y, w, h)
    ax.imshow(img, interpolation="bilinear")
    ax.set_axis_off()
    ax.set_xlim(0, img.shape[1]); ax.set_ylim(img.shape[0], 0)
    print(f"fig3 render: {img.shape[1] / (w / 25.4):.0f} dpi at {w:.0f} mm, "
          f"{h:.1f} mm tall")

    # three skylights, ringed at their photogrammetric rim centroids; the
    # label always sits above the ring so it stays inside the frame
    for j, (a, b) in enumerate(zip(hx, hy)):
        ax.add_patch(Circle((a, b), 95, fill=False, ec=C["skylight"], lw=1.0))
        ax.annotate(f"S{j+1}", (a, b), xytext=(8, 6), textcoords="offset points",
                    fontsize=7, fontweight="bold", color=C["skylight"],
                    ha="left", va="bottom", path_effects=halo(1.6))
    # along-tube stations, so the perspective render carries a distance scale
    for t in ticks:
        i = int(np.argmin(np.abs(cs - t)))
        if not (0 <= tx[i] < img.shape[1] and 0 <= ty[i] < img.shape[0]):
            print(f"fig3 render: station s = {t} m is outside the frame")
            continue
        ax.plot(tx[i], ty[i], "o", ms=2.4, mfc="w", mec="k", mew=0.5)
        ax.annotate(f"$s$ = {t} m", (tx[i], ty[i]), xytext=(5, -5),
                    textcoords="offset points", fontsize=6, color="0.1",
                    ha="left", va="top", path_effects=halo(1.5))
    return ax, h


def repeatability_panel(ax):
    """Flight 1 minus flight 2 ceiling at every one-metre station, from
    analysis/flight2_station_offsets.py; the summary printed on the panel is
    the archived analysis_out_v9/flight2_repeatability.json."""
    rows = list(csv.DictReader(open(OUT10 / "flight2_station_offsets.csv")))
    s = np.array([float(r_["s"]) for r_ in rows])
    dz = np.array([float(r_["dz"]) if r_["dz"] else np.nan for r_ in rows])
    rep = json.loads((OUT / "flight2_repeatability.json").read_text())
    ok = np.isfinite(dz)
    assert int(ok.sum()) == rep["n_stations"], (ok.sum(), rep["n_stations"])
    ax.axhline(0.0, color="0.5", lw=0.6, ls="--")
    ax.axhspan(-rep["rms"], rep["rms"], color=C["envelope"], alpha=0.10, lw=0,
               label=f"±RMS ({rep['rms']:.2f} m)")
    ax.plot(s[ok], dz[ok], lw=0.8, color=C["tube"])
    ax.axhline(rep["median"], color=C["skylight"], lw=0.8,
               label=f"Median {rep['median']:.2f} m")
    ax.set_xlim(-5, 310)
    ax.set_ylim(-1.6, 1.6)
    ax.set_xlabel("Distance along tube $s$ (m)")
    ax.set_ylabel("Ceiling difference,\nflight 1 minus flight 2 (m)")
    ax.legend(loc="upper left", ncol=2, handletextpad=0.5, columnspacing=1.2,
              handlelength=1.6, borderaxespad=0.2)
    print(f"fig3 c: {int(ok.sum())} stations, median {rep['median']:.3f} m, "
          f"RMS {rep['rms']:.3f} m, 5-95% {rep['p5']:.2f} to {rep['p95']:.2f} m")


def fig3():
    """a the skylight constellation, b the registered model, c the
    flight-to-flight repeatability, d a long section of the registered model."""
    rep = json.loads((OUT / "registration_report.json").read_text())
    T = np.array(json.loads((ROOT / "analysis_out" /
                             "transform_ed_slice.json").read_text())["matrix"])
    fig = plt.figure(figsize=(180 * MM, 170 * MM))
    axa = mm_axes(fig, 17, 106, 66, 44)
    ned.reg_constellation(axa, rep, T)
    axa.legend(loc="upper left", fontsize=5.5, ncols=1, handletextpad=0.4,
               borderaxespad=0.2, labelspacing=0.25)
    axb, hb = showcase_overlay(fig, 96, 103, 80)
    if 103 + hb > 166.0:
        raise ValueError(f"fig3 b reaches {103 + hb:.1f} mm, no room for its letter")
    cax = mm_axes(fig, 100, 95, 40, 2.8)
    cb = matplotlib.colorbar.Colorbar(cax, cmap=msf.TAU_CMAP, norm=msf.TAU_NORM,
                                      orientation="horizontal",
                                      ticks=[0, 4, 8, 12, 16])
    cb.set_label("Overburden of the interior in b and d (m)", fontsize=6,
                 labelpad=2)
    cb.ax.tick_params(labelsize=6, pad=1.5)
    cb.ax.xaxis.set_ticks_position("bottom")
    cb.ax.xaxis.set_label_position("top")
    cb.outline.set_linewidth(0.4)
    axc = mm_axes(fig, 20, 64, 154, 28)
    repeatability_panel(axc)
    axd = mm_axes(fig, 20, 8, 154, 44)
    msf.panel_c(axd)
    axd.set_ylabel(ELEV_LABEL.replace(" (", "\n("))
    axd.set_xlabel("Along-tube distance $s$ (m)")
    for t in axd.texts:
        if t.get_text().startswith("vertical exaggeration"):
            t.set_position((1.0, -0.30)); t.set_va("top")
    capitalise_labels(fig)
    clamp_fonts(fig)
    place_letters(fig, [(axa, "a", None, (0, 0)), (axb, "b", None, (0, 0)),
                        (axc, "c", None, (0, 0)), (axd, "d", None, (0, 0))])
    save(fig, "fig3_merged_model")


# ================================================================== figure 4
def fig4():
    rows, s, x, y, tau, span, klass = mf.load_roof()
    sig = float(rows[0]["sigma_tau"])
    assert abs(sig - PN["reg_sigma_tau_m"]) < 1e-6, (sig, PN["reg_sigma_tau_m"])
    it = klass == "intact"
    excl = np.isin(klass, ["multipass", "inconsistent", "low_coverage", "no_dem"])
    sky = klass == "skylight"
    bands = mf.skylight_bands(s, klass)
    mrows = list(csv.DictReader(open(OUT / "morphometry.csv")))
    ok = [r for r in mrows if r["area"]]
    ms_ = np.array([float(r["s"]) for r in ok])
    A = np.array([float(r["area"]) for r in ok])
    W = np.array([float(r["width"]) for r in ok])
    Hh = np.array([float(r["height"]) for r in ok])
    env = REV["envelope"]
    beta, rho = env["beta"], env["rho_mid"]
    print(f"fig4: n intact {it.sum()}, excluded {excl.sum()}, skylight {sky.sum()}, "
          f"sigma band {sig:.2f} m, tau median {np.median(tau[it]):.2f}, "
          f"beta {beta}, rho {rho}, skylight bands {bands}")

    fig = plt.figure(figsize=(180 * MM, 170 * MM))
    xl, wl = 18.0, 100.0          # left column (along-tube panels)
    xr, wr = 134.0, 42.0          # right column
    smax = float(s.max())

    def bands_on(ax, ymax=1.0):
        for lo, hi in bands:
            ax.axvspan(lo, hi, ymax=ymax, color=C["skylight"], alpha=0.15, lw=0)

    # a: plan in the principal-axis frame
    XY = np.column_stack([x, y])
    ctr = XY.mean(0)
    _, _, Vt = np.linalg.svd(XY - ctr, full_matrices=False)
    uv = (XY - ctr) @ Vt.T
    u, v = uv[:, 0], uv[:, 1]
    if u[-1] < u[0]:
        u, v = -u, -v
    axa = mm_axes(fig, xl, 134, wl, 30)
    axa.plot(u, v, color=C["tube"], lw=1.4)
    for lo, hi in bands:
        m = (s >= lo) & (s <= hi)
        axa.plot(u[m], v[m], color=C["skylight"], lw=3.0, solid_capstyle="butt")
    for j, h in enumerate(HOLES):
        i = np.argmin(np.hypot(x - h["centroid"][0], y - h["centroid"][1]))
        axa.annotate(f"S{j+1}", (u[i], v[i]), xytext=(0, -9), textcoords="offset points",
                     fontsize=6.5, fontweight="bold", color=C["skylight"], ha="center")
    for tick in range(0, int(smax) + 1, 50):
        j = np.argmin(np.abs(s - tick))
        axa.plot(u[j], v[j], "o", ms=2.2, mfc="w", mec=C["tube"], mew=0.6)
        axa.annotate(f"{tick}", (u[j], v[j]), fontsize=6, color="0.3",
                     xytext=(0, 5), textcoords="offset points", ha="center")
    axa.set_aspect("equal")
    axa.set_xlim(u.min() - 8, u.max() + 8)
    axa.set_ylim(v.min() - 8, v.max() + 14)
    axa.set_xlabel("Distance along the principal axis (m)")
    axa.set_ylabel("Across (m)")
    axa.text(u[0] - 1, v[0] - 4.5, "Entrance", fontsize=6, ha="left", va="top", color="0.3")

    # b: area, then width and height, two strips
    axb1 = mm_axes(fig, xl, 100, wl, 23)
    bands_on(axb1)
    axb1.plot(ms_, A, lw=0.8, color=C["tube"])
    axb1.set_ylabel("Area (m$^2$)")
    axb1.set_xlim(0, smax); axb1.tick_params(labelbottom=False)
    axb2 = mm_axes(fig, xl, 72, wl, 23)
    bands_on(axb2)
    axb2.plot(ms_, W, lw=0.8, color=C["tube"])
    axb2.plot(ms_, Hh, lw=0.8, color=C["intact"])
    axb2.set_ylabel("Extent (m)")
    axb2.set_xlim(0, smax); axb2.tick_params(labelbottom=False)
    j = np.argmax(W)
    axb2.annotate("width", (ms_[j], W[j]), xytext=(4, 0), textcoords="offset points",
                  fontsize=6, color=C["tube"], va="center")
    jh = np.argmin(np.abs(ms_ - 200))
    axb2.annotate("height", (ms_[jh], Hh[jh]), xytext=(0, -10), textcoords="offset points",
                  fontsize=6, color=C["intact"], ha="center")

    # c: overburden profile with the systematic band
    axc = mm_axes(fig, xl, 14, wl, 46)
    bands_on(axc, ymax=0.84)                # key sits in the top strip, off the bands
    # band on the station grid, broken where stations are not intact
    tau_g = np.full(len(s), np.nan)
    tau_g[it] = tau[it]
    axc.fill_between(s, tau_g - sig, tau_g + sig, color=C["intact"], alpha=0.18, lw=0,
                     label=f"systematic band ±{sig:.1f} m")
    axc.plot(s[it], tau[it], "o", ms=1.8, color=C["intact"], label="intact roof")
    if excl.any():
        axc.plot(s[excl], np.full(excl.sum(), -0.8), "|", ms=5, color=C["flagged"],
                 label="excluded by the consistency screen")
    axc.axhline(0, color="k", lw=0.5)
    axc.set_xlim(0, smax)
    axc.set_ylim(-1.6, 19.5)
    axc.set_yticks([0, 5, 10, 15])
    axc.set_ylabel("Overburden $\\tau$ (m)")
    axc.set_xlabel("Along-tube distance $s$ (m)")
    axc.legend(loc="upper right", ncol=3, handletextpad=0.5, borderaxespad=0.2,
               columnspacing=1.2, handlelength=1.5)
    for ax in (axb1, axb2, axc):
        ax.spines["left"].set_position(("outward", 2))

    # d: distribution with the shielding-mass axis on top
    axd = mm_axes(fig, xr, 108, wr, 42)
    axd.hist(tau[it], bins=16, color=C["intact"], alpha=0.9)
    med = np.median(tau[it])
    axd.axvline(med, color="k", lw=0.8)
    axd.set_xlabel("Overburden $\\tau$ (m)")
    axd.set_ylabel("Stations")
    axd.spines["top"].set_visible(True)
    ax2 = axd.twiny()
    ax2.set_xlim(np.array(axd.get_xlim()) * rho / 10.0)     # m -> g cm^-2 at rho
    ax2.set_xlabel("Shielding mass (g cm$^{-2}$)")
    ax2.set_xticks([0, 1000, 2000, 3000])
    comma_fmt(ax2, "x")
    ax2.tick_params(labelsize=6)
    ax2.spines["right"].set_visible(False)
    print(f"fig4 d: shielding axis at rho {rho:.0f} kg m^-3 -> {rho / 10:.0f} g cm^-2 per m; "
          f"median tau {med:.2f} m -> {med * rho / 10:.0f} g cm^-2 "
          f"(paper_numbers shielding.median {PN['shielding']['median']})")

    # e: overburden against span with the clamped-strip curves
    axe = mm_axes(fig, xr, 14, wr, 78)
    axe.scatter(span[it], tau[it], s=7, color=C["intact"], lw=0, label="intact roof")
    axe.scatter(span[sky], np.zeros(sky.sum()), s=11, marker="v", color=C["skylight"],
                lw=0, label="skylight (failed)")
    Ls = np.linspace(3, 33, 120)
    ymax = 16.5
    for st_mpa, ls, name in ((0.5, "--", "$\\sigma_t$ = 0.5 MPa"), (1.0, "-.", "1 MPa"),
                             (2.5, ":", "2.5 MPa")):
        yy = beta * rho * 9.81 * Ls ** 2 / (st_mpa * 1e6)
        axe.plot(Ls, yy, "k", ls=ls, lw=0.7)
        m_in = (yy <= ymax - 0.3) & (Ls <= 32.5)
        xe, ye = Ls[m_in][-1], yy[m_in][-1]
        if ye >= ymax - 1.5:            # leaves through the top: label to its left
            axe.annotate(name, (xe, ye), xytext=(-3, -2), textcoords="offset points",
                         fontsize=6, ha="right", va="top", color="0.2")
        elif st_mpa < 2:                # leaves through the right edge: label above
            axe.annotate(name, (xe, ye), xytext=(0, 3), textcoords="offset points",
                         fontsize=6, ha="right", va="bottom", color="0.2")
        else:                           # flattest curve: label below, clear of the points
            axe.annotate(name, (xe, ye), xytext=(-2, -4), textcoords="offset points",
                         fontsize=6, ha="right", va="top", color="0.2")
    axe.set_xlim(0, 33)
    axe.set_ylim(-0.8, ymax)
    axe.set_xlabel("Local span $L$ (m)")
    axe.set_ylabel("Overburden $\\tau$ (m)")
    axe.legend(loc="lower right", bbox_to_anchor=(1.03, 0.03), handletextpad=0.3,
               handlelength=1.0, borderaxespad=0.0, labelspacing=0.3)
    print("fig4 e: curve check at L = 20 m, 1 MPa: tau_min = "
          f"{beta * rho * 9.81 * 20 ** 2 / 1e6:.2f} m")

    capitalise_labels(fig)
    # reading order: left column a-c (plan, geometry, overburden),
    # right column d, e (distribution, span)
    place_letters(fig, [(axa, "a", "L", (0, 0)), ([axb1, axb2], "b", "L", (0, 0)),
                        (axc, "c", "L", (0, 0)),
                        ([axd, ax2], "d", "R", (0, 0)), (axe, "e", "R", (0, 0))])
    save(fig, "fig4_geometry_cover")


# ================================================================== figure 5
def fig5():
    """fig_planetary from make_figures, re-laid at 180 x 118 mm."""
    real_figure, real_close = mf.plt.figure, mf.plt.close
    holder = {}

    def sized_figure(*a, **k):
        k["figsize"] = (180 * MM, 118 * MM)
        holder["fig"] = real_figure(*a, **k)
        return holder["fig"]

    mf.plt.figure, mf.plt.close = sized_figure, lambda *a, **k: None
    try:
        mf.fig_planetary()
    finally:
        mf.plt.figure, mf.plt.close = real_figure, real_close
    fig = holder["fig"]
    assert round(PN["reg_sigma_tau_m"], 1) == 1.5, PN["reg_sigma_tau_m"]   # dashed line label
    axes = fig.axes
    letters = []
    for ax in axes:
        t = ax.get_title(loc="left")
        if t[:1] in "abc" and t[1:3] == "  ":
            letters.append((ax, t[0]))
        ax.set_title("", loc="left")
    axa, axb, axc = [a for a, _ in letters]
    axb.set_ylabel("Maximum stable span $L_{\\max}$ (m)")
    axa.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else ""))
    axc.set_xlabel("Per-station overburden error, RMS (m)")
    fig.subplots_adjust(left=0.16, right=0.985, top=0.95, bottom=0.10, wspace=0.30, hspace=0.45)
    # the key of a hid the top of the impact-melt curve: move it above the axes
    pos = axa.get_position()
    dh = 11.0 / 118.0
    axa.set_position([pos.x0, pos.y0, pos.width, pos.height - dh])
    axa.legend(loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False, fontsize=6,
               handlelength=1.6, handletextpad=0.4, borderaxespad=0.0, labelspacing=0.25)
    capitalise_labels(fig)
    clamp_fonts(fig)
    place_letters(fig, [(axa, "a", "T", (0, 0)), (axb, "b", None, (0, 0)), (axc, "c", "T", (0, 0))])
    save(fig, "fig5_planetary")


# ================================================================== QA
# Machado, Oliveira & Fernandes (2009) deuteranopia matrix, severity 1.0, on
# linear RGB (the "3x3 matrix" option of the figure brief).
DEUTER = np.array([[0.367322, 0.860646, -0.227968],
                   [0.280085, 0.672501, 0.047413],
                   [-0.011820, 0.042940, 0.968881]])


def simulate_deuteranopia(png_in: Path, png_out: Path) -> None:
    from PIL import Image
    a = np.asarray(Image.open(png_in).convert("RGB")).astype(float) / 255.0
    lin = np.where(a <= 0.04045, a / 12.92, ((a + 0.055) / 1.055) ** 2.4)
    sim = np.clip(lin @ DEUTER.T, 0, 1)
    srgb = np.where(sim <= 0.0031308, sim * 12.92, 1.055 * sim ** (1 / 2.4) - 0.055)
    Image.fromarray((srgb * 255).round().astype(np.uint8)).save(png_out)


def qa() -> None:
    """Rasterise every figure at 110 and 300 dpi, simulate deuteranopia, check
    fonts and page size; the images are then inspected by eye (rule 0)."""
    import subprocess
    qdir = OUT10 / "figure_qa_nature"
    qdir.mkdir(exist_ok=True)
    for pdf in sorted(FIGS.glob("fig?_*.pdf")):
        info = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True).stdout
        size = [l for l in info.splitlines() if l.startswith("Page size")][0]
        pts = [float(v) for v in size.split(":")[1].split("pts")[0].split("x")]
        wmm, hmm = pts[0] * 25.4 / 72, pts[1] * 25.4 / 72
        fonts = subprocess.run(["pdffonts", str(pdf)], capture_output=True, text=True).stdout
        names = sorted({l.split()[0].split("+")[-1] for l in fonts.splitlines()[2:] if l.strip()})
        ok_size = abs(wmm - 180) < 0.2 and hmm <= 170.2
        ok_font = all(n.startswith("NimbusSans") for n in names) and names
        print(f"{pdf.name}: {wmm:.1f} x {hmm:.1f} mm {'OK' if ok_size else 'FAIL'}; "
              f"fonts {names} {'OK' if ok_font else 'FAIL'}")
        for dpi in (110, 300):
            stem = qdir / f"{pdf.stem}_{dpi}"
            subprocess.run(["pdftoppm", "-r", str(dpi), "-png", "-singlefile", str(pdf), str(stem)],
                           check=True)
            if dpi == 110:
                simulate_deuteranopia(stem.with_suffix(".png"), qdir / f"{pdf.stem}_deut.png")
    print("QA rasters in", qdir)


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3", "4", "5", "qa"]
    for k in which:
        {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "qa": qa}[k]()
