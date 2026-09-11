"""Extended Data items for the Nature initial submission (plan/13, Step 4).

Ten items, each a PDF plus a 300 dpi PNG in analysis_out_v10/figures_nature/:

    ed1_site_photos              skylights, interior, FPV frames (rasters)
    ed2_platform_architecture    platform photo, behaviour tree (vector redraw)
    ed3_planner_behaviour        planner screenshots (rasters)
    ed4_deployments              placeholder only (--only ed4); the real item
                                 comes from make_nature_ed4_deployments.py
    ed5_registration_consistency registration validation + consistency screen
    ed6_centreline_morphometry   centreline, sections, sinuosity, aspect ratio
    ed7_dem_validation           ArcticDEM check + orbital-quality profile
    ed8_survey1970               the 1970 survey and the road crossing
    ed_table1_missions           mission summary of the two long flights
    ed_table2_planner_params     platform configuration + planner parameters

Nature style for Extended Data (not restyled by the journal): 180 mm wide,
<= 170 mm tall, Nimbus Sans (Helvetica clone) throughout, 7 pt text, 6 pt
ticks, bold 8 pt panel letters, no titles inside panels (the legend carries
the description), RGB. Tables are rendered as images with booktabs rules.

Data: every number comes from analysis_out_v9/ (paper run), the legacy
analysis_out/ transforms that the v9 registration panel already compares
against, analysis_out_ed/mission_evidence.json (mission table) and
analysis_out/planner_params.json (planner table). Plotting code follows the
v9 panel functions in make_ed_figures.py, make_dem_validation_figure.py,
survey1970_compare.figure() and make_figures.fig_morphometry(); the loaders
and geometry helpers are imported from those modules, not copied.

Run (repo root):
    env -u PYTHONPATH PYTHONPATH=src .venv/bin/python \
        analysis/make_nature_ed_figures.py [--only ed5,ed_table1]
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# The v9 (orthometric) run is the paper run; the imported modules read their
# inputs from these variables at import time, so they are fixed before import.
os.environ.setdefault("ANALYSIS_OUT", str(ROOT / "analysis_out_v9"))
os.environ.setdefault("ROOF_DEM_PCD", "maps/aerial_isn16_ortho_dem.pcd")
os.environ.setdefault("EDFIG_SHELL_PCD", str(ROOT / "maps/flf_30cm_ed_ortho.pcd"))
os.environ.setdefault("EDFIG_SECT_PCD", str(ROOT / "maps/flf_10cm_ed_ortho.pcd"))
os.environ.setdefault("EDFIG_AERIAL_PCD", str(ROOT / "maps/aerial_isn16_ortho_10cm.pcd"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402
from matplotlib.textpath import TextToPath  # noqa: E402
from PIL import Image  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import behaviour_tree as bt  # noqa: E402
import make_ed_figures as edf  # noqa: E402
import orbital_dem_test as odt  # noqa: E402
import survey1970_compare as s70  # noqa: E402
from make_figures import skylight_bands  # noqa: E402

OUT = Path(os.environ["ANALYSIS_OUT"])
OUT_LEGACY = ROOT / "analysis_out"          # 4-DOF baseline + Ed transforms
OUT_ED = ROOT / "analysis_out_ed"           # mission_evidence.json
FIGS = ROOT / "analysis_out_v10" / "figures_nature"
RASTERS = Path(os.environ.get("ED_RASTER_DIR", str(
    ROOT.parent / "_Nature__Autonomous_aerial_reconnaissance_of_a_basaltic_"
    "lava_tube_reveals_interior_morphology_and_roof_thickness_distribution_"
    "for_planetary_subsurface" / "figures")))
Image.MAX_IMAGE_PIXELS = None

MM = 1.0 / 25.4
W_MM = 180.0
ELEV = "Elevation (m a.s.l., ISH2004)"
if edf.ELEV_LABEL != "elevation (m a.s.l., ISH2004)":
    raise RuntimeError(f"unexpected datum label {edf.ELEV_LABEL!r}; "
                       "this script is written for the v9 orthometric run")
C = edf.C
CLASS_COLOUR = edf.CLASS_COLOUR
FONT = "Nimbus Sans"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [FONT],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0,
    "xtick.major.pad": 2.0, "ytick.major.pad": 2.0,
    "axes.labelpad": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.handlelength": 1.5,
    "legend.borderaxespad": 0.3, "legend.handletextpad": 0.5,
    "lines.linewidth": 0.9,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": FONT,
    "mathtext.it": f"{FONT}:italic", "mathtext.bf": f"{FONT}:bold",
    "mathtext.default": "it", "mathtext.cal": FONT,
    "figure.dpi": 100, "savefig.dpi": 300,
})
LETTER_KW = dict(fontsize=8, fontweight="bold", ha="left", va="bottom")


# ----------------------------------------------------------------------
# generic helpers
# ----------------------------------------------------------------------
def new_fig(h_mm: float):
    if h_mm > 170.0 + 1e-6:
        raise ValueError(f"figure height {h_mm} mm exceeds 170 mm")
    return plt.figure(figsize=(W_MM * MM, h_mm * MM))


def ax_mm(fig, x, y, w, h):
    """Axes at (x, y) mm from the top-left corner, w x h mm."""
    fw, fh = fig.get_size_inches() / MM
    return fig.add_axes([x / fw, 1 - (y + h) / fh, w / fw, h / fh])


def load_png(name: str) -> np.ndarray:
    im = Image.open(RASTERS / name)
    if im.mode == "RGBA":            # flatten onto white; the PDF is RGB
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        im = bg
    return np.asarray(im.convert("RGB"))


def image_panel(fig, name, x, y, w, h=None, letter=None, min_dpi=300.0,
                max_dpi=600.0):
    """Place a raster at exact size (mm) and check its resolution. Rasters
    above max_dpi at the final size are resampled (Lanczos) to keep the PDF
    within upload limits; nothing is ever upsampled."""
    img = load_png(name)
    ph, pw = img.shape[:2]
    if h is None:
        h = w * ph / pw
    dpi = pw / (w * MM)
    if dpi < min_dpi:
        raise ValueError(f"{name}: {dpi:.0f} dpi at {w:.0f} mm (< {min_dpi})")
    if dpi > max_dpi:
        f = max_dpi / dpi
        img = np.asarray(Image.fromarray(img).resize(
            (int(round(pw * f)), int(round(ph * f))), Image.LANCZOS))
        dpi = img.shape[1] / (w * MM)
    ax = ax_mm(fig, x, y, w, h)
    ax.imshow(img, interpolation="none", aspect="auto")
    ax.set_axis_off()
    if letter:
        # letters sit 0.8 mm above the raster, on white, never on the image
        ax.text(0.0, 1.0 + 0.8 / h, letter, transform=ax.transAxes,
                clip_on=False, **LETTER_KW)
    print(f"    {name}: {w:.1f} x {h:.1f} mm, {dpi:.0f} dpi")
    return ax, h


def place_letters(fig, items, dy_mm=0.8):
    """Bold panel letters at the top-left of each panel's tight bbox
    (labels included), 0.8 mm above it. Call after layout is final."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fw, fh = fig.get_size_inches()
    for ax, ch in items:
        axes = ax if isinstance(ax, (list, tuple)) else [ax]
        bb = None
        for a in axes:
            b = a.get_tightbbox(r).transformed(fig.transFigure.inverted())
            bb = b if bb is None else matplotlib.transforms.Bbox.union([bb, b])
        fig.text(bb.x0, bb.y1 + dy_mm * MM / fh, ch, **LETTER_KW)


def freeze_layout(fig):
    fig.canvas.draw()
    fig.set_layout_engine("none")


def save(fig, stem: str):
    FIGS.mkdir(parents=True, exist_ok=True)
    pdf, png = FIGS / f"{stem}.pdf", FIGS / f"{stem}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300, facecolor="white")
    Image.open(png).convert("RGB").save(png, dpi=(300, 300))   # RGB, no alpha
    plt.close(fig)
    w, h = fig.get_size_inches() / MM
    print(f"  {stem}: {w:.0f} x {h:.0f} mm -> {pdf.name}, {png.name}")


def layout(fig, w_pad=1.0, h_pad=4.5, top_mm=4.0, hspace=0.0, wspace=0.02):
    fh = fig.get_size_inches()[1] / MM
    fig.set_layout_engine("constrained")
    fig.get_layout_engine().set(w_pad=w_pad * MM, h_pad=h_pad * MM,
                                hspace=hspace, wspace=wspace,
                                rect=(0, 0, 1, 1 - top_mm / fh))


# ----------------------------------------------------------------------
# table renderer (booktabs look)
# ----------------------------------------------------------------------
_T2P = TextToPath()


def text_width_pt(s: str, size=7.0, bold=False) -> float:
    fp = FontProperties(family=FONT, size=size,
                        weight="bold" if bold else "normal")
    ismath = ("$" in s)
    w, _, _ = _T2P.get_text_width_height_descent(s, fp, ismath=ismath)
    return w


def wrap(s: str, width_pt: float, size=7.0) -> list[str]:
    words, lines, cur = s.split(" "), [], ""
    for w in words:
        t = (cur + " " + w).strip()
        if cur and text_width_pt(t, size) > width_pt:
            lines.append(cur); cur = w
        else:
            cur = t
    if cur:
        lines.append(cur)
    return lines


def draw_table(fig, x, y, w, header, rows, col_w, align=None,
               row_h=3.0, size=7.0, letter=None, section_rows=()):
    """Booktabs table at (x, y) mm, w mm wide. rows: list of cell lists;
    a row whose index is in section_rows is an italic full-width heading.
    Returns the table height in mm."""
    col_w = np.asarray(col_w, float) / np.sum(col_w) * w
    align = align or ["left"] * len(header)
    pad = 1.0                                       # mm inside a cell
    # wrap every cell first so the height is known
    wrapped = []
    for i, row in enumerate(rows):
        if i in section_rows:
            wrapped.append([[row[0]]] + [[""]] * (len(header) - 1)); continue
        wrapped.append([wrap(c, (cw - 2 * pad) * MM * 72, size)
                        for c, cw in zip(row, col_w)])
    n_lines = [max(len(c) for c in r) for r in wrapped]
    h = row_h * (1 + sum(n_lines)) + 2.0            # header + rows + rules
    ax = ax_mm(fig, x, y, w, h)
    ax.set_xlim(0, w); ax.set_ylim(h, 0); ax.set_axis_off()
    xs = np.concatenate([[0], np.cumsum(col_w)])

    def cell(cx, cy, txt, al, **kw):
        if al == "left":
            ax.text(cx + pad, cy, txt, ha="left", va="center", fontsize=size, **kw)
        elif al == "right":
            ax.text(cx - pad, cy, txt, ha="right", va="center", fontsize=size, **kw)
        else:
            ax.text(cx, cy, txt, ha="center", va="center", fontsize=size, **kw)

    def rule(yy, lw):
        ax.plot([0, w], [yy, yy], color="k", lw=lw, solid_capstyle="butt")

    yy = 0.3
    rule(yy, 0.75)
    yc = yy + row_h / 2 + 0.2
    for j, hd in enumerate(header):
        cx = xs[j] if align[j] == "left" else (xs[j + 1] if align[j] == "right"
                                                else (xs[j] + xs[j + 1]) / 2)
        cell(cx, yc, hd, align[j], fontweight="bold")
    yy += row_h + 0.4
    rule(yy, 0.5)
    for i, r in enumerate(wrapped):
        if i in section_rows:
            rule(yy, 0.5) if i > 0 else None
            cell(xs[0], yy + row_h / 2 + 0.2, r[0][0], "left", style="italic")
            yy += row_h
            continue
        for j, lines in enumerate(r):
            cx = xs[j] if align[j] == "left" else (
                xs[j + 1] if align[j] == "right" else (xs[j] + xs[j + 1]) / 2)
            for k, ln in enumerate(lines):
                cell(cx, yy + row_h * (k + 0.5) + 0.2, ln, align[j])
        yy += row_h * n_lines[i]
    yy += 0.4
    rule(yy, 0.75)
    if letter:
        ax.text(0.0, 1.0 + 0.8 / h, letter, transform=ax.transAxes,
                clip_on=False, **LETTER_KW)
    return h


# ----------------------------------------------------------------------
# ED 1: site photographs (rasters)
# ----------------------------------------------------------------------
def ed1_site_photos():
    fig = new_fig(170.0)
    # each raster at most as wide as 300 dpi allows (skylights.png 1740 px,
    # fpv_combo.png 1538 px), left-aligned, 4 mm gaps, letters above
    y = 4.0
    _, h = image_panel(fig, "skylights.png", 0, y, 147.0, letter="a")
    y += h + 5.0
    _, h = image_panel(fig, "lava_tube_merged.png", 0, y, 180.0, letter="b")
    y += h + 5.0
    _, h = image_panel(fig, "fpv_combo.png", 0, y, 130.0, letter="c")
    y += h
    fig.set_size_inches(W_MM * MM, y * MM)   # trim to content
    _refit_axes(fig, 170.0, y)
    save(fig, "ed1_site_photos")


def _refit_axes(fig, old_h, new_h):
    """Axes were placed against the old page height; rescale to the new."""
    for ax in fig.axes:
        l, b, w, h = ax.get_position().bounds
        top = (1 - b - h) * old_h
        ax.set_position([l, 1 - (top + h * old_h) / new_h, w, h * old_h / new_h])


# ----------------------------------------------------------------------
# ED 2: platform photo, behaviour tree (vector); platform table in ED Table 2
# ----------------------------------------------------------------------
def platform_rows():
    """ED Table 1 rows. Numeric cells that exist in a JSON are read from
    it; the descriptive cells and the surface-survey figures are transcribed
    from the manuscript (and were checked against the Metashape report)."""
    pp = json.load(open(OUT_LEGACY / "planner_params.json"))
    datum = json.load(open(OUT / "datum_transfer.json"))
    dz = abs(datum["nearest_neighbour"]["dz_new_minus_old"]["median"])
    vox = pp["mapping"]["tsdf_voxel_size_m"]
    ray = pp["mapping"]["max_ray_length_m"]
    rows = [
        ["Aerial platform and autonomy stack"],
        ["Airframe", "in-house quadrotor, about 0.8 m rotor span, roll cage"],
        ["Ranging sensor", "Ouster spinning 3D lidar (360° azimuth)"],
        ["Inertial sensing", "VectorNav IMU"],
        ["Onboard computer", "Intel NUC-class single-board computer (full autonomy stack)"],
        ["Flight controller", "Pixhawk Cube Orange (attitude stabilisation)"],
        ["Power", "lithium high-voltage (LiHV) battery"],
        ["Context camera", "first-person-view camera, separately powered, not used by autonomy"],
        ["State estimation", "DLIO lidar-inertial odometry, onboard, no external fix"],
        ["Mapping", f"hashed TSDF/ESDF volumetric map, {vox:g} m voxels, {ray:g} m ray integration"],
        ["Planning", "STAGE graph-based exploration (Supplementary Methods)"],
        ["Supervision", "behaviour tree (Extended Data Fig. 2b), onboard"],
        ["Surface photogrammetric survey"],
        ["Aircraft, camera", "senseFly eBee X RTK fixed wing, S.O.D.A. camera (10.6 mm, 5472 × 3648 px)"],
        ["Acquisition", "28 August 2025, 135 m above ground, 964 images, 2.08 km$^2$, no snow"],
        ["Processing", "Agisoft Metashape Professional 2.1.3; 7.43 million tie points, 0.48 px reprojection error"],
        ["Products", "dense cloud 6.99 × 10$^8$ points (3.08 cm ground sampling distance); DEM at 15.6 points m$^{\\mathrm{-2}}$; orthomosaic"],
        ["Georeferencing", "onboard GNSS with IceCORS real-time correction, camera stations to 9.5 cm; three ground control points (Trimble R8s), 14.0 cm total RMSE, 7.6 cm in height"],
        ["Reference frame", f"ISN2016 / Lambert 2016, orthometric heights (ISH2004 geoid); delivered ellipsoidal heights differ by {dz:.2f} m at the site"],
    ]
    return rows, (0, 12)


def ed2_platform_architecture():
    """a, the platform photograph (85 mm, so its baked labels print at about
    5 pt); the node key beside it. b, the behaviour tree redrawn as vectors
    from the Groot screenshot by behaviour_tree.py, full width."""
    fig = new_fig(170.0)
    y0 = 4.0
    _, ha = image_panel(fig, "shafterx2.png", 0, y0, 85.0, letter="a", max_dpi=350.0)
    st = bt.Style()
    axk = ax_mm(fig, 92.0, y0 + 2.0, 88.0, 16.0)
    axk.set_xlim(0, 88.0); axk.set_ylim(16.0, 0); axk.set_aspect("equal"); axk.set_axis_off()
    axk.text(0.0, 1.0, "Node types in b", fontsize=7, fontweight="bold", va="center")
    bt.draw_key(axk, 0.0, 5.0, st)
    y = y0 + ha + 5.0
    probe = bt.build(bt.TREE, st)
    tw, th, _ = bt.layout(probe, st)
    if tw > 180.0:
        raise ValueError(f"behaviour tree {tw:.1f} mm wide (> 180)")
    hb = th + 1.0
    axt = ax_mm(fig, 0, y, 180.0, hb)
    bt.render(axt, 180.0, hb, st, key=False)
    axt.text(0.0, 1.0 + 0.8 / hb, "b", transform=axt.transAxes, clip_on=False, **LETTER_KW)
    print(f"    behaviour tree: {tw:.1f} x {th:.1f} mm (vector)")
    y += hb
    fig.set_size_inches(W_MM * MM, y * MM)
    _refit_axes(fig, 170.0, y)
    save(fig, "ed2_platform_architecture")


# ----------------------------------------------------------------------
# ED 3: planner behaviour (rasters)
# ----------------------------------------------------------------------
def ed3_planner_behaviour():
    fig = new_fig(170.0)
    y = 4.0
    _, h = image_panel(fig, "vertices_information_gain_instance.png", 0, y,
                       180.0, letter="a")
    y += h + 6.0
    ib = load_png("inaccessible_frontier_instance.png").shape
    ic = load_png("global_graph_views.png").shape
    rb, rc_ = ib[1] / ib[0], ic[1] / ic[0]
    hrow = (180.0 - 4.0) / (rb + rc_)
    _, _ = image_panel(fig, "inaccessible_frontier_instance.png", 0, y,
                       rb * hrow, letter="b")
    _, _ = image_panel(fig, "global_graph_views.png", rb * hrow + 4.0, y,
                       rc_ * hrow, letter="c")
    y += hrow
    fig.set_size_inches(W_MM * MM, y * MM)
    _refit_axes(fig, 170.0, y)
    save(fig, "ed3_planner_behaviour")


# ----------------------------------------------------------------------
# ED 4: placeholder
# ----------------------------------------------------------------------
def ed4_deployments():
    fig = new_fig(60.0)
    ax = ax_mm(fig, 0, 0, 180.0, 60.0)
    ax.set_axis_off()
    ax.add_patch(plt.Rectangle((0.002, 0.005), 0.996, 0.99, fill=False,
                               ec="0.6", lw=0.5, ls="--",
                               transform=ax.transAxes))
    ax.text(0.5, 0.5, "Extended Data Fig. 4: image to be supplied\n"
            "(deployments in two underground mines and a road tunnel)",
            ha="center", va="center", fontsize=9, color="0.5",
            transform=ax.transAxes)
    save(fig, "ed4_deployments")


def _tick_offset(tick, s_max, dy):
    """Station labels: first to the right, last to the right of the end
    (the centreline dives there and an overhead label crosses the curve),
    others above."""
    if tick == 0:
        return (2, dy)
    if tick + 50 > s_max:
        return (5, -2)
    return (0, dy)


def _tick_ha(tick, s_max):
    return "left" if (tick == 0 or tick + 50 > s_max) else "center"


# ----------------------------------------------------------------------
# ED 5: registration validation and consistency screen
# ----------------------------------------------------------------------
def reg_constellation(ax, rep, T):
    matched_aer = {p[1] for p in rep["landmark"]["inliers"]}
    cluster = []
    for h in rep["aerial_holes"]:
        cx, cy, _ = h["centroid"]
        if h["id"] not in matched_aer:
            continue
        cluster.append((cy, cx))
        ang = 90.0 - np.degrees(h.get("orientation", 0.0))
        ax.add_patch(Ellipse((cy, cx), 2 * max(h["semi_major"], 0.8),
                             2 * max(h["semi_minor"], 0.8), angle=ang,
                             fill=False, lw=1.0, edgecolor=C["surface"]))
        ax.annotate(f"S{sorted(matched_aer).index(h['id']) + 1}", (cy, cx),
                    xytext=(9, 0), textcoords="offset points", ha="left",
                    va="center", fontsize=7)
    th = np.array([h["centroid"] for h in rep["tube_holes"]])
    th_a = (T @ np.column_stack([th, np.ones(len(th))]).T).T[:, :3]
    ax.scatter(th_a[:, 1], th_a[:, 0], marker="+", s=40, lw=1.0,
               color=C["tube"], label="Interior rim centroid (registered)",
               zorder=5)
    ax.scatter([], [], marker="o", facecolor="none", edgecolor=C["surface"],
               s=40, label="Surface-model opening")
    ax.plot([], [], "x", color=C["flagged"], ms=4, ls="none",
            label="Rejected 4th detection (inset)")
    cl = np.array(cluster)
    ax.set_xlim(cl[:, 0].min() - 8, cl[:, 0].max() + 8)
    ax.set_ylim(cl[:, 1].min() - 10, cl[:, 1].max() + 15)
    ax.set_aspect("equal")
    ax.set_xlabel("Northing (m, local frame)")
    ax.set_ylabel("Easting (m)")
    ax.legend(loc="upper left", ncols=1)
    axi = ax.inset_axes([0.70, 0.05, 0.28, 0.30])
    rej = [h for h in rep["aerial_holes"] if h["id"] not in matched_aer]
    axi.scatter(cl[:, 0], cl[:, 1], s=6, color=C["surface"], lw=0)
    for h in rej:
        axi.plot([h["centroid"][1]], [h["centroid"][0]], "x", ms=4,
                 color=C["flagged"])
    axi.set_aspect("equal"); axi.margins(0.45)
    axi.set_xticks([]); axi.set_yticks([])
    for sp in axi.spines.values():
        sp.set_visible(True); sp.set_color("0.6"); sp.set_linewidth(0.5)


def reg_floor_hist(ax, val):
    """Floor seen through S3 against the interior floor (recomputed from the
    point clouds exactly as make_ed_figures.fig_registration does)."""
    from scipy.spatial import cKDTree
    from lava_pcd.io.pcd_reader import BinaryPcdReader
    aer_holes = json.load(open(OUT_LEGACY / "aerial_holes.json"))
    s3 = np.array(aer_holes["holes"][2]["centroid"])

    def load_near(p):
        parts = []
        with BinaryPcdReader(p) as r:
            for c in r.chunks():
                q = c[:, :3].astype(np.float64)
                parts.append(q[np.hypot(q[:, 0] - s3[0], q[:, 1] - s3[1]) < 5.0])
        return np.vstack(parts)

    aer = load_near(edf.AERIAL_PCD)
    tub = load_near(edf.SECT_PCD)
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
    ax.axvline(0, color="k", lw=0.6)
    ax.axvline(st["median_m"], color=C["skylight"], lw=0.9, ls="--")
    ax.annotate(f"median {st['median_m']:.2f} m\nRMS {st['rms_m']:.2f} m",
                (0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
                fontsize=6, color="0.25")
    ax.set_xlabel("Floor offset at S3 (m)")
    ax.set_ylabel("Points")


def reg_ceiling_two_solutions(ax):
    from scipy.spatial import cKDTree
    base_rows = list(csv.DictReader(open(OUT_LEGACY / "roof_thickness.csv")))
    ed_rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    g = lambda rows, k: np.array([float(r[k]) if r[k] else np.nan for r in rows])
    bx, by, bz = g(base_rows, "x"), g(base_rows, "y"), g(base_rows, "z_ceil")
    es, ex, ey = g(ed_rows, "s"), g(ed_rows, "x"), g(ed_rows, "y")
    ez, ed_dem = g(ed_rows, "z_ceil"), g(ed_rows, "z_dem")
    dist, idx = cKDTree(np.c_[bx, by]).query(np.c_[ex, ey])
    bz_m = np.where(dist < 2.0, bz[idx], np.nan)
    bz_m = bz_m + float(edf._DATUM.get("dz_applied_m", 0.0))
    ax.plot(es, ed_dem, lw=0.9, color=C["surface"], label="Surface DEM")
    ax.plot(es, bz_m, lw=0.9, color=C["dlio"],
            label="Ceiling, rim-centroid solve (4-DOF)")
    ax.plot(es, ez, lw=0.9, color=C["tube"], label="Ceiling, slice registration")
    ax.set_ylabel(ELEV.replace(" (", "\n("))
    ax.set_xlabel("Distance along tube $s$ (m)")
    ax.set_ylim(None, np.nanmax(ed_dem) + 6.0)
    ax.legend(loc="upper left", ncols=3, columnspacing=1.0)


def reg_two_slam(ax, chain):
    zones = chain["ceiling_offset_validation"]
    mids, means, stds = [], [], []
    for k, v in zones.items():
        lo, hi = k[1:].split("-")
        mids.append((float(lo) + float(hi)) / 2)
        means.append(v["mean"]); stds.append(v["std"])
    ax.errorbar(mids, means, yerr=stds, fmt="o", ms=3, lw=0.9, capsize=2.5,
                color=C["envelope"])
    ax.axhline(0, color="0.5", lw=0.6, ls="--")
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel("Distance along tube $s$ (m)")
    # sign as in chain_flf_to_ed.py: p99 ceiling of DLIO minus FAST-LIO
    ax.set_ylabel("Ceiling offset,\nDLIO minus FAST-LIO (m)")


def _class_marker(cls):
    """The two rare classes get a larger, black-edged marker so they stay
    visible under a deuteranopia simulation (pink vs green collapses)."""
    if cls in ("inconsistent", "multipass"):
        return dict(s=12, lw=0.4, edgecolor="k", zorder=4)
    return dict(s=4, lw=0)


def cons_column_stat(ax, s, ghost, klass):
    for cls, col in CLASS_COLOUR.items():
        m = klass == cls
        if m.any():
            ax.scatter(s[m], np.clip(ghost[m], 0, 1), color=col,
                       **_class_marker(cls))
    ax.axhline(0.05, color="k", lw=0.7, ls="--")
    ax.annotate("5% screen threshold", (0.99, 0.055),
                xycoords=("axes fraction", "data"), ha="right", va="bottom",
                fontsize=6, color="0.25")
    ax.set_ylabel("Returns above\nsurface (fraction)")
    ax.set_xlabel("Distance along tube $s$ (m)")


def cons_classification(ax, s, x, y, klass):
    to_uv = edf.principal_frame(np.column_stack([x, y]))
    u, v = to_uv(np.column_stack([x, y]))
    for cls, col in CLASS_COLOUR.items():
        m = klass == cls
        if m.any():
            ax.scatter(u[m], v[m], color=col, label=cls, **_class_marker(cls))
    for tick in range(0, int(s.max()) + 1, 50):
        j = int(np.argmin(np.abs(s - tick)))
        ax.annotate(f"{tick} m", (u[j], v[j]), fontsize=6, color="0.35",
                    xytext=_tick_offset(tick, s.max(), 6),
                    textcoords="offset points",
                    ha=_tick_ha(tick, s.max()))
    ax.set_aspect("equal")
    ax.margins(x=0.02, y=0.12)
    ax.set_xlabel("Along principal axis (m)")
    ax.set_ylabel("Across (m)")
    ax.legend(loc="lower left", ncols=4, markerscale=2.0, columnspacing=1.0)


def ed5_registration_consistency():
    rep = json.load(open(OUT / "registration_report.json"))
    T = np.array(json.load(open(OUT_LEGACY / "transform_ed_slice.json"))["matrix"])
    val = json.load(open(OUT / "registration_validation.json"))
    chain = json.load(open(OUT_LEGACY / "transform_flf_ed_chain.json"))
    s, x, y, z_dem, z_ceil, tau, ghost, klass = edf.load_roof()

    fig = new_fig(168.0)
    gs = fig.add_gridspec(4, 2, height_ratios=[1.25, 1.0, 0.75, 1.15],
                          width_ratios=[1.55, 1])
    ax_a = fig.add_subplot(gs[0, 0]); reg_constellation(ax_a, rep, T)
    ax_b = fig.add_subplot(gs[0, 1]); reg_floor_hist(ax_b, val)
    ax_c = fig.add_subplot(gs[1, :]); reg_ceiling_two_solutions(ax_c)
    ax_d = fig.add_subplot(gs[2, 0]); reg_two_slam(ax_d, chain)
    ax_e = fig.add_subplot(gs[2, 1]); cons_column_stat(ax_e, s, ghost, klass)
    ax_f = fig.add_subplot(gs[3, :]); cons_classification(ax_f, s, x, y, klass)
    for a in (ax_c, ax_d, ax_e):
        a.set_xlim(-5, 310)
    layout(fig)
    freeze_layout(fig)
    place_letters(fig, [(ax_a, "a"), (ax_b, "b"), (ax_c, "c"), (ax_d, "d"),
                        (ax_e, "e"), (ax_f, "f")])
    save(fig, "ed5_registration_consistency")


# ----------------------------------------------------------------------
# ED 6: centreline, sections, sinuosity, aspect ratio
# ----------------------------------------------------------------------
def ed6_centreline_morphometry():
    from scipy.spatial import cKDTree
    s, P, Tn = edf.load_centreline()
    _, s_all, x, y, _, _, _, klass = edf.load_roof()
    shell = edf.load_xyz(edf.SHELL_PCD, keep=350_000)
    bands = skylight_bands(s_all, klass)
    to_uv = edf.principal_frame(P[:, :2])
    su, sv = to_uv(shell)
    cu, cv = to_uv(P)

    fig = new_fig(170.0)
    gs = fig.add_gridspec(4, 6, height_ratios=[2.2, 1.0, 1.1, 1.0])

    # (a) plan view, (b) side view
    ax_a = fig.add_subplot(gs[0, :])
    ax_a.scatter(su, sv, s=0.2, color=C["shell"], lw=0, rasterized=True)
    ax_a.plot(cu, cv, lw=1.2, color=C["tube"], label="L1-medial centreline")
    sky = klass == "skylight"
    # markers only: a line would bridge the intact gaps between the three
    # apertures and read as one continuous skylight reach
    ax_a.scatter(cu[sky], cv[sky], s=7, color=C["skylight"], lw=0, zorder=4,
                 label="Skylight stations")
    for tick in range(0, int(s.max()) + 1, 50):
        j = int(np.argmin(np.abs(s - tick)))
        ax_a.annotate(f"{tick} m", (cu[j], cv[j]), fontsize=6, color="0.35",
                      xytext=_tick_offset(tick, s.max(), 7),
                      textcoords="offset points",
                      ha=_tick_ha(tick, s.max()))
    ax_a.set_aspect("equal")
    ax_a.set_xlim(cu.min() - 25, cu.max() + 25)
    ax_a.set_ylim(cv.min() - 12, cv.max() + 12)
    ax_a.set_ylabel("Across (m)")
    ax_a.legend(loc="upper left")

    ax_b = fig.add_subplot(gs[1, :], sharex=ax_a)
    ax_b.scatter(su, shell[:, 2], s=0.2, color=C["shell"], lw=0, rasterized=True)
    ax_b.plot(cu, P[:, 2], lw=1.2, color=C["tube"])
    ax_b.set_ylim(P[:, 2].min() - 15, P[:, 2].max() + 15)
    ax_b.set_xlabel("Along principal axis (m)")
    ax_b.set_ylabel(ELEV.replace(" (", "\n("))

    # (c-h) cross-section gallery from the 10 cm map
    cloud = edf.load_xyz(edf.SECT_PCD)
    tree = cKDTree(cloud)
    stations = [25, 75, 125, 175, 225, 275]
    morpho = {float(r["s"]): r for r in
              csv.DictReader(open(OUT / "morphometry.csv"))}
    sec_axes, secs = [], []
    for st in stations:
        j = int(np.argmin(np.abs(s - st)))
        stn, T = P[j], Tn[j]
        idx = tree.query_ball_point(stn, 20.0)
        d = cloud[idx] - stn
        sec = d[np.abs(d @ T) <= 0.375]
        uvec = np.cross(T, [0, 0, 1.0]); uvec = uvec / max(np.linalg.norm(uvec), 1e-6)
        vvec = np.cross(T, uvec)
        if vvec[2] < 0:
            vvec = -vvec
        secs.append((j, sec @ uvec, sec @ vvec))
    # one common frame for the six sections so their sizes compare directly
    hw = max(np.percentile(np.abs(u), 99.5) for _, u, _ in secs) + 1.0
    vlo = min(np.percentile(v, 0.5) for _, _, v in secs) - 1.0
    vhi = max(np.percentile(v, 99.5) for _, _, v in secs) + 1.0
    for i, (j, su_, sv_) in enumerate(secs):
        ax = fig.add_subplot(gs[2, i])
        ax.scatter(su_, sv_, s=0.4, color=C["tube"], lw=0, rasterized=True)
        row = morpho.get(round(float(s[j]), 4)) or \
            morpho.get(min(morpho, key=lambda k: abs(k - s[j])))
        area = f"{float(row['area']):.0f}" if row and row["area"] else "-"
        ax.set_aspect("equal")
        ax.set_xlim(-hw, hw); ax.set_ylim(vlo, vhi)
        ax.text(1.0, 1.02, f"$s$ = {s[j]:.0f} m, $A$ = {area} m$^2$",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=6,
                color="0.25")
        if i == 0:
            ax.set_ylabel("Height (m)")
        else:
            ax.tick_params(labelleft=False)
        ax.set_xlabel("Across (m)")
        sec_axes.append(ax)

    # (i) sinuosity over 50 m windows, (j) aspect ratio (as fig_morphometry)
    rows = list(csv.DictReader(open(OUT / "morphometry.csv")))
    ok = [r for r in rows if r["area"]]
    ms = np.array([float(r["s"]) for r in ok])
    eta = np.array([float(r["eta"]) for r in ok])
    Wn = 50
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    sinu = np.full(len(P), np.nan)
    for k in range(len(P)):
        lo, hi = max(0, k - Wn // 2), min(len(P) - 1, k + Wn // 2)
        if cum[hi] - cum[lo] >= 0.9 * Wn:
            sinu[k] = (cum[hi] - cum[lo]) / max(np.linalg.norm(P[hi] - P[lo]), 1e-6)
    ax_i = fig.add_subplot(gs[3, :3])
    ax_i.plot(np.arange(len(P)), sinu, lw=0.9, color=C["tube"])
    for lo, hi in bands:
        ax_i.axvspan(lo, hi, color=C["skylight"], alpha=0.15, lw=0)
    ax_i.set_ylabel("Sinuosity\n(50 m window)")
    ax_i.set_xlabel("Distance along tube $s$ (m)")
    ax_j = fig.add_subplot(gs[3, 3:])
    for lo, hi in bands:
        ax_j.axvspan(lo, hi, color=C["skylight"], alpha=0.15, lw=0)
    ax_j.plot(ms, eta, lw=0.9, color=C["tube"])
    ax_j.set_ylabel("Aspect ratio $\\eta$")
    ax_j.set_xlabel("Distance along tube $s$ (m)")
    for a in (ax_i, ax_j):
        a.set_xlim(-5, 310)

    layout(fig, h_pad=4.0, wspace=0.03)
    freeze_layout(fig)
    place_letters(fig, [(ax_a, "a"), (ax_b, "b")] +
                  [(ax, "cdefgh"[i]) for i, ax in enumerate(sec_axes)] +
                  [(ax_i, "i"), (ax_j, "j")])
    save(fig, "ed6_centreline_morphometry")


# ----------------------------------------------------------------------
# ED 7: surface model validation and orbital transfer test
# ----------------------------------------------------------------------
def ed7_dem_validation():
    import rasterio
    from pyproj import Transformer
    ORIGIN = np.array([479158.0, 7089826.0])
    val = json.loads((OUT / "national_dem_check.json").read_text())["arcticdem"]
    z, xmin, ymin, res = odt.load_dem()
    # the DEM cache is x-major: z[ix, iy] (HANDOFF_VERIFICATION rule)
    step = 4
    ix, iy = np.meshgrid(np.arange(0, z.shape[0], step),
                         np.arange(0, z.shape[1], step), indexing="ij")
    zz = z[ix, iy]
    ok = np.isfinite(zz)
    xs = ORIGIN[0] + xmin + ix[ok] * res
    ys = ORIGIN[1] + ymin + iy[ok] * res
    tr = Transformer.from_crs(32627, 3413, always_xy=True)
    bdx, bdy = val["best_shift_m"]
    X, Y = tr.transform(xs + bdx, ys + bdy)
    with rasterio.open(ROOT / "maps/national_dem_site.tif") as r:
        v = np.array([q[0] for q in r.sample(np.c_[X, Y])], dtype=np.float64)
        if r.nodata is not None:
            v[v == r.nodata] = np.nan
    v[np.abs(v) > 1e4] = np.nan
    m = np.isfinite(v)
    d = zz[ok][m] - v[m]
    d0 = d - np.median(d)

    fig = new_fig(150.0)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], width_ratios=[1.25, 1])
    ax_a = fig.add_subplot(gs[0, 0])
    sc = ax_a.scatter(xs[m] - ORIGIN[0], ys[m] - ORIGIN[1], c=d0, s=2.0,
                      cmap="RdBu_r", vmin=-1, vmax=1, lw=0, rasterized=True)
    ax_a.set_aspect("equal")
    ax_a.set_anchor("W")
    ax_a.set_xlabel("x (m, local frame)")
    ax_a.set_ylabel("y (m)")
    cax = ax_a.inset_axes([1.04, 0.1, 0.05, 0.8])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("Photogrammetric surface minus ArcticDEM (m)")
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5, length=2)

    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.hist(d0, bins=np.arange(-1.5, 1.55, 0.1), color=C["tube"], alpha=0.9)
    ax_b.axvline(0, color="k", lw=0.6)
    ax_b.annotate(f"std {val['std_raw_m']:.2f} m\n"
                  f"doming {val['doming_amplitude_m']:.2f} m\n"
                  f"tilt {val['tilt_m_per_km']:.2f} m km$^{{\\mathrm{{-1}}}}$",
                  (0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
                  fontsize=6, color="0.25")
    ax_b.set_xlabel("Elevation difference (m)")
    ax_b.set_ylabel("Samples")

    rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    g = lambda k: np.array([float(r[k]) if r[k] else np.nan for r in rows])
    s, x, y, zceil, tau0 = g("s"), g("x"), g("y"), g("z_ceil"), g("tau")
    it = np.array([r["class"] for r in rows]) == "intact"
    ax_c = fig.add_subplot(gs[1, :])
    o = np.argsort(s[it])
    ax_c.fill_between(s[it][o], tau0[it][o] - 1.5, tau0[it][o] + 1.5,
                      color=C["intact"], alpha=0.15, lw=0,
                      label="Baseline ± $\\sigma_\\tau$")
    ax_c.plot(s[it][o], tau0[it][o], color=C["intact"], lw=1.1,
              label="Drone surface model (baseline)")
    for name, grid, sig, col, ls in (
            ("HiRISE-like 1 m / 0.3 m", 1.0, 0.3, "#CC79A7", (0, (2, 1))),
            ("LROC-NAC-like 5 m / 1 m", 5.0, 1.0, "#D55E00", "-")):
        rng = np.random.default_rng(104)
        zd = odt.degrade(z, res, grid, sig, 0.0, rng)
        zt = odt.sample(zd, xmin, ymin, res, x, y)
        taud = zt - zceil
        ax_c.plot(s[it][o], taud[it][o], lw=0.7, color=col, alpha=0.85, ls=ls,
                  label=name)
    ax_c.axhline(0, color="k", lw=0.5)
    ax_c.set_xlim(-5, 310)
    ax_c.set_xlabel("Distance along tube $s$ (m)")
    ax_c.set_ylabel("Overburden $\\tau$ (m)")
    ax_c.legend(ncols=2, loc="upper left")
    layout(fig, wspace=0.04)
    freeze_layout(fig)
    place_letters(fig, [([ax_a, cb.ax], "a"), (ax_b, "b"), (ax_c, "c")])
    save(fig, "ed7_dem_validation")


# ----------------------------------------------------------------------
# ED 8: the 1970 survey and the road crossing
# ----------------------------------------------------------------------
def ed8_survey1970():
    from matplotlib.colors import LinearSegmentedColormap
    sv = json.load(open(OUT / "survey1970.json"))
    rc = json.load(open(OUT / "road_crossing.json"))
    G = np.load(OUT / "survey1970_georef.npz")
    Cl, s, sky, ms, mw, mh = s70.load_ours()
    road = s70.osm_road()
    shell = edf.load_xyz(edf.SHELL_PCD, keep=300_000)
    fig3 = s70.load_gray("ellis-4_upright.png")
    sim3 = s70.Sim(float(G["fig3_sim_s"]), G["fig3_sim_R"], G["fig3_sim_t"])
    xlim, ylim = (1130, 1440), (540, 930)
    # scan ink in a neutral dark tone (the v9 dark red paired with the green
    # medial line; Nature asks for no red/green pairs)
    ink = LinearSegmentedColormap.from_list(
        "ink", [(0.15, 0.15, 0.15, 1.0), (0.15, 0.15, 0.15, 0.0)])
    cols = {"A_Munger_1954": "#CC79A7", "B_Cambridgeshire_1969": "#E69F00",
            "C_SMCC_1970": "#D55E00"}
    labels = {"A_Munger_1954": "Munger 1954",
              "B_Cambridgeshire_1969": "Cambridgeshire 1969",
              "C_SMCC_1970": "SMCC 1970 (Fig. 4C)"}
    # line styles carry the distinction where the colours collapse under a
    # deuteranopia simulation (pink/green and orange/vermillion pairs)
    lss = {"A_Munger_1954": (0, (3, 1, 1, 1)),
           "B_Cambridgeshire_1969": (0, (3, 1.5)), "C_SMCC_1970": "-"}
    GREEN = "#009E73"

    fig = new_fig(170.0)
    gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.0, 1.7],
                          width_ratios=[1.15, 1.15, 1])
    # (a) Fig. 3 plan warped onto the shell + centreline
    ax_a = fig.add_subplot(gs[0:2, 0])
    ax_a.scatter(shell[:, 0], shell[:, 1], s=0.2, color=C["shell"], lw=0,
                 rasterized=True)
    plan = np.full_like(fig3, 255)
    plan[470:960, 380:2470] = fig3[470:960, 380:2470]
    w3, ext = s70.warp_raster(plan, sim3, xlim, ylim, res=0.3)
    ax_a.imshow(w3, extent=ext, origin="lower", cmap=ink, vmin=110, vmax=235,
                interpolation="bilinear", zorder=2)
    P3 = G["fig3_medial_local"]
    ax_a.plot(P3[:, 0], P3[:, 1], color=GREEN, lw=0.7, ls=":", zorder=3,
              label="Fig. 3 medial (digitised)")
    ax_a.plot(road[:, 0], road[:, 1], color="0.25", lw=0.9, ls="--",
              label="Road (Route 39)")
    ax_a.plot(Cl[:, 0], Cl[:, 1], color=C["tube"], lw=1.1, label="Lidar centreline")
    ax_a.scatter(sky[:, 0], sky[:, 1], s=22, facecolor="none",
                 edgecolor=C["skylight"], lw=1.0, zorder=5, label="Skylights S1-S3")
    ax_a.set_xlim(*xlim); ax_a.set_ylim(*ylim); ax_a.set_aspect("equal")
    ax_a.set_xlabel("Local east (m)"); ax_a.set_ylabel("Local north (m)")
    ax_a.legend(loc="lower left", fontsize=5.5, frameon=True, framealpha=0.9,
                edgecolor="none", borderpad=0.3)

    # (b) the three Fig. 4 plans + map-sheet medial
    ax_b = fig.add_subplot(gs[0:2, 1])
    ax_b.scatter(shell[:, 0], shell[:, 1], s=0.2, color=C["shell"], lw=0,
                 rasterized=True)
    ax_b.plot(road[:, 0], road[:, 1], color="0.25", lw=0.9, ls="--")
    for k in cols:
        P = G[f"fig4_{k}"]
        ax_b.plot(P[:, 0], P[:, 1], color=cols[k], lw=0.9, ls=lss[k],
                  label=f"{labels[k]}, RMS {sv['fig4_1971']['plans'][k]['fit']['rms_m']:.1f} m")
    Pm = G["map_medial_local"]
    ax_b.plot(Pm[:, 0], Pm[:, 1], color=GREEN, lw=0.9, marker="o", ms=1.8,
              markevery=8, label="SMCC 1970 map sheet, tie-point fit")
    ax_b.plot(Cl[:, 0], Cl[:, 1], color=C["tube"], lw=1.1, label="Lidar centreline")
    dp = G["map_danger"]
    ax_b.plot(dp[0], dp[1], marker="*", ms=7, color="k", ls="none",
              label="1970 'point of greatest danger'")
    ax_b.set_xlim(*xlim); ax_b.set_ylim(*ylim); ax_b.set_aspect("equal")
    ax_b.set_xlabel("Local east (m)")
    ax_b.tick_params(labelleft=False)
    ax_b.legend(loc="lower left", fontsize=5.5, frameon=True, framealpha=0.9,
                edgecolor="none", borderpad=0.3)

    # (c) offset vs distance
    ax_c = fig.add_subplot(gs[0, 2])
    ss, o, e = G["offsets_s"], G["offsets_m"], G["offsets_env"]
    ax_c.fill_between(ss, 0, np.maximum(e, 0), color="0.85",
                      label="Ellis 1970 stated 2% of distance")
    ax_c.plot(ss, o, "o-", ms=2.5, lw=0.9, color=GREEN,
              label="Map sheet 1970 vs lidar")
    for k in cols:
        P = G[f"fig4_{k}"]
        d4 = [s70.nearest_dist(Cl[int(np.argmin(np.abs(s - st)))][None], P)[0][0]
              for st in ss]
        ax_c.plot(ss, d4, lw=0.7, color=cols[k], alpha=0.9, ls=lss[k])
    ax_c.axvspan(0, 90, color=C["skylight"], alpha=0.08, lw=0)
    ax_c.set_ylim(0, 45)
    ax_c.text(45, 1.01, "tie-point reach", ha="center", va="bottom",
              fontsize=5.5, color=C["skylight"],
              transform=matplotlib.transforms.blended_transform_factory(
                  ax_c.transData, ax_c.transAxes))
    ax_c.set_xlabel("Distance along the lidar centreline (m)")
    ax_c.set_ylabel("Horizontal offset (m)")
    h, l = ax_c.get_legend_handles_labels()
    h += [Line2D([], [], color=cols[k], lw=0.7, ls=lss[k]) for k in cols]
    l += [labels[k] for k in cols]
    ax_c.legend(h, l, fontsize=5.5, loc="upper right", frameon=True,
                framealpha=0.95, edgecolor="none", borderpad=0.3)

    # (d) width comparison
    ax_d = fig.add_subplot(gs[1, 2])
    rows = list(csv.DictReader(open(OUT / "survey1970_widths.csv")))
    ws = np.array([float(r["s"]) for r in rows])
    wl = np.array([float(r["width_lidar_m"]) for r in rows])
    w70 = np.array([float(r["width_1970_m"]) for r in rows])
    use = np.array([r["usable"] == "True" for r in rows])
    ax_d.plot(ws, wl, color=C["tube"], lw=1.0, label="Lidar (p1-p99 extent)")
    ax_d.plot(ws[use], w70[use], ".", ms=3, color=GREEN, label="1970 map sheet walls")
    wc = sv["width_comparison"]["ratio_1970_over_lidar"]
    ax_d.set_ylim(0, None)
    ax_d.annotate(f"1970/lidar median {wc['median']:.2f}\n"
                  f"(p10-p90 {wc['p10']:.2f}-{wc['p90']:.2f})",
                  (0.03, 0.72), xycoords="axes fraction", ha="left", va="top",
                  fontsize=6, color="0.25")
    ax_d.set_xlabel("Distance along the lidar centreline (m)")
    ax_d.set_ylabel("Passage width (m)")
    ax_d.legend(loc="upper left", fontsize=5.5)

    # (e) road-crossing long section
    ax = fig.add_subplot(gs[2, :])
    cells = rc["cells"]
    sa = np.array([c["s"] for c in cells])
    rrows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    rs = np.array([float(r["s"]) for r in rrows])
    rdem = np.array([float(r["z_dem"]) for r in rrows])
    rce = np.array([float(r["z_ceil"]) for r in rrows])
    m = (rs >= 250) & (rs < sa.min())
    ax.plot(np.r_[rs[m], sa], np.r_[rdem[m], [c["dem_cache"] for c in cells]],
            color=C["surface"], lw=1.1, label="Surface DEM (0.5 m cache)")
    ax.plot(sa, [c["dem_geotiff"] for c in cells], color=C["surface"], lw=0.7,
            ls=":", label="Surface DEM (6 cm GeoTIFF)")
    mc = rs >= 250
    ax.plot(rs[mc], rce[mc], color=C["tube"], lw=1.1, label="Ceiling, surveyed stations")
    mcol = {"flf": C["tube"], "tube": C["dlio"], "slf": "#56B4E9"}
    mlab = {"flf": "FAST-LIO", "tube": "DLIO", "slf": "second flight"}
    mmk = {"flf": "o", "tube": "s", "slf": "^"}
    for mk in ("flf", "tube", "slf"):
        zc = np.array([c["maps"][mk]["ceiling_pipeline"] for c in cells], float)
        n = np.array([c["maps"][mk]["n_returns"] for c in cells])
        zs = np.array([c["maps"][mk]["ceiling_boot_std"] for c in cells], float)
        okc = np.isfinite(zc) & (sa > rc["survey_end_s"] - 1)
        ax.errorbar(sa[okc], zc[okc], yerr=np.nan_to_num(zs[okc]), fmt=mmk[mk],
                    ms=2.2, lw=0.6, color=mcol[mk], label=f"Ceiling, {mlab[mk]}")
        sp = okc & (n < 30)
        ax.plot(sa[sp], zc[sp], "o", ms=4.5, mfc="none", mec=mcol[mk], lw=0.5)
        ns = np.array([not (c["maps"][mk].get("roof_sampled", True) and
                            c["maps"][mk].get("pipeline_consistent", True))
                       for c in cells]) & okc
        ax.plot(sa[ns], zc[ns], "x", ms=4.5, color=mcol[mk], lw=0.7)
        fl = np.array([c["maps"][mk]["floor_p1"] for c in cells], float)
        ax.plot(sa[okc], fl[okc], ".", ms=2, color=mcol[mk], alpha=0.5)
    r0, r1 = rc["road"]["crest_band_s"]
    ax.axvspan(r0, r1, color="0.85", lw=0, label="Road embankment crest")
    crest_z = rc["road"]["embankment"]["cache"]["crest_m"]
    xr = (250, sa.max() + 5)
    ax.axhspan(crest_z - 12.0, crest_z - 5.5, xmin=(r0 - xr[0]) / (xr[1] - xr[0]),
               xmax=(r1 - xr[0]) / (xr[1] - xr[0]), color=C["skylight"],
               alpha=0.18, lw=0, label="1970: roof 5.5-12 m under road")
    ax.axvline(rc["survey_end_s"], color="k", lw=0.5, ls=":")
    ax.set_xlim(*xr)
    ax.set_ylim(155.0, 184)
    ax.text(rc["survey_end_s"] - 0.8, 183.6, "survey end", fontsize=6,
            va="top", ha="right")
    ax.set_xlabel("Distance along the tube axis (m; beyond 301 m the "
                  "return-following axis)")
    ax.set_ylabel(ELEV)
    ax.legend(fontsize=5.5, loc="lower left", ncol=3, columnspacing=0.8,
              handlelength=1.4, handletextpad=0.5)
    ax_e = ax

    layout(fig, h_pad=4.0, wspace=0.04)
    freeze_layout(fig)
    place_letters(fig, [(ax_a, "a"), (ax_b, "b"), (ax_c, "c"), (ax_d, "d"),
                        (ax_e, "e")])
    save(fig, "ed8_survey1970")


# ----------------------------------------------------------------------
# ED tables
# ----------------------------------------------------------------------
def mission_rows():
    ev = json.load(open(OUT_ED / "mission_evidence.json"))["bags"]
    f1, f2 = ev["first_long_flight"], ev["second_long_flight"]

    def n_iter(b):
        hits = b["rosout"]["keyword_hits"]
        it = [int(m.group(1)) for kw in hits for _, _, msg in hits[kw]
              for m in [re.search(r"Planning iter -> (\d+)", msg)] if m]
        return max(it)

    def col(b):
        bt = b["bt_transitions"]
        return [
            f"{b['pose']['duration_s']:.0f}",
            f"{b['pose']['path_length_m']:.0f}",
            f"{b['pose']['max_dist_from_start_m']:.1f}",
            f"{b['pose']['z_range_m']:.1f}",
            f"{n_iter(b)}",
            f"{bt['n_total']:,}",
            f"{bt['failures_by_node'].get('RunLocalPlanner', 0)}",
            f"{bt['by_node']['CheckPlanValidity']} / "
            f"{bt['failures_by_node'].get('CheckPlanValidity', 0)}",
            f"{b['execution_complete']['n']}",
            f"{b['homing_path']['n_waypoints'][0]} → "
            f"{b['homing_path']['n_waypoints'][-1]}",
        ]

    names = ["Recorded mission duration (s)",
             "Flown path length, odometry (m)",
             "Furthest excursion from launch (m)",
             "Vertical range of trajectory (m)",
             "Planning iterations",
             "Behaviour tree transitions (1 Hz tick)",
             "Local-planner failures (single recovered burst)",
             "Segment-validity checks run / failed",
             "Planned segments executed to completion",
             "Return-route waypoints maintained (first → last)"]
    c1, c2 = col(f1), col(f2)
    return [[n, a, b] for n, a, b in zip(names, c1, c2)]


def planner_rows():
    pp = json.load(open(OUT_LEGACY / "planner_params.json"))
    mp, rb, pl = pp["mapping"], pp["robot"], pp["planner"]
    sm = pl["sensor_model"]
    rmax = float(re.search(r"max_range ([0-9.]+) m", sm).group(1))
    rres = float(re.search(r"ray resolution ([0-9.]+) rad", sm).group(1))
    bb = rb["bounding_box_m"]
    return [
        ["Vehicle bounding box", "", f"{bb[0]:g} × {bb[1]:g} × {bb[2]:g} m"],
        ["Clearance inflation margin", "δ", f"{rb['clearance_extension_m']:g} m"],
        ["Path safety extension", "", f"{rb['safety_extension_m']:g} m"],
        ["Gain-model sensor range", "$r_{\\max}$", f"{rmax:g} m"],
        ["Gain-model ray spacing", "", f"{rres:.2f} rad"],
        ["Local sampling space", "", f"±{pl['local_bounded_space_m']:g} m cube"],
        ["Maximum local graph size", "", f"{pl['num_vertices_max']} vertices"],
        ["Maximum edge length", "", f"{pl['edge_length_max_m']:.1f} m"],
        ["Executed path truncation", "", f"{pl['traverse_length_max_m']:g} m"],
        ["Gain weight (applied to the normalised gain ĝ)", "$w_g$",
         f"{pl['unknown_voxel_gain']:g}"],
        ["Heading-change penalty weight", "$w_\\psi$", f"{pl['path_direction_penalty']:g}"],
        ["Path-length penalty weight", "$w_d$", f"{pl['path_length_penalty']:g}"],
        ["Cruise / homing speed limit", "",
         f"{pl['v_max_mps']:.1f} / {pl['v_homing_mps']:.1f} m s$^{{\\mathrm{{-1}}}}$"],
        ["Yaw rate limit", "", f"{pl['yaw_rate_max_radps']:g} rad s$^{{\\mathrm{{-1}}}}$"],
        ["Exploration time budget", "", f"{int(pl['exploration_time_budget_s']):,} s"],
        ["Onboard map voxel size", "$\\ell$", f"{mp['tsdf_voxel_size_m']:g} m"],
        ["Map ray integration length", "$r_{\\mathrm{map}}$", f"{mp['max_ray_length_m']:g} m"],
        ["ESDF truncation distance", "", f"{mp['esdf_max_distance_m']:.1f} m"],
    ]


def ed_table1_missions():
    rows = mission_rows()
    fig = new_fig(60.0)
    h = draw_table(fig, 0, 1.0, 180.0, ["", "Flight 1", "Flight 2"], rows,
                   col_w=[100, 40, 40], align=["left", "center", "center"],
                   row_h=3.4)
    fig.set_size_inches(W_MM * MM, (h + 2.0) * MM)
    _refit_axes(fig, 60.0, h + 2.0)
    save(fig, "ed_table1_missions")


def ed_table2_planner_params():
    """Two booktabs tables on one page: the platform and survey
    configuration (formerly ED Fig. 2c) above the planner parameters."""
    fig = new_fig(170.0)
    prow, secs = platform_rows()
    y = 1.0
    h1 = draw_table(fig, 0, y, 180.0, ["Component", "Configuration"], prow,
                    col_w=[28, 152], section_rows=secs)
    y += h1 + 5.0
    rows = planner_rows()
    h2 = draw_table(fig, 0, y, 180.0, ["Parameter", "Symbol", "Value"], rows,
                    col_w=[100, 30, 50], align=["left", "center", "center"],
                    row_h=3.4)
    y += h2 + 1.0
    fig.set_size_inches(W_MM * MM, y * MM)
    _refit_axes(fig, 170.0, y)
    save(fig, "ed_table2_planner_params")


# ----------------------------------------------------------------------
ITEMS = {
    "ed1": ed1_site_photos,
    "ed2": ed2_platform_architecture,
    "ed3": ed3_planner_behaviour,
    "ed4": ed4_deployments,
    "ed5": ed5_registration_consistency,
    "ed6": ed6_centreline_morphometry,
    "ed7": ed7_dem_validation,
    "ed8": ed8_survey1970,
    "ed_table1": ed_table1_missions,
    "ed_table2": ed_table2_planner_params,
}

if __name__ == "__main__":
    only = None
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1].split(",")
    FIGS.mkdir(parents=True, exist_ok=True)
    print("Nature Extended Data items ->", FIGS)
    for k, fn in ITEMS.items():
        if only and k not in only:
            continue
        if k == "ed4" and not only:
            # analysis/make_nature_ed4_deployments.py writes the real ED4
            # (thesis figures) to the same file; the placeholder here is
            # only produced on an explicit --only ed4
            print("[ed4] skipped (placeholder; see make_nature_ed4_deployments.py)")
            continue
        print(f"[{k}]")
        fn()
