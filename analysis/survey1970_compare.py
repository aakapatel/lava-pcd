"""The lidar survey against the 1970 Shepton Mallet Caving Club theodolite
survey (Ellis 1970 map sheet; Ellis 1971 Fig. 3 and Fig. 4).

Scans (pdftoppm, 300 dpi) are rotated upright and the two halves of the
folded map sheet are stitched. Every pixel coordinate used below was read
off gridded crops of those rasters by hand (helper crops, not committed) or
measured from the pixel profile (scale-bar ticks, wall crossings); they are
recorded here so the run is reproducible from the PDFs alone.

Georeferencing:
  map sheet   scale from its 0-100 m bar; rotation + translation from four
              tie points (entrance line midpoint, three roof-collapse arrow
              tips) matched to our s=0 station and the three skylight
              centroids; a free-scale similarity is fitted as a check.
  Fig. 3      scale from its 0-200 m bar; rotation + translation by a shape
              fit (point-to-polyline ICP) of its digitised main-tube medial
              line to our centreline over the shared reach, anchored at the
              entrance. The collapse marks cannot be told apart at the
              reproduction scale, so this is not an independent tie-point
              fit and is reported as such.
  Fig. 4      scale from its shared bar; each of the three plans (Munger
              1954, Cambridgeshire 1969, SMCC 1970) fitted like Fig. 3.
Medial lines: a hand seed polyline snapped to the drawn walls along the
local perpendicular (midpoint = medial, separation = width).

Outputs (ANALYSIS_OUT): survey1970.json, survey1970_offsets.csv,
survey1970_widths.csv, survey1970_georef.npz (georeferenced rasters and
polylines for the figure).  With --figure (after road_crossing.py):
figures_for_manuscript/ed_survey1970_panel.pdf.

Run:  ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
      .venv/bin/python analysis/survey1970_compare.py [--figure]
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

Image.MAX_IMAGE_PIXELS = None
ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
SCANS = OUT / "scans1970"
BIRGIR = ROOT.parent / "Birgir_data_and_papers_08_sep"
ELLIS_PDF = BIRGIR / "Ellis_1971_The_survey_Raufarholshellir.pdf"
MAP_PDF = BIRGIR / "Map_Raufarholshellir.pdf"
ORIGIN = np.array([479158.0, 7089826.0])

# ----------------------------------------------------------------------
# Digitisation record (pixel coordinates on the upright 300 dpi rasters)
# ----------------------------------------------------------------------
# Map sheet (Ellis 1970): page 1 rotated +90 deg, page 2 rotated -90 deg.
# Stitch: two matched features (page 2 -> page 1 pixel coordinates) give a
# similarity for page 2; the canvas is page 1 shifted down by CANVAS_OY.
MAP_STITCH_PAIRS = dict(A2=(200.0, 748.0), A1=(3150.0, 300.0),
                        B2=(200.0, 3645.0), B1=(3186.0, 3200.0))
CANVAS_OY = 460
MAP_BAR_M = dict(x0=133.0, x100=1258.0)          # metres bar, 0 and 100 m
MAP_BAR_FT = dict(x0=132.0, x300ft=1160.0)       # feet bar, 0 and 300 ft
MAP_ENTRANCE_LINE = [(1063.0, 2440.0), (1040.0, 2550.0)]   # drawn mouth
MAP_COLLAPSE_TIPS = [(1352.0, 2703.0), (1551.0, 2897.0), (1835.0, 3103.0)]
MAP_DANGER_TIP = (5107.0, 2690.0)                # "Point of greatest danger"
MAP_NORTH_ARROW = [(3375.0, 1588.0), (3971.0, 2076.0)]     # TRUE NORTH
# dash centres read on the faint pencilled road line (two edge lines ~40 px
# apart are drawn; readings alternate between them)
MAP_ROAD_DASHES = [(2350, 2225), (2600, 2290), (3040, 2290), (3135, 2307),
                   (3310, 2355), (3515, 2405), (4225, 2542), (4550, 2555),
                   (4950, 2640), (5265, 2645)]
MAP_MEDIAL_SEED = [
    (1052, 2495), (1180, 2560), (1280, 2660), (1400, 2760), (1500, 2865),
    (1600, 3010), (1700, 3095), (1800, 3115), (1900, 3145), (2000, 3150),
    (2100, 3125), (2200, 3120), (2300, 3155), (2400, 3210), (2500, 3260),
    (2600, 3300), (2700, 3370), (2800, 3370), (2900, 3365), (3000, 3370),
    (3100, 3345), (3200, 3355), (3300, 3345), (3400, 3305), (3500, 3260),
    (3600, 3240), (3700, 3210), (3800, 3160), (3850, 3100), (3900, 3040),
    (3950, 2960), (4000, 2870), (4050, 2760), (4100, 2680), (4200, 2640),
    (4300, 2635), (4400, 2600), (4500, 2610), (4600, 2665), (4700, 2740),
    (4800, 2775), (4900, 2785), (5000, 2725), (5075, 2640), (5100, 2560),
    (5100, 2500), (5150, 2400), (5250, 2300), (5350, 2200), (5400, 2150)]
# side-passage junctions where the wall snap can jump (sheet px along x)
MAP_JUNCTION_X = [(2420, 2700), (3780, 4150)]

# Fig. 3 (Ellis 1971, PDF page 4 rotated -90 deg)
FIG3_BAR = {0: 1013.0, 50: 1178.0, 100: 1339.0, 150: 1502.0, 200: 1664.0}
FIG3_ENTRANCE = (435.0, 560.0)
FIG3_NORTH_ARROW = [(186.0, 789.0), (526.0, 1074.0)]
FIG3_SECTION_ENTRANCE_X = 100.0          # extended section, left end
FIG3_SECTION_ROAD_BRACKET_X = (1235.0, 1575.0)   # "Þrengsli Road" arrows
FIG3_REFPOINT1 = dict(circle=(610.0, 650.0), tip=(590.0, 668.0))
_c = [(85, 100), (130, 150), (165, 205), (200, 250), (240, 290), (290, 320),
      (350, 345), (410, 370), (470, 400), (540, 420), (620, 425), (700, 415),
      (770, 400), (830, 370), (880, 330), (930, 295), (990, 290), (1050, 320),
      (1100, 345), (1150, 330), (1200, 280), (1250, 230), (1300, 180),
      (1370, 140), (1440, 125), (1500, 150), (1560, 175), (1620, 190),
      (1700, 190), (1770, 175), (1850, 150), (1930, 135), (2000, 145)]
FIG3_MEDIAL_SEED = [(350 + 1.05 * x, 450 + 1.05 * y) for x, y in _c]

# Fig. 4 (PDF page 5 rotated -90 deg); shared bar, labels 0..200 m
FIG4_BAR = {0: 1795.0, 100: 2205.0, 200: 2600.0}
FIG4_NORTH_ARROW = [(2006.0, 636.0), (2521.0, 1096.0)]
_k = 3508.0 / 1600.0
FIG4_SEEDS = {
    "A_Munger_1954": [(60, 130), (75, 200), (100, 260), (140, 300), (190, 335),
                      (250, 345), (320, 320), (380, 290), (430, 275), (470, 300),
                      (510, 350), (560, 340), (610, 290), (660, 240), (720, 220),
                      (800, 215), (880, 225), (950, 235), (1030, 225),
                      (1120, 215), (1200, 215), (1260, 225)],
    "B_Cambridgeshire_1969": [(60, 380), (80, 430), (110, 480), (140, 520),
                              (180, 560), (220, 590), (270, 600), (320, 590),
                              (370, 565), (430, 540), (490, 555), (540, 580),
                              (580, 575), (620, 540), (660, 500), (710, 475),
                              (760, 480), (810, 510), (860, 540), (920, 545),
                              (980, 540), (1040, 530), (1100, 535),
                              (1150, 555), (1200, 565)],
    "C_SMCC_1970": [(60, 700), (80, 740), (110, 780), (150, 810), (200, 835),
                    (250, 850), (300, 865), (360, 875), (420, 870), (480, 860),
                    (530, 835), (580, 815), (630, 800), (680, 795), (720, 760),
                    (760, 730), (810, 715), (860, 725), (920, 745), (980, 755),
                    (1040, 760), (1100, 765), (1160, 770), (1220, 775),
                    (1290, 760), (1340, 745)],
}
FIG4_SEEDS = {k: [(x * _k, y * _k) for x, y in v] for k, v in FIG4_SEEDS.items()}

ELLIS_HORIZONTAL_FRAC = 0.02      # stated: < 2% of traverse distance
ELLIS_VERTICAL_FRAC = 0.002       # stated: < 0.2%


# ----------------------------------------------------------------------
# rasters
# ----------------------------------------------------------------------
def rasterise() -> None:
    SCANS.mkdir(parents=True, exist_ok=True)
    if not (SCANS / "ellis-4.png").exists():
        subprocess.run(["pdftoppm", "-r", "300", "-f", "4", "-l", "5", "-png",
                        str(ELLIS_PDF), str(SCANS / "ellis")], check=True)
    if not (SCANS / "map-1.png").exists():
        subprocess.run(["pdftoppm", "-r", "300", "-png", str(MAP_PDF),
                        str(SCANS / "map")], check=True)
    for name, rot in (("ellis-4", -90), ("ellis-5", -90), ("map-1", 90),
                      ("map-2", -90)):
        p = SCANS / f"{name}_upright.png"
        if not p.exists():
            Image.open(SCANS / f"{name}.png").rotate(rot, expand=True).save(p)


def load_gray(name: str) -> np.ndarray:
    return np.array(Image.open(SCANS / name).convert("L"))


def stitch_map() -> tuple[np.ndarray, dict]:
    p = SCANS / "map_stitched.png"
    A2, A1 = np.array(MAP_STITCH_PAIRS["A2"]), np.array(MAP_STITCH_PAIRS["A1"])
    B2, B1 = np.array(MAP_STITCH_PAIRS["B2"]), np.array(MAP_STITCH_PAIRS["B1"])
    d2, d1 = B2 - A2, B1 - A1
    s = np.linalg.norm(d1) / np.linalg.norm(d2)
    th = np.arctan2(d1[1], d1[0]) - np.arctan2(d2[1], d2[0])
    R = s * np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    t = A1 - R @ A2 + np.array([0.0, CANVAS_OY])
    info = dict(page2_scale=float(s), page2_rotation_deg=float(np.degrees(th)),
                page2_translation_px=t.tolist(), canvas_oy=CANVAS_OY)
    if p.exists():
        return load_gray("map_stitched.png"), info
    p1 = load_gray("map-1_upright.png").astype(np.float32)
    p2 = load_gray("map-2_upright.png").astype(np.float32)
    W, H = 6520, p1.shape[0] + CANVAS_OY + 20
    canvas = np.full((H, W), 255, np.float32)
    canvas[CANVAS_OY:CANVAS_OY + p1.shape[0], :p1.shape[1]] = p1
    Rinv = np.linalg.inv(R)
    M = np.array([[Rinv[1, 1], Rinv[1, 0]], [Rinv[0, 1], Rinv[0, 0]]])
    off = Rinv @ (-t); off = np.array([off[1], off[0]])
    warp = ndi.affine_transform(p2, M, offset=off, output_shape=(H, W),
                                order=1, cval=255)
    canvas = np.minimum(canvas, warp)
    Image.fromarray(canvas.astype(np.uint8)).save(p)
    return canvas.astype(np.uint8), info


# ----------------------------------------------------------------------
# medial line from a seed polyline
# ----------------------------------------------------------------------
def refine_medial(img, seed, step, search, thr, bg_size, n_iter=2):
    """Snap a seed polyline to the drawn walls: at each resampled point search
    along the local perpendicular for the nearest dark pixel on each side."""
    img = img.astype(float)
    dark = (img - ndi.median_filter(img, size=bg_size)) < -thr
    dark = ndi.binary_dilation(dark, iterations=1)
    seed = np.asarray(seed, float)
    d = np.r_[0, np.cumsum(np.hypot(*np.diff(seed, axis=0).T))]
    t = np.arange(0, d[-1], step)
    P = np.column_stack([np.interp(t, d, seed[:, 0]), np.interp(t, d, seed[:, 1])])
    flagged = np.zeros(len(P), bool)
    width = np.full(len(P), np.nan)
    for _ in range(n_iter):
        T = np.gradient(P, axis=0); T /= np.linalg.norm(T, axis=1)[:, None]
        N = np.column_stack([-T[:, 1], T[:, 0]])
        left = np.full(len(P), np.nan); right = np.full(len(P), np.nan)
        for i, (p, n) in enumerate(zip(P, N)):
            for sgn, arr in ((1, right), (-1, left)):
                for r in range(4, search):
                    q = p + sgn * r * n
                    x, y = int(round(q[0])), int(round(q[1]))
                    if not (0 <= y < dark.shape[0] and 0 <= x < dark.shape[1]):
                        break
                    if dark[y, x]:
                        arr[i] = r; break

        def clean(a):
            m = ndi.median_filter(np.nan_to_num(a, nan=np.nanmedian(a)), size=5)
            bad = np.isnan(a) | (np.abs(a - m) > 0.5 * m + 8)
            a = a.copy(); a[bad] = m[bad]; return a, bad
        left, bl = clean(left); right, br = clean(right)
        P = P + ((right - left) / 2.0)[:, None] * N
        width = left + right
        flagged = bl | br
    return P, width, flagged, t


def smooth_poly(P, size):
    return np.column_stack([ndi.median_filter(P[:, 0], size=size, mode="nearest"),
                            ndi.median_filter(P[:, 1], size=size, mode="nearest")])


def arclen(P):
    return np.r_[0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]


# ----------------------------------------------------------------------
# transforms: image px (x right, y down) -> local metres (x east, y north)
# ----------------------------------------------------------------------
class Sim:
    """p_local = s * R @ (x, -y) + t"""
    def __init__(self, s, R, t):
        self.s, self.R, self.t = float(s), np.asarray(R, float), np.asarray(t, float)

    def __call__(self, p):
        p = np.atleast_2d(np.asarray(p, float))
        v = np.column_stack([p[:, 0], -p[:, 1]])
        return self.s * (self.R @ v.T).T + self.t

    @property
    def rotation_deg(self):
        return float(np.degrees(np.arctan2(self.R[1, 0], self.R[0, 0])))

    def as_dict(self):
        return dict(scale_m_per_px=self.s, px_per_m=1.0 / self.s,
                    rotation_deg=self.rotation_deg, translation_m=self.t.tolist(),
                    convention="p_local = s * R @ (x_px, -y_px) + t")


def fit_similarity(src, dst, scale=None) -> tuple[Sim, np.ndarray]:
    src = np.asarray(src, float); dst = np.asarray(dst, float)
    v = np.column_stack([src[:, 0], -src[:, 1]])
    mv, md = v.mean(0), dst.mean(0)
    A, B = v - mv, dst - md
    U, S, Vt = np.linalg.svd(A.T @ B)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1; R = Vt.T @ U.T
    s = (S.sum() / (A ** 2).sum()) if scale is None else scale
    t = md - s * R @ mv
    sim = Sim(s, R, t)
    return sim, sim(src) - dst


def icp_to_polyline(src_pts, target, scale, init: Sim, n_iter=40,
                    anchor=None, anchor_w=1.0):
    """Rigid (fixed-scale) point-to-nearest-target fit of src_pts (px) onto a
    target polyline (m). Optional anchor: (src_px, dst_m) pair with weight."""
    sim = init
    for _ in range(n_iter):
        P = sim(src_pts)
        # nearest points on the target polyline (dense resample)
        j = np.argmin(((P[:, None, :] - target[None, :, :]) ** 2).sum(-1), axis=1)
        Q = target[j]
        S_, D_ = src_pts, Q
        w = np.ones(len(S_))
        if anchor is not None:
            S_ = np.vstack([S_, [anchor[0]] * int(anchor_w)])
            D_ = np.vstack([D_, [anchor[1]] * int(anchor_w)])
        sim_new, _ = fit_similarity(S_, D_, scale=scale)
        if np.allclose(sim_new.t, sim.t, atol=1e-3) and abs(sim_new.rotation_deg - sim.rotation_deg) < 1e-4:
            sim = sim_new; break
        sim = sim_new
    P = sim(src_pts)
    j = np.argmin(((P[:, None, :] - target[None, :, :]) ** 2).sum(-1), axis=1)
    res = np.hypot(*(P - target[j]).T)
    return sim, res


def dense(target, step=0.5):
    d = arclen(target); t = np.arange(0, d[-1], step)
    return np.column_stack([np.interp(t, d, target[:, 0]), np.interp(t, d, target[:, 1])])


def nearest_dist(P, poly):
    j = np.argmin(((P[:, None, :] - poly[None, :, :]) ** 2).sum(-1), axis=1)
    return np.hypot(*(P - poly[j]).T), j


def north_bearing(sim: Sim, arrow):
    """Bearing (deg east of true north) that the sheet's drawn north arrow
    takes in the local frame after georeferencing."""
    a, b = np.asarray(arrow, float)
    v = np.array([b[0] - a[0], -(b[1] - a[1])])
    u = sim.R @ v
    return float(np.degrees(np.arctan2(u[0], u[1])))


# ----------------------------------------------------------------------
def load_ours():
    rows = list(csv.DictReader(open(OUT / "centreline.csv")))
    C = np.array([[float(r["x"]), float(r["y"])] for r in rows])
    s = np.array([float(r["s"]) for r in rows])
    holes = json.load(open(OUT / "aerial_holes.json"))["holes"]
    sky = np.array([h["centroid"][:2] for h in holes[:3]])
    morpho = list(csv.DictReader(open(OUT / "morphometry.csv")))
    ms = np.array([float(r["s"]) for r in morpho])
    mw = np.array([float(r["width"]) if r["width"] else np.nan for r in morpho])
    mh = np.array([float(r["height"]) if r["height"] else np.nan for r in morpho])
    return C, s, sky, ms, mw, mh


def osm_road():
    d = json.load(open(OUT / "osm_way_620673065.json"))
    from pyproj import Transformer
    nodes = {e["id"]: (e["lon"], e["lat"]) for e in d["elements"] if e["type"] == "node"}
    way = [e for e in d["elements"] if e["type"] == "way"][0]
    tr = Transformer.from_crs(4326, 32627, always_xy=True)
    pts = np.array([tr.transform(*nodes[n]) for n in way["nodes"]]) - ORIGIN
    return pts


def stats(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return dict(n=int(len(v)), median=float(np.median(v)), mean=float(v.mean()),
                min=float(v.min()), max=float(v.max()),
                p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)))


# ----------------------------------------------------------------------
def main() -> None:
    rasterise()
    C, s, sky, ms, mw, mh = load_ours()
    Cd = dense(C, 0.25)
    road = osm_road()
    out: dict = {}

    # ---------------- map sheet ----------------
    canvas, stitch = stitch_map()
    px_per_m = (MAP_BAR_M["x100"] - MAP_BAR_M["x0"]) / 100.0
    px_per_m_ft = (MAP_BAR_FT["x300ft"] - MAP_BAR_FT["x0"]) / (300 * 0.3048)
    ent = np.mean(MAP_ENTRANCE_LINE, axis=0)
    src = np.vstack([ent, MAP_COLLAPSE_TIPS])
    dst = np.vstack([C[0], sky])
    sim_map, res = fit_similarity(src, dst, scale=1.0 / px_per_m)
    sim_map_free, res_free = fit_similarity(src, dst, scale=None)
    resn = np.hypot(*res.T)
    # medial line, width, arc length
    Pm, Wm, Fm, tm = refine_medial(canvas, MAP_MEDIAL_SEED, step=25, search=140,
                                   thr=35, bg_size=41)
    junction = np.zeros(len(Pm), bool)
    for a, b in MAP_JUNCTION_X:
        junction |= (Pm[:, 0] > a) & (Pm[:, 0] < b)
    Pm_loc = sim_map(Pm)
    arc_m = arclen(Pm_loc)
    width_m = Wm / px_per_m
    danger = sim_map(MAP_DANGER_TIP)[0]
    rx = np.array([p[0] for p in MAP_ROAD_DASHES], float)
    ry = np.array([p[1] for p in MAP_ROAD_DASHES], float)
    road_c = np.polyfit(rx, ry, 1)
    road_rms_px = float(np.sqrt(np.mean((ry - np.polyval(road_c, rx)) ** 2)))
    xx = np.linspace(1600, 5700, 60)
    road70 = sim_map(np.column_stack([xx, np.polyval(road_c, xx)]))
    # danger point relative to us
    d_danger_c, j = nearest_dist(danger[None], Cd)
    j_arc = int(np.argmin(np.hypot(*(Pm_loc - danger).T)))
    # where the 1970 medial crosses the modern road (OSM polyline)
    roadd = dense(road, 0.5)
    dm, jr = nearest_dist(Pm_loc, roadd)
    k = int(np.argmin(dm))
    # our last station in 1970 arc terms
    dl, jl = nearest_dist(C[-1][None], Pm_loc)
    # 1970 north arrow in the local frame
    map_north = north_bearing(sim_map, MAP_NORTH_ARROW)
    out["map_sheet_1970"] = dict(
        source=str(MAP_PDF), stitch=stitch,
        scale=dict(px_per_m_metres_bar=float(px_per_m), px_per_m_feet_bar=float(px_per_m_ft),
                   px_per_m_free_fit=float(1.0 / sim_map_free.s),
                   sheet_scale_1_to=float(25.4 / 300 * px_per_m * 1000 / 1000 * 1000)),
        tie_points=dict(
            names=["entrance line midpoint -> s=0 station", "collapse 1 -> S1",
                   "collapse 2 -> S2", "collapse 3 -> S3"],
            sheet_px=src.tolist(), local_m=dst.tolist(),
            residual_m=resn.round(3).tolist(), rms_m=float(np.sqrt((resn ** 2).mean())),
            residual_m_free_scale=np.hypot(*res_free.T).round(3).tolist()),
        transform=sim_map.as_dict(),
        north_arrow_bearing_in_local_frame_deg=map_north,
        road_line=dict(fit_px="y = a x + b", a=float(road_c[0]), b=float(road_c[1]),
                       dash_rms_px=road_rms_px, dash_rms_m=road_rms_px / px_per_m,
                       local_polyline=road70[::6].round(2).tolist()),
        danger_point=dict(sheet_px=list(MAP_DANGER_TIP), local_m=danger.round(2).tolist(),
                          arc_along_1970_medial_m=float(arc_m[j_arc]),
                          straight_from_s0_m=float(np.hypot(*(danger - C[0]))),
                          nearest_our_station_s=float(s[np.argmin(np.hypot(*(C - danger).T))]),
                          dist_to_our_centreline_m=float(d_danger_c[0]),
                          dist_to_modern_road_m=float(nearest_dist(danger[None], roadd)[0][0])),
        medial=dict(n=int(len(Pm)), total_arc_m=float(arc_m[-1]),
                    our_last_station_nearest_1970_arc_m=float(arc_m[jl[0]]),
                    our_last_station_offset_m=float(dl[0]),
                    crosses_modern_road_at_local=Pm_loc[k].round(2).tolist(),
                    crosses_modern_road_at_1970_arc_m=float(arc_m[k]),
                    min_dist_to_modern_road_m=float(dm[k])),
    )
    # sheet scale 1:N from px/m at 300 dpi: 1 m -> px_per_m px -> px_per_m*25.4/300 mm
    out["map_sheet_1970"]["scale"]["sheet_scale_1_to"] = float(1000.0 / (px_per_m * 25.4 / 300))

    # offsets vs distance from the entrance (our s every 10 m), 2% envelope
    offs = []
    for st in np.arange(0, s.max() + 1e-6, 10.0):
        i = int(np.argmin(np.abs(s - st)))
        d, jj = nearest_dist(C[i][None], Pm_loc)
        offs.append(dict(s=float(s[i]), offset_m=float(d[0]),
                         arc_1970_m=float(arc_m[jj[0]]),
                         envelope_2pct_m=ELLIS_HORIZONTAL_FRAC * float(arc_m[jj[0]]),
                         junction=bool(junction[jj[0]])))
    with open(OUT / "survey1970_offsets.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(offs[0].keys())); w.writeheader(); w.writerows(offs)
    o = np.array([r["offset_m"] for r in offs]); e = np.array([r["envelope_2pct_m"] for r in offs])
    ss = np.array([r["s"] for r in offs])
    beyond = ss > 90
    out["offset_1970_vs_lidar"] = dict(
        note="nearest horizontal distance from each lidar station (every 10 m) "
             "to the georeferenced 1970 map-sheet medial line; the first 90 m "
             "are constrained by the tie points, beyond that the comparison is "
             "independent",
        all=stats(o), beyond_tie_points=stats(o[beyond]),
        within_2pct_envelope_frac=float(np.mean(o <= np.maximum(e, 1.0))),
        within_2pct_envelope_frac_beyond_tie_points=float(np.mean((o <= np.maximum(e, 1.0))[beyond])),
        max_offset_s=float(ss[np.argmax(o)]), max_offset_m=float(o.max()),
        offset_at_s=[dict(s=float(a), offset_m=float(b)) for a, b in zip(ss[::5], o[::5])])

    # widths: 1970 wall separation at the nearest medial point vs our width
    wid = []
    for i in range(len(ms)):
        j0 = int(np.argmin(np.abs(s - ms[i])))
        d, jj = nearest_dist(C[j0][None], Pm_loc)
        j1 = jj[0]
        ok = (not Fm[j1]) and (not junction[j1]) and np.isfinite(width_m[j1]) and np.isfinite(mw[i])
        wid.append(dict(s=float(ms[i]), width_lidar_m=float(mw[i]),
                        width_1970_m=float(width_m[j1]), ratio_1970_over_lidar=float(width_m[j1] / mw[i]) if ok else np.nan,
                        usable=bool(ok)))
    with open(OUT / "survey1970_widths.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(wid[0].keys())); w.writeheader(); w.writerows(wid)
    r = np.array([x["ratio_1970_over_lidar"] for x in wid])
    wl = np.array([x["width_lidar_m"] for x in wid]); w70 = np.array([x["width_1970_m"] for x in wid])
    use = np.isfinite(r)
    out["width_comparison"] = dict(
        n_stations=int(use.sum()), ratio_1970_over_lidar=stats(r[use]),
        width_lidar_m=stats(wl[use]), width_1970_m=stats(w70[use]),
        note="1970 width = wall separation of the map-sheet plan measured "
             "perpendicular to its medial line at the point nearest each lidar "
             "station; lidar width = morphometry.csv p1-p99 horizontal extent")
    out["height_comparison"] = dict(
        done=False,
        reason="the elevation ticks of the Fig. 3 extended section and the "
               "map-sheet elevation are unreadable at the reproduction scale "
               "(datum 'roughly 160 m'); no ceiling-floor heights digitised")

    # ---------------- Fig. 3 ----------------
    fig3 = load_gray("ellis-4_upright.png")
    xb = np.array(list(FIG3_BAR.keys()), float); pb = np.array(list(FIG3_BAR.values()))
    f3_px_per_m = float(np.polyfit(xb, pb, 1)[0])
    P3, W3, F3, t3 = refine_medial(fig3, FIG3_MEDIAL_SEED, step=10, search=60, thr=40, bg_size=31)
    P3s = smooth_poly(P3, 9)
    a3 = arclen(P3s) / f3_px_per_m
    # initial guess: entrance -> s=0, rotation from the map-sheet fit
    init = Sim(1.0 / f3_px_per_m, sim_map.R, np.zeros(2))
    init = Sim(init.s, init.R, C[0] - init(FIG3_ENTRANCE)[0])
    use3 = a3 <= 300
    sim3, res3 = icp_to_polyline(P3s[use3], Cd, 1.0 / f3_px_per_m, init,
                                 anchor=(np.array(FIG3_ENTRANCE), C[0]), anchor_w=5)
    P3_loc = sim3(P3s)
    sky3, _ = nearest_dist(sky, P3_loc)
    # road bracket on the extended section
    br = (np.array(FIG3_SECTION_ROAD_BRACKET_X) - FIG3_SECTION_ENTRANCE_X) / f3_px_per_m
    out["fig3_1971"] = dict(
        source=f"{ELLIS_PDF} page 4", px_per_m_bar=f3_px_per_m,
        bar_fit_residual_px=np.round(pb - np.polyval(np.polyfit(xb, pb, 1), xb), 1).tolist(),
        fit=dict(method="fixed-scale ICP of the digitised main-tube medial (first "
                        "300 m) to the lidar centreline, entrance anchored",
                 rms_m=float(np.sqrt((res3 ** 2).mean())), n_points=int(use3.sum()),
                 entrance_residual_m=float(np.hypot(*(sim3(FIG3_ENTRANCE)[0] - C[0]))),
                 skylight_to_fig3_medial_m=sky3.round(2).tolist()),
        transform=sim3.as_dict(),
        north_arrow_bearing_in_local_frame_deg=north_bearing(sim3, FIG3_NORTH_ARROW),
        section_road_bracket_from_entrance_m=br.round(1).tolist(),
        medial_total_arc_m=float(a3[-1]),
        refpoint_1=dict(tip_local=sim3(FIG3_REFPOINT1["tip"])[0].round(2).tolist(),
                        nearest_lidar_s=float(s[np.argmin(np.hypot(*(C - sim3(FIG3_REFPOINT1["tip"])[0]).T))])))

    # ---------------- Fig. 4 ----------------
    fig4 = load_gray("ellis-5_upright.png")
    xb = np.array(list(FIG4_BAR.keys()), float); pb = np.array(list(FIG4_BAR.values()))
    f4_px_per_m = float(np.polyfit(xb, pb, 1)[0])
    out["fig4_1971"] = dict(source=f"{ELLIS_PDF} page 5", px_per_m_bar=f4_px_per_m, plans={})
    fig4_loc = {}
    for key, seed in FIG4_SEEDS.items():
        P4, W4, F4, t4 = refine_medial(fig4, seed, step=12, search=70, thr=40, bg_size=31)
        P4s = smooth_poly(P4, 9)
        a4 = arclen(P4s) / f4_px_per_m
        init = Sim(1.0 / f4_px_per_m, sim3.R, np.zeros(2))
        init = Sim(init.s, init.R, C[0] - init(P4s[0])[0])
        use4 = a4 <= 300
        sim4, res4 = icp_to_polyline(P4s[use4], Cd, 1.0 / f4_px_per_m, init)
        P4_loc = sim4(P4s)
        d_sky, _ = nearest_dist(sky, P4_loc)
        # offsets of our stations to this plan (every 10 m)
        offs4 = []
        for st in np.arange(0, s.max() + 1e-6, 10.0):
            i = int(np.argmin(np.abs(s - st)))
            d, _ = nearest_dist(C[i][None], P4_loc); offs4.append(float(d[0]))
        out["fig4_1971"]["plans"][key] = dict(
            fit=dict(method="fixed-scale ICP to the lidar centreline (first 300 m of "
                            "the plan), entrance cut at the sheet margin so no anchor",
                     rms_m=float(np.sqrt((res4 ** 2).mean())), n_points=int(use4.sum())),
            transform=sim4.as_dict(),
            north_arrow_bearing_in_local_frame_deg=north_bearing(sim4, FIG4_NORTH_ARROW),
            offset_to_lidar_stations_m=stats(offs4),
            skylight_to_plan_medial_m=d_sky.round(2).tolist(),
            medial_total_arc_m=float(a4[-1]))
        fig4_loc[key] = P4_loc
    out["assumptions"] = dict(
        entrance_tie="the sheet's drawn mouth line midpoint is matched to our s=0 "
                     "station (1.8 m of roof above it); the 1970 prime survey point "
                     "was outside the entrance, so this tie carries a few metres "
                     "of definitional slack",
        ellis_accuracy=dict(horizontal_frac=ELLIS_HORIZONTAL_FRAC, vertical_frac=ELLIS_VERTICAL_FRAC,
                            surface_height_m=1.5))
    (OUT / "survey1970.json").write_text(json.dumps(out, indent=2, default=float))
    np.savez_compressed(OUT / "survey1970_georef.npz",
                        map_medial_local=Pm_loc, map_medial_arc=arc_m, map_width_m=width_m,
                        map_flag=Fm | junction, map_road70=road70, map_danger=danger,
                        map_sim_R=sim_map.R, map_sim_t=sim_map.t, map_sim_s=sim_map.s,
                        fig3_medial_local=P3_loc, fig3_sim_R=sim3.R, fig3_sim_t=sim3.t, fig3_sim_s=sim3.s,
                        **{f"fig4_{k}": v for k, v in fig4_loc.items()},
                        offsets_s=ss, offsets_m=o, offsets_env=e)
    print(json.dumps({k: out[k] for k in ("offset_1970_vs_lidar", "width_comparison")}, indent=1, default=float))
    print("map tie residuals m", resn.round(2), "rms", round(float(np.sqrt((resn ** 2).mean())), 2))
    print("map free-scale px/m", round(1 / sim_map_free.s, 3), "bar", round(px_per_m, 3))
    print("danger point local", danger.round(1), "1970 arc", round(float(arc_m[j_arc]), 1))
    print("fig3 rms", round(float(np.sqrt((res3 ** 2).mean())), 2))
    for k, v in out["fig4_1971"]["plans"].items():
        print("fig4", k, "rms", round(v["fit"]["rms_m"], 2), "offset median", round(v["offset_to_lidar_stations_m"]["median"], 1))


# ----------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------
def warp_raster(img: np.ndarray, sim: Sim, xlim, ylim, res=0.25):
    """Resample a scan into the local frame (x east, y north) on a grid."""
    nx = int((xlim[1] - xlim[0]) / res); ny = int((ylim[1] - ylim[0]) / res)
    xs = xlim[0] + (np.arange(nx) + 0.5) * res
    ys = ylim[0] + (np.arange(ny) + 0.5) * res
    XX, YY = np.meshgrid(xs, ys)
    q = np.column_stack([XX.ravel(), YY.ravel()]) - sim.t
    v = (np.linalg.inv(sim.R) @ q.T).T / sim.s
    px, py = v[:, 0], -v[:, 1]
    vals = ndi.map_coordinates(img.astype(np.float32), [py, px], order=1, cval=255.0)
    return vals.reshape(ny, nx), [xlim[0], xlim[1], ylim[0], ylim[1]]


def figure() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    sys.path.insert(0, str(ROOT / "analysis"))
    from make_ed_figures import C as COL, load_xyz, MAPS  # noqa: E402  (shared style)
    plt.rcParams.update({
        "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
        "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
    sv = json.load(open(OUT / "survey1970.json"))
    rc = json.load(open(OUT / "road_crossing.json"))
    G = np.load(OUT / "survey1970_georef.npz")
    C, s, sky, ms, mw, mh = load_ours()
    road = osm_road()
    shell = load_xyz(Path(os.environ.get("EDFIG_SHELL_PCD", str(MAPS / "flf_30cm_ed_ortho.pcd"))), keep=300_000)
    canvas, _ = stitch_map()
    fig3 = load_gray("ellis-4_upright.png")
    sim_map = Sim(float(G["map_sim_s"]), G["map_sim_R"], G["map_sim_t"])
    sim3 = Sim(float(G["fig3_sim_s"]), G["fig3_sim_R"], G["fig3_sim_t"])
    xlim, ylim = (1130, 1440), (540, 930)
    ink = LinearSegmentedColormap.from_list("ink", [(0.55, 0.0, 0.0, 1.0), (0.55, 0.0, 0.0, 0.0)])

    fig = plt.figure(figsize=(7.1, 9.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.75, 1.0, 1.05], hspace=0.42, wspace=0.28)

    # (a) Fig. 3 plan warped onto the shell + centreline
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(shell[:, 0], shell[:, 1], s=0.25, color=COL["shell"], lw=0, rasterized=True)
    plan = np.full_like(fig3, 255); plan[470:960, 380:2470] = fig3[470:960, 380:2470]   # plan only
    w3, ext = warp_raster(plan, sim3, xlim, ylim, res=0.3)
    ax.imshow(w3, extent=ext, origin="lower", cmap=ink, vmin=110, vmax=235, interpolation="bilinear", zorder=2)
    P3 = G["fig3_medial_local"]
    ax.plot(P3[:, 0], P3[:, 1], color="#009E73", lw=0.7, ls=":", zorder=3, label="Fig. 3 medial (digitised)")
    ax.plot(road[:, 0], road[:, 1], color="0.25", lw=1.2, ls="--", label="road (Route 39)")
    ax.plot(C[:, 0], C[:, 1], color=COL["tube"], lw=1.4, label="lidar centreline")
    ax.scatter(sky[:, 0], sky[:, 1], s=28, facecolor="none", edgecolor=COL["skylight"], lw=1.2, zorder=5, label="skylights S1-S3")
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal")
    ax.set_xlabel("local east (m)"); ax.set_ylabel("local north (m)")
    ax.set_title(f"a  Ellis 1971 Fig. 3 over the lidar shell\n    (shape fit, RMS {sv['fig3_1971']['fit']['rms_m']:.1f} m)", loc="left")
    ax.legend(frameon=True, framealpha=0.9, edgecolor="none", loc="lower left", fontsize=6.5)

    # (b) the three Fig. 4 plans + map-sheet medial
    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(shell[:, 0], shell[:, 1], s=0.25, color=COL["shell"], lw=0, rasterized=True)
    ax.plot(road[:, 0], road[:, 1], color="0.25", lw=1.2, ls="--")
    cols = {"A_Munger_1954": "#CC79A7", "B_Cambridgeshire_1969": "#E69F00", "C_SMCC_1970": "#D55E00"}
    labels = {"A_Munger_1954": "Munger 1954", "B_Cambridgeshire_1969": "Cambridgeshire 1969", "C_SMCC_1970": "SMCC 1970 (Fig. 4C)"}
    for k in cols:
        P = G[f"fig4_{k}"]
        ax.plot(P[:, 0], P[:, 1], color=cols[k], lw=1.1, label=f"{labels[k]}, RMS {sv['fig4_1971']['plans'][k]['fit']['rms_m']:.1f} m")
    Pm = G["map_medial_local"]
    ax.plot(Pm[:, 0], Pm[:, 1], color="#009E73", lw=1.1, label=f"SMCC 1970 map sheet, tie-point fit")
    ax.plot(C[:, 0], C[:, 1], color=COL["tube"], lw=1.4, label="lidar centreline")
    dp = G["map_danger"]
    ax.plot(dp[0], dp[1], marker="*", ms=9, color="k", ls="none", label="1970 'point of greatest danger'")
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal")
    ax.set_xlabel("local east (m)")
    ax.set_title("b  Older plans on the same frame\n    (medial lines)", loc="left")
    ax.legend(frameon=True, framealpha=0.9, edgecolor="none", loc="lower left", fontsize=6)

    # (c) offset vs distance
    ax = fig.add_subplot(gs[1, 0])
    ss, o, e = G["offsets_s"], G["offsets_m"], G["offsets_env"]
    ax.fill_between(ss, 0, np.maximum(e, 0), color="0.85", label="Ellis 1970 stated 2% of distance")
    ax.plot(ss, o, "o-", ms=3, lw=1, color="#009E73", label="map sheet 1970 vs lidar")
    for k in cols:
        P = G[f"fig4_{k}"]
        d4 = [nearest_dist(C[int(np.argmin(np.abs(s - st)))][None], P)[0][0] for st in ss]
        ax.plot(ss, d4, lw=0.8, color=cols[k], alpha=0.9)
    ax.axvspan(0, 90, color=COL["skylight"], alpha=0.08, lw=0)
    ax.set_ylim(0, 45)
    ax.text(64, 44, "tie-point\nreach", ha="center", va="top", fontsize=6, color=COL["skylight"],
            bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.9))
    ax.set_xlabel("distance along the lidar centreline (m)"); ax.set_ylabel("horizontal offset (m)")
    ax.set_title("c  1970 traverses against the lidar centreline", loc="left")
    h, l = ax.get_legend_handles_labels()
    from matplotlib.lines import Line2D
    h += [Line2D([], [], color=cols[k], lw=0.8) for k in cols]; l += [labels[k] for k in cols]
    ax.legend(h, l, frameon=False, fontsize=6, loc="upper right")

    # (d) width comparison
    ax = fig.add_subplot(gs[1, 1])
    rows = list(csv.DictReader(open(OUT / "survey1970_widths.csv")))
    ws = np.array([float(r["s"]) for r in rows]); wl = np.array([float(r["width_lidar_m"]) for r in rows])
    w70 = np.array([float(r["width_1970_m"]) for r in rows]); use = np.array([r["usable"] == "True" for r in rows])
    ax.plot(ws, wl, color=COL["tube"], lw=1.2, label="lidar (p1-p99 extent)")
    ax.plot(ws[use], w70[use], ".", ms=3, color="#009E73", label="1970 map sheet walls")
    wc = sv["width_comparison"]["ratio_1970_over_lidar"]
    ax.set_xlabel("distance along the lidar centreline (m)"); ax.set_ylabel("passage width (m)")
    ax.set_title(f"d  Width: 1970/lidar median {wc['median']:.2f} (p10-p90 {wc['p10']:.2f}-{wc['p90']:.2f})", loc="left")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")

    # (e) road-crossing long section (same content as fig_road_crossing.pdf panel a)
    ax = fig.add_subplot(gs[2, :])
    cells = rc["cells"]
    sa = np.array([c["s"] for c in cells])
    rrows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    rs = np.array([float(r["s"]) for r in rrows]); rdem = np.array([float(r["z_dem"]) for r in rrows]); rce = np.array([float(r["z_ceil"]) for r in rrows])
    m = (rs >= 250) & (rs < sa.min())
    ax.plot(np.r_[rs[m], sa], np.r_[rdem[m], [c["dem_cache"] for c in cells]], color=COL["surface"], lw=1.4, label="surface DEM (0.5 m cache)")
    ax.plot(sa, [c["dem_geotiff"] for c in cells], color=COL["surface"], lw=0.8, ls=":", label="surface DEM (6 cm GeoTIFF)")
    mc = rs >= 250
    ax.plot(rs[mc], rce[mc], color=COL["tube"], lw=1.4, label="ceiling, surveyed stations")
    mcol = {"flf": COL["tube"], "tube": COL["dlio"], "slf": "#56B4E9"}
    mlab = {"flf": "FAST-LIO", "tube": "DLIO", "slf": "second flight"}
    for mk in ("flf", "tube", "slf"):
        z = np.array([c["maps"][mk]["ceiling_pipeline"] for c in cells], float)
        n = np.array([c["maps"][mk]["n_returns"] for c in cells])
        zs = np.array([c["maps"][mk]["ceiling_boot_std"] for c in cells], float)
        ok = np.isfinite(z) & (sa > rc["survey_end_s"] - 1)
        ax.errorbar(sa[ok], z[ok], yerr=np.nan_to_num(zs[ok]), fmt="o", ms=2.6, lw=0.7, color=mcol[mk], label=f"ceiling, {mlab[mk]}")
        sp = ok & (n < 30); ax.plot(sa[sp], z[sp], "o", ms=5, mfc="none", mec=mcol[mk], lw=0.6)
        ns = np.array([not (c["maps"][mk].get("roof_sampled", True) and c["maps"][mk].get("pipeline_consistent", True)) for c in cells]) & ok
        ax.plot(sa[ns], z[ns], "x", ms=5, color=mcol[mk], lw=0.8)
        fl = np.array([c["maps"][mk]["floor_p1"] for c in cells], float)
        ax.plot(sa[ok], fl[ok], ".", ms=2, color=mcol[mk], alpha=0.5)
    r0, r1 = rc["road"]["crest_band_s"]
    ax.axvspan(r0, r1, color="0.85", lw=0, label="road embankment crest")
    crest_z = rc["road"]["embankment"]["cache"]["crest_m"]
    xr = (250, sa.max() + 5)
    ax.axhspan(crest_z - 12.0, crest_z - 5.5, xmin=(r0 - xr[0]) / (xr[1] - xr[0]), xmax=(r1 - xr[0]) / (xr[1] - xr[0]),
               color=COL["skylight"], alpha=0.18, lw=0, label="1970: roof 5.5-12 m under road")
    ax.axvline(rc["survey_end_s"], color="k", lw=0.6, ls=":")
    ax.set_xlim(*xr)
    ax.text(rc["survey_end_s"] - 0.8, ax.get_ylim()[1] - 0.2, "survey end", fontsize=6.5, va="top", ha="right")
    ax.set_xlabel("distance along the tube axis (m; beyond 301 m the return-following axis)")
    ax.set_ylabel("elevation (m a.s.l., ISH2004)")
    ax.set_title("e  Road crossing: surface, ceiling and floor per 2.5 m cell (open circles: < 30 returns; crosses: ceiling not sampled)", loc="left")
    ax.set_ylim(155.0, 184)
    ax.legend(frameon=False, fontsize=5.5, loc="lower left", ncol=2, bbox_to_anchor=(0.0, 0.0), columnspacing=0.8, handlelength=1.4, handletextpad=0.5, borderaxespad=0.2)
    FIGS = OUT / "figures_for_manuscript"; FIGS.mkdir(exist_ok=True)
    fig.savefig(FIGS / "ed_survey1970_panel.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote", FIGS / "ed_survey1970_panel.pdf")


if __name__ == "__main__":
    if "--figure" in sys.argv:
        figure()
    else:
        main()
