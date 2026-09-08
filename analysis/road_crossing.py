"""Overburden at the road crossing from the far-range lidar returns.

The surveyed centreline ends at s = 301 m heading due west, 25-35 m short of
the Route 39 embankment. The interior maps carry sparse far-range returns
beyond the last vehicle position (DLIO map to x = 1223, FAST-LIO to 1248,
second flight to 1161). This script
  1. extends the tube axis beyond s = 301 two ways: the end tangent of the
     last 20 m of centreline, and a return-following axis (per 2.5 m x-bin
     centroid of the returns west of x = 1280 in the band y 800-840);
  2. in 2.5 m cells along the return-following axis from s = 280 to the
     last cell with returns, computes per map the number of returns, the
     pipeline ceiling (run_morphometry.section_metrics: top of the innermost
     radial layer over bins within ~45 deg of up, 0.75 m half-slab, 30 m
     search radius; falls back to the slab maximum), a wider-window p99
     ceiling with a bootstrap standard deviation, the p1 floor, the surface
     DEM (0.5 m per-cell-maximum cache, x-major z[ix, iy]; and Birgir's
     6 cm GeoTIFF sampled at the cell), the road-crest flag and the
     overburden DEM - ceiling. Cells with fewer than 30 returns are flagged;
  3. locates the 1970 "point of greatest danger" (survey1970.json) on the
     axis and reports the overburden across the crossing reach;
  4. measures the embankment height from the DEM and compares with the
     1970 figures (5.5-12 m roof to road surface incl. ~1.5 m banking; error
     budget 0.2% of traverse distance + 1.5 m surface heights).
Outputs: ANALYSIS_OUT/road_crossing.json, road_crossing_cells.csv,
         fig_road_crossing.pdf
Run:  ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
      .venv/bin/python analysis/road_crossing.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
MAPS = ROOT / "maps"
sys.path.insert(0, str(ROOT / "analysis"))
import run_morphometry as rm  # noqa: E402  (section_metrics, pipeline constants)
from lava_pcd.io.pcd_reader import BinaryPcdReader  # noqa: E402

ORIGIN = np.array([479158.0, 7089826.0])
DEM_NPZ = OUT / "dem_grid_aerial_isn16_ortho_dem.npz"
TIF = Path(os.environ.get("ISN16_DEM_TIF", str(
    ROOT.parent / "Birgir_data_and_papers_08_sep" / "20250828_Raufarholshellir_DEM_ISN16.tif")))
ORTHO = ROOT.parent / "Birgir_data_and_papers_08_sep" / "20250828_Raufarholshellir_ORTHO_ISN16.tif"
MAP_FILES = {"flf": MAPS / "flf_10cm_ed_ortho.pcd",     # offline FAST-LIO, primary
             "tube": MAPS / "tube_10cm_ed_ortho.pcd",   # onboard DLIO
             "slf": MAPS / "slf_10cm_ed_ortho.pcd"}     # second flight
CELL = 2.5            # cell length along the axis (m)
HALF_CELL = CELL / 2
TRANSVERSE = 15.0     # transverse half-window for the cell statistics (m)
MIN_RETURNS = 30
N_BOOT = 300
S_START = 280.0
X_BIN_WEST = 1280.0   # return-following axis starts here
Y_BAND = (795.0, 845.0)
ROAD_HALF_WIDTH = 5.0  # crest ~10 m wide in the DEM
# Ellis 1970 special report
ELLIS = dict(roof_to_road_surface_m=(5.5, 12.0), banking_m=1.5, rock_min_m=4.0,
             rock_worst_case_m=5.5, height_error_m=1.0, roof_fallen_m=6.0,
             crossing_length_m=100.0, vertical_frac=0.002, surface_height_m=1.5)


def load_xyz(path: Path) -> np.ndarray:
    parts = []
    with BinaryPcdReader(path) as r:
        for c in r.chunks():
            parts.append(c[:, :3].astype(np.float64))
    X = np.vstack(parts)
    return X[np.isfinite(X).all(1)]


def dem_at(dem, x, y):
    ix = ((np.asarray(x) - dem["xmin"]) / dem["res"]).astype(int)
    iy = ((np.asarray(y) - dem["ymin"]) / dem["res"]).astype(int)
    z = dem["z"]; out = np.full(np.shape(ix), np.nan)
    ok = (ix >= 0) & (ix < z.shape[0]) & (iy >= 0) & (iy < z.shape[1])
    out[ok] = z[ix[ok], iy[ok]]
    return out


def geotiff_sampler():
    import rasterio
    from pyproj import Transformer
    tr = Transformer.from_crs(32627, 8088, always_xy=True)
    r = rasterio.open(TIF)

    def f(x, y):
        X, Y = tr.transform(np.asarray(x) + ORIGIN[0], np.asarray(y) + ORIGIN[1])
        rr, cc = rasterio.transform.rowcol(r.transform, X, Y)
        rr, cc = np.atleast_1d(rr), np.atleast_1d(cc)
        out = np.full(len(rr), np.nan)
        for i, (a, b) in enumerate(zip(rr, cc)):
            v = r.read(1, window=rasterio.windows.Window(int(b), int(a), 1, 1))
            out[i] = np.nan if v.size == 0 or v[0, 0] == r.nodata else float(v[0, 0])
        return out
    return f, float(r.res[0])


def osm_road():
    from pyproj import Transformer
    d = json.load(open(OUT / "osm_way_620673065.json"))
    nodes = {e["id"]: (e["lon"], e["lat"]) for e in d["elements"] if e["type"] == "node"}
    way = [e for e in d["elements"] if e["type"] == "way"][0]
    tr = Transformer.from_crs(4326, 32627, always_xy=True)
    return np.array([tr.transform(*nodes[n]) for n in way["nodes"]]) - ORIGIN


def dense(P, step=0.25):
    d = np.r_[0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]
    t = np.arange(0, d[-1], step)
    return np.column_stack([np.interp(t, d, P[:, 0]), np.interp(t, d, P[:, 1])]), t


def main() -> None:
    rows = list(csv.DictReader(open(OUT / "centreline.csv")))
    g = lambda k: np.array([float(r[k]) for r in rows])
    s, X, Y, Z = g("s"), g("x"), g("y"), g("z")
    roof = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    d = np.load(DEM_NPZ)
    dem = dict(z=d["z"], xmin=float(d["xmin"]), ymin=float(d["ymin"]), res=float(d["res"]))
    tif_at, tif_res = geotiff_sampler()
    road = osm_road()
    roadd, _ = dense(road, 0.25)
    sv = json.load(open(OUT / "survey1970.json"))
    G = np.load(OUT / "survey1970_georef.npz")
    danger = np.array(sv["map_sheet_1970"]["danger_point"]["local_m"])
    med70 = G["map_medial_local"]

    clouds = {k: load_xyz(p) for k, p in MAP_FILES.items()}
    trees2 = {k: cKDTree(v[:, :2]) for k, v in clouds.items()}
    trees3 = {k: cKDTree(v) for k, v in clouds.items()}

    # ---- 1. axis beyond s = 301 ----
    last = s >= s.max() - 20
    ctr = np.column_stack([X[last], Y[last]]).mean(0)
    _, _, Vt = np.linalg.svd(np.column_stack([X[last], Y[last]]) - ctr)
    tdir = Vt[0] / np.linalg.norm(Vt[0])
    if np.dot(tdir, [X[-1] - X[-2], Y[-1] - Y[-2]]) < 0:
        tdir = -tdir
    end = np.array([X[-1], Y[-1]])
    tang_axis = np.array([end + tdir * k for k in np.arange(0, 80, CELL)])
    # return-following axis: per 2.5 m x-bin centroid of returns in the band
    bins = {}
    for k, P in clouds.items():
        m = (P[:, 0] < X_BIN_WEST) & (P[:, 1] > Y_BAND[0]) & (P[:, 1] < Y_BAND[1])
        Q = P[m]
        edges = np.arange(X_BIN_WEST, Q[:, 0].min() - CELL, -CELL) if len(Q) else []
        rowsb = []
        for e0 in edges:
            mm = (Q[:, 0] <= e0) & (Q[:, 0] > e0 - CELL)
            if mm.sum() >= 5:
                rowsb.append(dict(x=float(e0 - HALF_CELL), y=float(np.median(Q[mm, 1])),
                                  n=int(mm.sum()), zmin=float(Q[mm, 2].min()), zmax=float(Q[mm, 2].max())))
        bins[k] = rowsb
    # union axis: all maps' bins, y = return-weighted median per x
    allb = {}
    for k, rb in bins.items():
        for r in rb:
            allb.setdefault(r["x"], []).append((r["y"], r["n"]))
    xs = sorted(allb.keys(), reverse=True)
    ys = [float(np.average([a for a, _ in allb[x]], weights=[n for _, n in allb[x]])) for x in xs]
    from scipy import ndimage as ndi
    ys = ndi.median_filter(np.array(ys), size=3, mode="nearest")
    ext = np.column_stack([xs, ys])
    # join to the centreline: keep bins west of the last station
    ext = ext[ext[:, 0] < X[-1] - 1.0]
    axis = np.vstack([np.column_stack([X, Y]), ext])
    axis_d, axis_t = dense(axis, 0.25)
    s_axis = axis_t  # arc length from s = 0 (centreline is 1 m stations)
    # tangents (horizontal)
    T = np.gradient(axis_d, axis=0); T /= np.linalg.norm(T, axis=1)[:, None]

    # ---- 2. cells ----
    cells = []
    s_cells = np.arange(S_START, s_axis[-1], CELL)
    zc_prev = float(Z[-1])
    for sc in s_cells:
        i = int(np.argmin(np.abs(s_axis - sc)))
        p, t = axis_d[i], T[i]
        n2 = np.array([-t[1], t[0]])
        cell = dict(s=float(sc), x=float(p[0]), y=float(p[1]), tx=float(t[0]), ty=float(t[1]),
                    dem_cache=float(dem_at(dem, [p[0]], [p[1]])[0]),
                    dem_geotiff=float(tif_at([p[0]], [p[1]])[0]),
                    dist_to_road_m=float(np.min(np.hypot(*(roadd - p).T))),
                    on_centreline=bool(sc <= s.max()), maps={})
        cell["under_road_crest"] = bool(cell["dist_to_road_m"] <= ROAD_HALF_WIDTH)
        for k, P in clouds.items():
            idx = trees2[k].query_ball_point(p, np.hypot(TRANSVERSE, HALF_CELL) + 0.5)
            Q = P[idx] if len(idx) else np.zeros((0, 3))
            if len(Q):
                dxy = Q[:, :2] - p
                along = dxy @ t; across = dxy @ n2
                m = (np.abs(along) <= HALF_CELL) & (np.abs(across) <= TRANSVERSE)
                Q = Q[m]
            n = int(len(Q))
            rec = dict(n_returns=n, ceiling_pipeline=np.nan, ceiling_p99=np.nan, ceiling_max=np.nan,
                       floor_p1=np.nan, ceiling_boot_std=np.nan, sparse=bool(n < MIN_RETURNS),
                       pipeline_fallback_max=None)
            if n >= 5:
                z = Q[:, 2]
                rec["ceiling_p99"] = float(np.percentile(z, 99)); rec["ceiling_max"] = float(z.max())
                rec["floor_p1"] = float(np.percentile(z, 1))
                rng = np.random.default_rng(int(sc * 10))
                rec["ceiling_boot_std"] = float(np.std([np.percentile(rng.choice(z, n), 99) for _ in range(N_BOOT)]))
                # pipeline ceiling: section_metrics at a 3-D station (z = cell median)
                st = np.array([p[0], p[1], float(np.median(z))])
                T3 = np.array([t[0], t[1], 0.0])
                sm = rm.section_metrics(P, st, T3, trees3[k], P)
                if sm is not None:
                    rec["ceiling_pipeline"] = float(sm["z_ceil"])
                    rec["pipeline_n_slab"] = int(sm["n_pts"])
                    rec["pipeline_fallback_max"] = bool(abs(sm["z_ceil"] - Q[:, 2].max()) < 1e-6 and n < 30)
            for key in ("ceiling_pipeline", "ceiling_p99"):
                rec["overburden_" + key.split("_")[1] + "_cache"] = float(cell["dem_cache"] - rec[key]) if np.isfinite(rec[key]) else np.nan
                rec["overburden_" + key.split("_")[1] + "_geotiff"] = float(cell["dem_geotiff"] - rec[key]) if np.isfinite(rec[key]) else np.nan
            cell["maps"][k] = rec
        cells.append(cell)
    # trim trailing cells with no returns in any map
    while cells and all(cells[-1]["maps"][k]["n_returns"] == 0 for k in MAP_FILES):
        cells.pop()
    # single-sided view: a cell whose p99 lies more than 3 m below the running
    # maximum of the ceiling has no roof returns (shadowed); flag it so the
    # 'overburden' there is not read as roof cover
    for k in MAP_FILES:
        run = -np.inf
        for c in cells:
            rec = c["maps"][k]
            z = rec["ceiling_p99"]
            if np.isfinite(z):
                rec["roof_sampled"] = bool(z >= run - 3.0)
                run = max(run, z)
            else:
                rec["roof_sampled"] = False
            # the pipeline's innermost-layer ceiling can lock onto a near
            # partial layer in a single-sided far-range cell; flag it when it
            # sits more than 3 m below the cell's p99
            zp = rec["ceiling_pipeline"]
            rec["pipeline_consistent"] = bool(np.isfinite(zp) and np.isfinite(z) and zp >= z - 3.0)

    # ---- 3. road crest band and embankment height along the axis ----
    dem_prof = np.array([c["dem_cache"] for c in cells]); sa = np.array([c["s"] for c in cells])
    under = np.array([c["under_road_crest"] for c in cells])
    # crest band from the OSM distance along the dense axis
    dist_axis = np.array([np.min(np.hypot(*(roadd - q).T)) for q in axis_d[::4]])
    s4 = s_axis[::4]
    band = s4[dist_axis <= ROAD_HALF_WIDTH]
    crest_band = [float(band.min()), float(band.max())] if len(band) else [np.nan, np.nan]
    # embankment: profile perpendicular to the road through the axis crossing
    jc = int(np.argmin(dist_axis)); pc = axis_d[::4][jc]
    jr = int(np.argmin(np.hypot(*(roadd - pc).T)))
    rt = roadd[min(jr + 8, len(roadd) - 1)] - roadd[max(jr - 8, 0)]; rt /= np.linalg.norm(rt)
    rn = np.array([-rt[1], rt[0]])
    u = np.arange(-40, 40.01, 0.5)
    prof_pts = roadd[jr] + u[:, None] * rn
    prof_cache = dem_at(dem, prof_pts[:, 0], prof_pts[:, 1])
    prof_tif = tif_at(prof_pts[:, 0], prof_pts[:, 1])
    def emb(prof):
        crest = np.nanmax(prof[np.abs(u) <= 6])
        flank = (np.abs(u) >= 15) & (np.abs(u) <= 30) & np.isfinite(prof)
        west = np.nanmedian(prof[(u <= -15) & (u >= -30)]); east = np.nanmedian(prof[(u >= 15) & (u <= 30)])
        # detrend with a line through the flanks (the flow surface slopes)
        cf = np.polyfit(u[flank], prof[flank], 1)
        det = prof - np.polyval(cf, u)
        above = np.abs(u) <= 12
        wid = u[above & (det > 0.3)]
        return dict(crest_m=float(crest), flow_surface_west_m=float(west), flow_surface_east_m=float(east),
                    height_above_mean_flow_m=float(crest - 0.5 * (west + east)),
                    height_above_west_m=float(crest - west), height_above_east_m=float(crest - east),
                    height_above_flank_trend_m=float(np.nanmax(det[np.abs(u) <= 6])),
                    crest_width_m=float(np.ptp(wid)) if len(wid) else np.nan)
    embank = dict(profile_through=pc.round(2).tolist(), road_normal=rn.round(4).tolist(),
                  cache=emb(prof_cache), geotiff=emb(prof_tif),
                  profile=dict(u_m=u.tolist(), dem_cache=np.round(prof_cache, 3).tolist(), dem_geotiff=np.round(prof_tif, 3).tolist()))

    # ---- 4. summaries: last stations, crest band, 1970 point ----
    def cells_in(lo, hi):
        return [c for c in cells if lo - 1e-6 <= c["s"] <= hi + 1e-6]
    def summarise(cs, key="ceiling_pipeline"):
        res = {}
        for k in MAP_FILES:
            ob = np.array([c["dem_cache"] - c["maps"][k][key] for c in cs], float)
            ob_t = np.array([c["dem_geotiff"] - c["maps"][k][key] for c in cs], float)
            n = np.array([c["maps"][k]["n_returns"] for c in cs])
            sd = np.array([c["maps"][k]["ceiling_boot_std"] for c in cs], float)
            rs = np.array([c["maps"][k].get("roof_sampled", True) for c in cs])
            if key == "ceiling_pipeline":
                rs &= np.array([c["maps"][k].get("pipeline_consistent", True) for c in cs])
            ok = np.isfinite(ob) & rs
            res[k] = dict(n_cells=int(ok.sum()), n_cells_roof_not_sampled=int((np.isfinite(ob) & ~rs).sum()),
                          n_cells_sparse=int((ok & (n < MIN_RETURNS)).sum()),
                          returns_total=int(n.sum()), returns_min=int(n[ok].min()) if ok.any() else 0,
                          returns_max=int(n.max()),
                          overburden_cache_min=float(ob[ok].min()) if ok.any() else np.nan,
                          overburden_cache_max=float(ob[ok].max()) if ok.any() else np.nan,
                          overburden_cache_median=float(np.median(ob[ok])) if ok.any() else np.nan,
                          overburden_geotiff_min=float(np.nanmin(ob_t[ok])) if ok.any() else np.nan,
                          overburden_geotiff_max=float(np.nanmax(ob_t[ok])) if ok.any() else np.nan,
                          ceiling_boot_std_median=float(np.nanmedian(sd[ok])) if ok.any() else np.nan,
                          ceiling_boot_std_max=float(np.nanmax(sd[ok])) if ok.any() else np.nan)
        return res
    last_return_s = {k: max([c["s"] for c in cells if c["maps"][k]["n_returns"] >= 5] or [np.nan]) for k in MAP_FILES}
    # 1970 danger point: on the axis?
    dd = np.hypot(*(axis_d - danger).T); jd = int(np.argmin(dd))
    # 1970 medial arc at the crest and at the danger point
    arc70 = G["map_medial_arc"]
    d_end = np.hypot(*(med70 - danger).T); j70 = int(np.argmin(d_end))
    danger_info = dict(local_m=danger.round(2).tolist(),
                       dist_to_extended_axis_m=float(dd[jd]), nearest_axis_s=float(s_axis[jd]),
                       beyond_last_return_by_map_m={k: float(np.hypot(*(danger - axis_d[int(np.argmin(np.abs(s_axis - last_return_s[k])))]).T)) if np.isfinite(last_return_s[k]) else np.nan for k in MAP_FILES},
                       arc_along_1970_medial_m=float(arc70[j70]),
                       dem_cache=float(dem_at(dem, [danger[0]], [danger[1]])[0]),
                       dem_geotiff=float(tif_at([danger[0]], [danger[1]])[0]),
                       overburden_measurable=bool(dd[jd] < 10 and any(np.isfinite(cells[int(np.argmin(np.abs(sa - s_axis[jd])))]["maps"][k]["ceiling_pipeline"]) for k in MAP_FILES)) if len(cells) else False)
    # 1970 error budget at the crossing
    dist70 = float(sv["map_sheet_1970"]["danger_point"]["arc_along_1970_medial_m"])
    budget = dict(traverse_distance_m=dist70, vertical_0p2pct_m=ELLIS["vertical_frac"] * dist70,
                  surface_height_m=ELLIS["surface_height_m"], stated_height_error_m=ELLIS["height_error_m"],
                  linear_sum_m=ELLIS["vertical_frac"] * dist70 + ELLIS["surface_height_m"],
                  quadrature_m=float(np.hypot(ELLIS["vertical_frac"] * dist70, ELLIS["surface_height_m"])))
    # overburden trend on the surveyed centreline (from roof_thickness.csv)
    tau_last = [(float(r["s"]), float(r["tau"])) for r in roof if float(r["s"]) >= 250 and r["class"] == "intact"]
    # 1970 road-line vs modern road
    road70 = G["map_road70"]
    sep = [float(np.min(np.hypot(*(roadd - q).T))) for q in road70]
    out = dict(
        inputs=dict(maps={k: str(v) for k, v in MAP_FILES.items()}, dem_cache=str(DEM_NPZ), geotiff=str(TIF),
                    geotiff_res_m=tif_res, osm_way=620673065, survey1970=str(OUT / "survey1970.json")),
        parameters=dict(cell_m=CELL, transverse_half_window_m=TRANSVERSE, min_returns=MIN_RETURNS, n_boot=N_BOOT,
                        pipeline_slab_half_m=rm.SLAB, pipeline_r_max_m=rm.R_MAX, pipeline_n_bins=rm.N_BINS,
                        x_bin_west=X_BIN_WEST, y_band=list(Y_BAND), road_half_width_m=ROAD_HALF_WIDTH),
        survey_end_s=float(s.max()), survey_end_xy=[float(X[-1]), float(Y[-1])],
        axis=dict(tangent_dir=tdir.round(4).tolist(),
                  tangent_axis_points=tang_axis[::4].round(2).tolist(),
                  return_bins_by_map=bins,
                  return_following_axis=ext.round(2).tolist(),
                  note="cells use the return-following axis joined to the centreline at s = 301"),
        road=dict(crest_band_s=crest_band, crossing_xy=pc.round(2).tolist(),
                  embankment=embank,
                  road_1970_vs_modern_separation_m=dict(min=float(np.min(sep)), max=float(np.max(sep)),
                                                        at_crossing=float(np.min(np.hypot(*(road70 - pc).T))))),
        last_return_s_by_map=last_return_s,
        overburden_last_surveyed=dict(s_tau=tau_last[-5:], max_s250_301=max(t for _, t in tau_last)),
        crest_band=summarise(cells_in(*crest_band)) if np.isfinite(crest_band[0]) else {},
        crest_band_p99=summarise(cells_in(*crest_band), "ceiling_p99") if np.isfinite(crest_band[0]) else {},
        beyond_survey_all=summarise([c for c in cells if c["s"] > s.max()]),
        crossing_reach_100m=summarise(cells_in(crest_band[0] - 50, crest_band[0] + 50)) if np.isfinite(crest_band[0]) else {},
        danger_point_1970=danger_info,
        ellis_1970=ELLIS, ellis_error_budget=budget,
        cells=cells,
    )
    (OUT / "road_crossing.json").write_text(json.dumps(out, indent=1, default=float))
    with open(OUT / "road_crossing_cells.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["s", "x", "y", "dem_cache", "dem_geotiff", "dist_to_road_m", "under_road_crest"] +
                   [f"{k}_{q}" for k in MAP_FILES for q in ("n", "roof_sampled", "ceil_pipeline", "ceil_p99", "boot_std", "floor_p1", "overburden_pipeline_cache", "overburden_pipeline_geotiff")])
        for c in cells:
            w.writerow([f"{c['s']:.1f}", f"{c['x']:.2f}", f"{c['y']:.2f}", f"{c['dem_cache']:.3f}", f"{c['dem_geotiff']:.3f}",
                        f"{c['dist_to_road_m']:.1f}", c["under_road_crest"]] +
                       [v for k in MAP_FILES for v in (c["maps"][k]["n_returns"], c["maps"][k]["roof_sampled"], f"{c['maps'][k]['ceiling_pipeline']:.3f}",
                                                        f"{c['maps'][k]['ceiling_p99']:.3f}", f"{c['maps'][k]['ceiling_boot_std']:.3f}",
                                                        f"{c['maps'][k]['floor_p1']:.3f}", f"{c['maps'][k]['overburden_pipeline_cache']:.3f}",
                                                        f"{c['maps'][k]['overburden_pipeline_geotiff']:.3f}")])
    print("crest band s", crest_band, "crossing", pc.round(1))
    print("embankment cache", {k: round(v, 2) for k, v in embank["cache"].items()})
    print("embankment geotiff", {k: round(v, 2) for k, v in embank["geotiff"].items()})
    print("last return s", last_return_s)
    print("crest band", json.dumps(out["crest_band"], indent=0, default=float))
    print("danger point", json.dumps(danger_info, indent=0, default=float))
    make_figure(out, axis_d, s_axis, clouds, road, med70, danger, dem, tang_axis)


def make_figure(out, axis_d, s_axis, clouds, road, med70, danger, dem, tang_axis):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "analysis"))
    from make_ed_figures import C as COL  # shared palette
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 7,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    cells = out["cells"]; sa = np.array([c["s"] for c in cells])
    fig = plt.figure(figsize=(7.1, 7.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.35], hspace=0.32)
    ax = fig.add_subplot(gs[0])
    rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    rs = np.array([float(r["s"]) for r in rows]); rdem = np.array([float(r["z_dem"]) for r in rows]); rce = np.array([float(r["z_ceil"]) for r in rows])
    m = (rs >= 250) & (rs < S_START)
    dem_s = np.r_[rs[m], sa]; dem_z = np.r_[rdem[m], [c["dem_cache"] for c in cells]]
    ax.plot(dem_s, dem_z, color=COL["surface"], lw=1.4, label="surface DEM (0.5 m cache)")
    ax.plot(sa, [c["dem_geotiff"] for c in cells], color=COL["surface"], lw=0.8, ls=":", label="surface DEM (6 cm GeoTIFF)")
    mc = rs >= 250
    ax.plot(rs[mc], rce[mc], color=COL["tube"], lw=1.4, label="ceiling, surveyed stations")
    mcol = {"flf": COL["tube"], "tube": COL["dlio"], "slf": "#56B4E9"}; mlab = {"flf": "FAST-LIO", "tube": "DLIO", "slf": "second flight"}
    for k in ("flf", "tube", "slf"):
        z = np.array([c["maps"][k]["ceiling_pipeline"] for c in cells], float); n = np.array([c["maps"][k]["n_returns"] for c in cells])
        sd = np.array([c["maps"][k]["ceiling_boot_std"] for c in cells], float); ok = np.isfinite(z) & (sa > out["survey_end_s"] - 1)
        ax.errorbar(sa[ok], z[ok], yerr=np.nan_to_num(sd[ok]), fmt="o", ms=2.8, lw=0.7, color=mcol[k], label=f"{mlab[k]} cell ceiling")
        sp = ok & (n < MIN_RETURNS); ax.plot(sa[sp], z[sp], "o", ms=5.5, mfc="none", mec=mcol[k], lw=0.6)
        ns = np.array([not (c["maps"][k].get("roof_sampled", True) and c["maps"][k].get("pipeline_consistent", True)) for c in cells]) & ok
        ax.plot(sa[ns], z[ns], "x", ms=5, color=mcol[k], lw=0.8)
        fl = np.array([c["maps"][k]["floor_p1"] for c in cells], float); ax.plot(sa[ok], fl[ok], ".", ms=2.5, color=mcol[k], alpha=0.5)
    r0, r1 = out["road"]["crest_band_s"]
    ax.axvspan(r0, r1, color="0.85", lw=0, label="road embankment crest")
    # 1970 band: roof 5.5-12 m below the road surface, drawn below the crest DEM
    crest_z = out["road"]["embankment"]["cache"]["crest_m"]
    ax.axhspan(crest_z - 12.0, crest_z - 5.5, xmin=(r0 - 250) / (sa.max() + 5 - 250), xmax=(r1 - 250) / (sa.max() + 5 - 250), color=COL["skylight"], alpha=0.18, lw=0, label="1970: roof 5.5-12 m below road")
    ax.axvline(out["survey_end_s"], color="k", lw=0.6, ls=":")
    ax.set_xlim(250, sa.max() + 5)
    ax.set_xlabel("distance along the tube axis (m; beyond 301 m the return-following axis)")
    ax.set_ylabel("elevation (m a.s.l., ISH2004)")
    ax.set_title("a  Long section at the road crossing; open circles: < 30 returns; crosses: ceiling not sampled", loc="left")
    ax.legend(frameon=False, fontsize=5.8, loc="center left", ncol=2, bbox_to_anchor=(0.0, 0.5),
              columnspacing=0.8, handlelength=1.6)
    ax.text(out["survey_end_s"] - 0.8, ax.get_ylim()[1] - 0.2, "survey end", fontsize=6.5, va="top", ha="right")

    ax = fig.add_subplot(gs[1])
    xlim, ylim = (1150, 1320), (760, 920)
    # ortho crop (25 cm)
    try:
        import rasterio
        from pyproj import Transformer
        tr = Transformer.from_crs(32627, 8088, always_xy=True)
        xs = np.arange(xlim[0], xlim[1], 0.25); ys = np.arange(ylim[0], ylim[1], 0.25)
        XX, YY = np.meshgrid(xs, ys); Xp, Yp = tr.transform(XX.ravel() + ORIGIN[0], YY.ravel() + ORIGIN[1])
        with rasterio.open(ORTHO) as r:
            rr, cc = rasterio.transform.rowcol(r.transform, Xp, Yp); rr = np.asarray(rr); cc = np.asarray(cc)
            r0_, c0_ = rr.min(), cc.min(); win = rasterio.windows.Window(c0_, r0_, cc.max() - c0_ + 1, rr.max() - r0_ + 1)
            tile = r.read([1, 2, 3], window=win)
        img = tile[:, rr - r0_, cc - c0_].reshape(3, len(ys), len(xs)).transpose(1, 2, 0)
        ax.imshow(img, extent=[xlim[0], xlim[1], ylim[0], ylim[1]], origin="lower", interpolation="bilinear")
    except Exception as exc:  # the ortho is 7 GB and external; the panel still works without it
        print("ortho not drawn:", exc)
    for k, c in (("tube", COL["dlio"]), ("flf", COL["tube"]), ("slf", "#56B4E9")):
        P = clouds[k]; mm = (P[:, 0] < 1290) & (P[:, 1] > 780) & (P[:, 1] < 900)
        ax.scatter(P[mm, 0], P[mm, 1], s=0.6, color=c, lw=0, alpha=0.6, rasterized=True, label=f"far-range returns, {mlab[k]}")
    rows = list(csv.DictReader(open(OUT / "centreline.csv")))
    ax.plot([float(r["x"]) for r in rows], [float(r["y"]) for r in rows], color="w", lw=2.6)
    ax.plot([float(r["x"]) for r in rows], [float(r["y"]) for r in rows], color=COL["tube"], lw=1.4, label="lidar centreline")
    ax.plot(tang_axis[:, 0], tang_axis[:, 1], color="w", lw=1.8); ax.plot(tang_axis[:, 0], tang_axis[:, 1], color="k", lw=0.8, ls="--", label="end-tangent axis")
    ext = np.array(out["axis"]["return_following_axis"])
    ax.plot(ext[:, 0], ext[:, 1], color="w", lw=2.2); ax.plot(ext[:, 0], ext[:, 1], color="k", lw=1.0, label="return-following axis")
    ax.plot(med70[:, 0], med70[:, 1], color="#009E73", lw=1.1, label="1970 map-sheet passage")
    ax.plot(road[:, 0], road[:, 1], color="0.2", lw=1.0, ls="--", label="road (OSM)")
    ax.plot(danger[0], danger[1], marker="*", ms=10, color="k", ls="none", label="1970 'point of greatest danger'")
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal")
    ax.set_xlabel("local east (m)"); ax.set_ylabel("local north (m)")
    ax.set_title("b  Plan: orthomosaic (28 Aug 2025), far-range returns and the 1970 crossing", loc="left")
    leg = ax.legend(frameon=False, fontsize=6.5, loc="center left", bbox_to_anchor=(1.02, 0.5))
    from matplotlib.collections import PathCollection
    for h in leg.legend_handles:
        if isinstance(h, PathCollection):
            h.set_sizes([18.0]); h.set_alpha(1.0)
    fig.savefig(OUT / "fig_road_crossing.pdf", bbox_inches="tight"); plt.close(fig)
    print("wrote", OUT / "fig_road_crossing.pdf")


if __name__ == "__main__":
    main()
