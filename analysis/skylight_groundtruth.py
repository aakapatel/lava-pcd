"""Ground truth at the skylight rims: Birgir's tape measurements of 15 Sep 2026
against the registered model.

Inputs (ANALYSIS_OUT/groundtruth_2026-09-15/):
  20260915_Raufarholshellir_crust.csv   seven rim points (entrance, S1-S3 N/S)
      with tape thickness of the roof lip (Thickn_1 = uppermost terrace,
      Thickn_2 = thicker terrace behind, read by eye), the rim-to-floor
      height at the northern end (Height_c) and lat/lon/ELEV digitised on
      the August 2025 orthomosaic and DEM.
  20260915_Raufarholshellir_measurements_notes.txt  the field notes.

For each rim point this script computes, in the analysis frame
(EPSG:32627 minus (479158, 7089826), orthometric ISH2004 heights):
  1. the surface elevation from the 0.5 m per-cell-maximum DEM cache and
     from Birgir's 6 cm GeoTIFF, against his ELEV column;
  2. the roof thickness at the rim from the interior lidar: the interior
     returns in a strip outward from the rim point (0.5 to 3.5 m beyond
     the point, 1.5 m half-width, away from the aperture centroid), whose
     highest returns (p99, below the surface) are the underside of the roof
     lip; tau_rim = z_DEM(rim) - that underside. Done for the three
     registered interior maps (offline FAST-LIO primary, onboard DLIO,
     second flight);
  3. the pipeline's own overburden at the intact stations that flank the
     aperture (nearest three intact stations to the rim point);
  4. the rim-to-floor height: z_DEM(rim) minus the lowest interior returns
     within 3 m of the rim point (p2; April 2025, snow cone possible), and
     minus the photogrammetric floor seen through the opening (aerial crop
     points within 3 m of the rim point and more than 1.5 m below the
     surface, p5; August 2025);
  5. the opening length |N - S| against twice the semi-major axis of the
     aerial and interior ellipse detections and the length of the
     aperture-class station run.
Outputs: ANALYSIS_OUT/skylight_groundtruth.json, skylight_groundtruth.csv,
         fig_skylight_groundtruth.pdf
Run:  ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
      .venv/bin/python analysis/skylight_groundtruth.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
from pyproj import Transformer
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
MAPS = ROOT / "maps"
sys.path.insert(0, str(ROOT / "src"))
from lava_pcd.io.pcd_reader import BinaryPcdReader  # noqa: E402

ORIGIN = np.array([479158.0, 7089826.0])
GT_DIR = OUT / "groundtruth_2026-09-15"
GT_CSV = GT_DIR / "20260915_Raufarholshellir_crust.csv"
DEM_NPZ = OUT / "dem_grid_aerial_isn16_ortho_dem.npz"
TIF = Path(os.environ.get("ISN16_DEM_TIF", str(
    ROOT.parent / "Birgir_data_and_papers_08_sep" / "20250828_Raufarholshellir_DEM_ISN16.tif")))
AERIAL_CROP = MAPS / "aerial_isn16_ortho_crop.pcd"
MAP_FILES = {"flf": MAPS / "flf_10cm_ed_ortho.pcd",     # offline FAST-LIO, primary
             "tube": MAPS / "tube_10cm_ed_ortho.pcd",   # onboard DLIO
             "slf": MAPS / "slf_10cm_ed_ortho.pcd"}     # second flight
SIGMA_TAU = 1.5          # systematic datum uncertainty of the profile (m)
STRIP_IN, STRIP_OUT, STRIP_HALF = 0.5, 3.5, 1.5   # rim strip geometry (m)
FLOOR_R = 3.0            # radius for the floor search (m)
N_FLANK = 3              # intact stations per side used for the pipeline tau
# aerial hole id per skylight (aerial_holes.json, v6 detections)
HOLE_ID = {"S1": 0, "S2": 1, "S3": 2}
# aperture-class station runs (roof_thickness.csv class == skylight)
SKY_RUN = {"S1": (29, 36), "S2": (54, 56), "S3": (89, 99)}


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


def fnum(v):
    return None if v is None or (isinstance(v, float) and not np.isfinite(v)) else float(v)


def main() -> None:
    rows = list(csv.DictReader(open(GT_CSV)))
    tr = Transformer.from_crs(4326, 32627, always_xy=True)
    pts = {}
    for r in rows:
        x, y = tr.transform(float(r["LONG"]), float(r["LAT"]))
        pts[r["Point"]] = dict(
            x=x - ORIGIN[0], y=y - ORIGIN[1], elev_birgir=float(r["ELEV"]),
            t1=float(r["Thickn_1"]) if r["Thickn_1"] else np.nan,
            t2=float(r["Thickn_2"]) if r["Thickn_2"] else np.nan,
            h_floor=float(r["Height_c"]) if r["Height_c"] else np.nan)

    d = np.load(DEM_NPZ)
    dem = dict(z=d["z"], xmin=float(d["xmin"]), ymin=float(d["ymin"]), res=float(d["res"]))
    try:
        tif_at, tif_res = geotiff_sampler()
    except Exception as e:  # noqa: BLE001
        print("GeoTIFF not sampled:", e)
        tif_at, tif_res = (lambda x, y: np.full(len(np.atleast_1d(x)), np.nan)), np.nan

    holes = json.load(open(OUT / "aerial_holes.json"))["holes"]
    tube_holes = json.load(open(OUT / "tube_holes.json"))["holes"]
    roof = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    rs = np.array([float(r["s"]) for r in roof])
    rx = np.array([float(r["x"]) for r in roof]); ry = np.array([float(r["y"]) for r in roof])
    rtau = np.array([float(r["tau"]) if r["tau"] != "nan" else np.nan for r in roof])
    rcls = np.array([r["class"] for r in roof])

    clouds = {k: load_xyz(p) for k, p in MAP_FILES.items()}
    trees = {k: cKDTree(v[:, :2]) for k, v in clouds.items()}
    aerial = load_xyz(AERIAL_CROP)
    atree = cKDTree(aerial[:, :2])

    results = []
    for name, p in pts.items():
        sky = name[:2] if name.startswith("S") else "E"
        px, py = p["x"], p["y"]
        z_dem = float(dem_at(dem, [px], [py])[0])
        z_tif = float(tif_at([px], [py])[0])
        rec = dict(point=name, skylight=sky, x=px, y=py, z_dem_cache=z_dem, z_dem_geotiff=z_tif,
                   elev_birgir=p["elev_birgir"], tape_upper=p["t1"], tape_lower=p["t2"],
                   tape_floor=p["h_floor"])
        # outward direction: away from the aperture centroid (entrance: use the
        # local centreline tangent reversed, i.e. toward the entrance trench)
        if sky in HOLE_ID:
            c = np.array(holes[HOLE_ID[sky]]["centroid"][:2])
            u = np.array([px, py]) - c
        else:
            k = int(np.argmin(np.hypot(rx - px, ry - py)))
            u = np.array([px, py]) - np.array([rx[k], ry[k]])
        u = u / (np.linalg.norm(u) + 1e-9)
        v = np.array([-u[1], u[0]])
        rec["outward_dir"] = [float(u[0]), float(u[1])]
        # 2. rim roof underside from each interior map
        for k, X in clouds.items():
            idx = trees[k].query_ball_point([px, py], r=STRIP_OUT + STRIP_HALF)
            if not idx:
                rec[f"tau_rim_{k}"] = np.nan; rec[f"n_rim_{k}"] = 0; continue
            Q = X[idx]
            rel = Q[:, :2] - np.array([px, py])
            a = rel @ u; b = rel @ v
            m = (a >= STRIP_IN) & (a <= STRIP_OUT) & (np.abs(b) <= STRIP_HALF) & (Q[:, 2] < z_dem - 0.2)
            rec[f"n_rim_{k}"] = int(m.sum())
            if m.sum() < 20:
                rec[f"tau_rim_{k}"] = np.nan; continue
            under = np.percentile(Q[m, 2], 99)
            rec[f"z_under_{k}"] = float(under)
            rec[f"tau_rim_{k}"] = float(z_dem - under)
        # 3. pipeline tau at the flanking intact stations
        dist = np.hypot(rx - px, ry - py)
        ok = (rcls == "intact") & np.isfinite(rtau)
        order = np.argsort(np.where(ok, dist, np.inf))[:N_FLANK]
        rec["flank_stations_s"] = [float(rs[i]) for i in order]
        rec["flank_dist_m"] = [float(dist[i]) for i in order]
        rec["tau_flank_mean"] = float(np.mean(rtau[order]))
        rec["tau_flank_min"] = float(np.min(rtau[order]))
        rec["tau_flank_max"] = float(np.max(rtau[order]))
        # 4. rim-to-floor height
        for k, X in clouds.items():
            idx = trees[k].query_ball_point([px, py], r=FLOOR_R)
            if len(idx) < 20:
                rec[f"floor_depth_{k}"] = np.nan; continue
            zf = np.percentile(X[idx, 2], 2)
            rec[f"floor_depth_{k}"] = float(z_dem - zf)
            rec[f"n_floor_{k}"] = len(idx)
        idx = atree.query_ball_point([px, py], r=FLOOR_R)
        A = aerial[idx] if idx else np.zeros((0, 3))
        A = A[A[:, 2] < z_dem - 1.5]
        rec["n_floor_photo"] = int(len(A))
        rec["floor_depth_photo"] = float(z_dem - np.percentile(A[:, 2], 5)) if len(A) >= 20 else np.nan
        results.append(rec)

    # 5. opening lengths
    openings = {}
    for sky, hid in HOLE_ID.items():
        n = pts[f"{sky}N"]; s = pts[f"{sky}S"]
        L = float(np.hypot(n["x"] - s["x"], n["y"] - s["y"]))
        a0, a1 = SKY_RUN[sky]
        openings[sky] = dict(
            rim_to_rim_tape_points_m=L,
            aerial_ellipse_2a_m=2 * float(holes[hid]["semi_major"]),
            aerial_ellipse_2b_m=2 * float(holes[hid]["semi_minor"]),
            interior_ellipse_2a_m=2 * float(tube_holes[hid]["semi_major"]),
            aperture_station_run_m=float(a1 - a0 + 1),
            centroid_to_midpoint_m=float(np.hypot(
                0.5 * (n["x"] + s["x"]) - holes[hid]["centroid"][0],
                0.5 * (n["y"] + s["y"]) - holes[hid]["centroid"][1])))

    # summaries
    def arr(key):
        return np.array([r[key] for r in results], dtype=float)
    tape = arr("tape_upper")
    d_elev = arr("z_dem_cache") - arr("elev_birgir")
    d_tif = arr("z_dem_geotiff") - arr("elev_birgir")
    summ = dict(
        n_points=len(results),
        elev_cache_minus_birgir_m=dict(mean=fnum(np.nanmean(d_elev)), std=fnum(np.nanstd(d_elev)),
                                       max_abs=fnum(np.nanmax(np.abs(d_elev)))),
        elev_geotiff_minus_birgir_m=dict(mean=fnum(np.nanmean(d_tif)), std=fnum(np.nanstd(d_tif)),
                                         max_abs=fnum(np.nanmax(np.abs(d_tif)))),
        geotiff_res_m=fnum(tif_res))
    for key in ["tau_rim_flf", "tau_rim_tube", "tau_rim_slf", "tau_flank_mean"]:
        v = arr(key); m = np.isfinite(v) & np.isfinite(tape)
        diff = v[m] - tape[m]
        summ[f"{key}_minus_tape_upper_m"] = dict(
            n=int(m.sum()), mean=fnum(diff.mean()) if m.any() else None,
            median=fnum(np.median(diff)) if m.any() else None,
            rms=fnum(np.sqrt(np.mean(diff ** 2))) if m.any() else None,
            max_abs=fnum(np.max(np.abs(diff))) if m.any() else None,
            within_sigma_tau=int(np.sum(np.abs(diff) <= SIGMA_TAU)))
    tf = arr("tape_floor")
    for key in ["floor_depth_flf", "floor_depth_tube", "floor_depth_slf", "floor_depth_photo"]:
        v = arr(key); m = np.isfinite(v) & np.isfinite(tf)
        diff = v[m] - tf[m]
        summ[f"{key}_minus_tape_floor_m"] = dict(
            n=int(m.sum()), points=[results[i]["point"] for i in np.where(m)[0]],
            values=[fnum(x) for x in diff],
            mean=fnum(diff.mean()) if m.any() else None,
            rms=fnum(np.sqrt(np.mean(diff ** 2))) if m.any() else None)
    summ["sigma_tau_m"] = SIGMA_TAU
    summ["tape_upper_range_m"] = [fnum(np.nanmin(tape)), fnum(np.nanmax(tape))]
    summ["note"] = ("tape_upper is the uppermost roof terrace at the rim; tape_lower the thicker "
                    "terrace behind it, read by eye (B.V.O. notes). tau_rim_* is z_DEM(rim) minus the "
                    "p99 of interior returns in a strip 0.5-3.5 m outward from the rim point. "
                    "tau_flank is the pipeline overburden at the nearest three intact stations. "
                    "floor_depth_* is z_DEM(rim) minus the p2 of interior returns (April 2025) or "
                    "the p5 of photogrammetric returns (August 2025) within 3 m of the rim point; the "
                    "photogrammetric value is not a floor at the rims (too few returns, rim walls) and is "
                    "reported for completeness only.")

    out = dict(inputs=dict(csv=str(GT_CSV), date_measured="2026-09-15",
                           measured_by="B. V. Oskarsson", maps={k: str(v.name) for k, v in MAP_FILES.items()},
                           aerial=AERIAL_CROP.name, dem=DEM_NPZ.name, geotiff=TIF.name),
               parameters=dict(strip_in_m=STRIP_IN, strip_out_m=STRIP_OUT, strip_half_m=STRIP_HALF,
                               floor_radius_m=FLOOR_R, n_flank=N_FLANK),
               points=[{k: (fnum(v) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()}
                       for r in results],
               openings=openings, summary=summ)
    json.dump(out, open(OUT / "skylight_groundtruth.json", "w"), indent=1)

    cols = ["point", "skylight", "x", "y", "elev_birgir", "z_dem_cache", "z_dem_geotiff",
            "tape_upper", "tape_lower", "tau_rim_flf", "tau_rim_tube", "tau_rim_slf",
            "tau_flank_mean", "tau_flank_min", "tau_flank_max",
            "tape_floor", "floor_depth_flf", "floor_depth_tube", "floor_depth_slf", "floor_depth_photo"]
    with open(OUT / "skylight_groundtruth.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols)
        for r in results:
            w.writerow([r.get(c, "") if not isinstance(r.get(c), float) else f"{r[c]:.3f}" for c in cols])

    # Supplementary Table 1 (LaTeX, booktabs) for the SI
    label = {"Entrance_cave": "Entrance", "S1S": "S1 south", "S1N": "S1 north",
             "S2S": "S2 south", "S2N": "S2 north", "S3S": "S3 south", "S3N": "S3 north"}
    order = ["Entrance_cave", "S1S", "S1N", "S2S", "S2N", "S3S", "S3N"]
    by = {r["point"]: r for r in results}
    def cell(v, nd=1):
        return "" if v is None or not np.isfinite(v) else f"{v:.{nd}f}"
    lines = [
        "%% Generated by analysis/skylight_groundtruth.py; do not edit by hand.",
        "\\begin{tabular}{@{}lcccccccc@{}}",
        "\\toprule",
        " & \\multicolumn{2}{c}{Tape (m)} & \\multicolumn{4}{c}{Lidar roof thickness at the rim (m)} & \\multicolumn{2}{c}{Rim to floor (m)} \\\\",
        "\\cmidrule(lr){2-3}\\cmidrule(lr){4-7}\\cmidrule(lr){8-9}",
        "Rim point & lip & terrace & offline map & second flight & onboard map & flanking stations & tape & lidar \\\\",
        "\\midrule"]
    for k in order:
        r = by[k]
        lines.append(" & ".join([label[k], cell(r["tape_upper"], 2), cell(r["tape_lower"], 2),
                                 cell(r["tau_rim_flf"]), cell(r["tau_rim_slf"]), cell(r["tau_rim_tube"]),
                                 cell(r["tau_flank_mean"]), cell(r["tape_floor"], 2),
                                 cell(r["floor_depth_flf"])]) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    (OUT / "skylight_groundtruth_table.tex").write_text("\n".join(lines) + "\n")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42})
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0))
    names = [r["point"] for r in results]; xi = np.arange(len(names))
    ax = axs[0]
    ax.axhspan(-SIGMA_TAU, SIGMA_TAU, color="0.9", lw=0, label=r"$\pm\sigma_\tau$ = 1.5 m")
    ax.axhline(0, color="k", lw=0.6)
    ax.plot(xi, arr("tau_rim_flf") - tape, "o", color="tab:blue", label="Rim strip, offline map")
    ax.plot(xi, arr("tau_rim_tube") - tape, "s", color="tab:cyan", ms=4, label="Rim strip, onboard map")
    ax.plot(xi, arr("tau_rim_slf") - tape, "^", color="tab:purple", ms=4, label="Rim strip, second flight")
    ax.plot(xi, arr("tau_flank_mean") - tape, "D", color="tab:green", ms=4, label="Pipeline, flanking stations")
    ax.set_xticks(xi); ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_ylabel("Lidar roof thickness minus tape (m)")
    ax.legend(fontsize=6, frameon=False)
    ax = axs[1]
    m = np.isfinite(tf)
    ax.plot(xi[m], tf[m], "k_", ms=12, mew=1.5, label="Tape, rim to floor (Sep 2026)")
    ax.plot(xi, arr("floor_depth_flf"), "o", color="tab:blue", label="Lidar floor, offline map (Apr 2025)")
    ax.plot(xi, arr("floor_depth_slf"), "^", color="tab:purple", ms=4, label="Lidar floor, second flight (Apr 2025)")
    # the photogrammetric floor is kept in the JSON only: within 3 m of a rim
    # point the aerial returns below the surface are rim walls, not floor
    ax.set_xticks(xi); ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_ylabel("Rim-to-floor height (m)")
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_skylight_groundtruth.pdf")
    fig.savefig(OUT / "fig_skylight_groundtruth.png", dpi=200)

    # console report
    print(f"{'point':14s} {'ELEV_B':>7s} {'DEM':>7s} {'TIF':>7s} | {'tape':>5s} {'rimFLF':>7s} {'rimDLIO':>7s} {'rimSLF':>7s} {'flank':>6s} | {'tapeF':>5s} {'lidF':>6s} {'photoF':>6s}")
    for r in results:
        g = lambda k: r.get(k, np.nan)
        print(f"{r['point']:14s} {g('elev_birgir'):7.2f} {g('z_dem_cache'):7.2f} {g('z_dem_geotiff'):7.2f} | "
              f"{g('tape_upper'):5.2f} {g('tau_rim_flf'):7.2f} {g('tau_rim_tube'):7.2f} {g('tau_rim_slf'):7.2f} {g('tau_flank_mean'):6.2f} | "
              f"{g('tape_floor'):5.2f} {g('floor_depth_flf'):6.2f} {g('floor_depth_photo'):6.2f}")
    print(json.dumps(openings, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
