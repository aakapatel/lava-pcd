"""v9 datum transfer, step D.3: measure the old-to-new surface transform.

Compares the ellipsoidal crop used since June (maps/aerial_crop.pcd, WGS84
heights, X/Y reprojected from lon/lat) with the same footprint cut from
Birgir's orthometric re-delivery (maps/aerial_isn16_ortho_crop.pcd, ISH2004
heights, X/Y reprojected from ISN2016 Lambert). Both are in the repo's local
frame (EPSG:32627 minus 479158/7089826).

Measurements (all written to ANALYSIS_OUT/datum_transfer.json):
  1. Nearest-neighbour vertical differences (new minus old) over the tube
     footprint: for a random sample of old points within FOOT_M of the v6
     centreline, the nearest new point in 2D within NN_MAX_M. Mean, median,
     std, percentiles, and a plane fit dz = a x + b y + c.
  2. Horizontal shift: per-cell-maximum DEMs from both clouds on the SAME
     0.5 m grid (the v6 cache grid origin), phase cross-correlation with a
     parabolic sub-cell peak fit, at 0.5 m and again at 0.25 m. Then the
     cell-wise difference of the two DEMs (this isolates the effect of the
     16% point-set difference from any grid-phase effect).
  3. Point density of both crops (points per m2 over the shared bbox).
  4. Per-cell-maximum DEM from the new crop against Birgir's 6.16 cm GeoTIFF
     (Metashape's interpolated DEM) sampled at the same cells: sign and size
     of the per-cell-max bias on the rough flow top.

Nothing here is used by the pipeline; it decides whether the plan's rule
(dz std < 5 cm, |dx|,|dy| < 0.5 m) allows a pure translation.

Run:  ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
      .venv/bin/python analysis/datum_transfer_check.py
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
REF = ROOT / "analysis_out_v6"
OLD = ROOT / "maps/aerial_crop.pcd"
NEW = ROOT / "maps/aerial_isn16_ortho_crop.pcd"
TIF = Path(os.environ.get(
    "ISN16_DEM_TIF",
    str(ROOT.parent / "Birgir_data_and_papers_08_sep"
        / "20250828_Raufarholshellir_DEM_ISN16.tif")))
ORIGIN = np.array([479158.0, 7089826.0])
FOOT_M = 40.0        # footprint half-width about the centreline
NN_MAX_M = 0.05      # 2D nearest-neighbour acceptance radius
N_SAMPLE = 400_000
RNG = np.random.default_rng(3)


def load_xyz(p: Path) -> np.ndarray:
    parts = []
    with BinaryPcdReader(p) as r:
        for c in r.chunks():
            parts.append(c[:, :3].astype(np.float64))
    return np.vstack(parts)


def raster_max(pts, xmin, ymin, res, nx, ny):
    ix = ((pts[:, 0] - xmin) / res).astype(np.int64)
    iy = ((pts[:, 1] - ymin) / res).astype(np.int64)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    z = np.full((nx, ny), -np.inf, dtype=np.float64)
    np.maximum.at(z, (ix[ok], iy[ok]), pts[ok, 2])
    z[~np.isfinite(z)] = np.nan
    return z


def phase_shift(a, b, res):
    """Sub-cell shift of b relative to a (b = a moved by +shift) by FFT
    cross-correlation of the mean-removed, NaN-filled rasters."""
    m = np.isfinite(a) & np.isfinite(b)
    aa = np.where(m, a - np.nanmean(a[m]), 0.0)
    bb = np.where(m, b - np.nanmean(b[m]), 0.0)
    F = np.fft.fft2(aa) * np.conj(np.fft.fft2(bb))
    cc = np.real(np.fft.ifft2(F))
    cc = np.fft.fftshift(cc)
    i, j = np.unravel_index(np.argmax(cc), cc.shape)
    ci, cj = cc.shape[0] // 2, cc.shape[1] // 2

    def para(y0, y1, y2):
        d = y0 - 2 * y1 + y2
        return 0.0 if d == 0 else 0.5 * (y0 - y2) / d

    di = para(cc[i - 1, j], cc[i, j], cc[i + 1, j]) if 0 < i < cc.shape[0] - 1 else 0.0
    dj = para(cc[i, j - 1], cc[i, j], cc[i, j + 1]) if 0 < j < cc.shape[1] - 1 else 0.0
    # positive lag means b sits at larger index than a
    return dict(dx_m=float(-(i - ci + di) * res), dy_m=float(-(j - cj + dj) * res),
                peak_cell=[int(i - ci), int(j - cj)])


def stats(v):
    v = np.asarray(v, dtype=np.float64)
    v = v[np.isfinite(v)]
    return dict(n=int(len(v)), mean=float(v.mean()), median=float(np.median(v)),
                std=float(v.std()), p5=float(np.percentile(v, 5)),
                p95=float(np.percentile(v, 95)), min=float(v.min()),
                max=float(v.max()))


def main() -> None:
    OUT.mkdir(exist_ok=True)
    old = load_xyz(OLD)
    new = load_xyz(NEW)
    bbox = [float(v) for v in (max(old[:, 0].min(), new[:, 0].min()),
                               min(old[:, 0].max(), new[:, 0].max()),
                               max(old[:, 1].min(), new[:, 1].min()),
                               min(old[:, 1].max(), new[:, 1].max()))]
    area = (bbox[1] - bbox[0]) * (bbox[3] - bbox[2])
    density = dict(old_pts=int(len(old)), new_pts=int(len(new)),
                   shared_bbox_local=bbox, shared_bbox_area_m2=float(area),
                   old_per_m2=float(len(old) / area),
                   new_per_m2=float(len(new) / area),
                   ratio_new_over_old=float(len(new) / len(old)))
    print("density", json.dumps(density, indent=1))

    # footprint: within FOOT_M of the v6 centreline
    cl = list(csv.DictReader(open(REF / "centreline.csv")))
    cxy = np.array([[float(r["x"]), float(r["y"])] for r in cl])
    ct = cKDTree(cxy)
    d_old, _ = ct.query(old[:, :2], workers=-1)
    d_new, _ = ct.query(new[:, :2], workers=-1)
    old_f = old[d_old < FOOT_M]
    new_f = new[d_new < FOOT_M]
    print(f"footprint points: old {len(old_f):,}, new {len(new_f):,}")

    # 1. nearest-neighbour vertical differences
    samp = old_f[RNG.choice(len(old_f), min(N_SAMPLE, len(old_f)), replace=False)]
    nt = cKDTree(new_f[:, :2])
    dist, idx = nt.query(samp[:, :2], workers=-1)
    ok = dist < NN_MAX_M
    dz = new_f[idx[ok], 2] - samp[ok, 2]
    xs, ys = samp[ok, 0], samp[ok, 1]
    A = np.c_[xs - xs.mean(), ys - ys.mean(), np.ones(ok.sum())]
    coef, *_ = np.linalg.lstsq(A, dz, rcond=None)
    resid = dz - A @ coef
    nn = dict(n_sampled=int(len(samp)), n_matched=int(ok.sum()),
              match_radius_m=NN_MAX_M, dz_new_minus_old=stats(dz),
              plane_fit=dict(slope_x_m_per_km=float(coef[0] * 1000),
                             slope_y_m_per_km=float(coef[1] * 1000),
                             intercept_at_centroid_m=float(coef[2]),
                             std_after_plane_m=float(resid.std())),
              nn_2d_distance_m=stats(dist[ok]))
    print("nn", json.dumps(nn, indent=1))

    # 2. horizontal shift on a common grid (v6 cache origin, footprint bbox)
    ref = np.load(REF / "dem_grid_full_surface_and_subsurface_merged.npz")
    gx0, gy0 = float(ref["xmin"]), float(ref["ymin"])
    shift = {}
    dem_diff = {}
    for res in (0.5, 0.25):
        # align grid origin to the v6 cache lattice
        xmin = gx0 + np.floor((bbox[0] - gx0) / res) * res
        ymin = gy0 + np.floor((bbox[2] - gy0) / res) * res
        nx = int(np.ceil((bbox[1] - xmin) / res)) + 1
        ny = int(np.ceil((bbox[3] - ymin) / res)) + 1
        zo = raster_max(old_f, xmin, ymin, res, nx, ny)
        zn = raster_max(new_f, xmin, ymin, res, nx, ny)
        m = np.isfinite(zo) & np.isfinite(zn)
        d = zn[m] - zo[m]
        shift[f"grid_{res}m"] = phase_shift(zo, zn, res)
        # brute-force refinement around zero: shift the NEW points by (sx, sy)
        # and re-rasterise; the true horizontal offset minimises std(diff)
        best = None
        grid_search = {}
        for sx in np.arange(-0.5, 0.51, 0.25):
            for sy in np.arange(-0.5, 0.51, 0.25):
                zs = raster_max(new_f - [sx, sy, 0.0], xmin, ymin, res, nx, ny)
                mm = np.isfinite(zo) & np.isfinite(zs)
                sd = float(np.std(zs[mm] - zo[mm]))
                grid_search[f"{sx:+.2f},{sy:+.2f}"] = round(sd, 4)
                if best is None or sd < best[0]:
                    best = (sd, float(sx), float(sy))
        shift[f"grid_{res}m"]["brute_force_best_shift_m"] = [best[1], best[2]]
        shift[f"grid_{res}m"]["brute_force_std_at_best_m"] = best[0]
        shift[f"grid_{res}m"]["brute_force_std_at_zero_m"] = grid_search["+0.00,+0.00"]
        dem_diff[f"grid_{res}m"] = dict(
            n_cells=int(m.sum()), new_minus_old=stats(d),
            frac_abs_gt_0p10=float(np.mean(np.abs(d - np.median(d)) > 0.10)),
            frac_abs_gt_0p30=float(np.mean(np.abs(d - np.median(d)) > 0.30)))
        if res == 0.5:
            zo5, zn5, xmin5, ymin5, nx5, ny5 = zo, zn, xmin, ymin, nx, ny
    print("shift", json.dumps(shift, indent=1))
    print("dem_diff", json.dumps(dem_diff, indent=1))

    # also: DEM difference at the v6 centreline stations (what run_roof sees)
    ix = ((cxy[:, 0] - xmin5) / 0.5).astype(int)
    iy = ((cxy[:, 1] - ymin5) / 0.5).astype(int)
    okc = (ix >= 0) & (ix < nx5) & (iy >= 0) & (iy < ny5)
    dst = zn5[ix[okc], iy[okc]] - zo5[ix[okc], iy[okc]]
    station_diff = stats(dst)
    station_diff["note"] = ("new minus old per-cell-max DEM at the v6 centreline "
                            "stations inside the aerial_crop bbox, same grid")

    # 4. per-cell-max (new crop) vs Birgir's GeoTIFF at the same cells
    tr = Transformer.from_crs(32627, 8088, always_xy=True)
    ii, jj = np.meshgrid(np.arange(nx5), np.arange(ny5), indexing="ij")
    cx = xmin5 + (ii.ravel() + 0.5) * 0.5 + ORIGIN[0]
    cy = ymin5 + (jj.ravel() + 0.5) * 0.5 + ORIGIN[1]
    X, Y = tr.transform(cx, cy)
    with rasterio.open(TIF) as r:
        rows, cols = rasterio.transform.rowcol(r.transform, X, Y)
        rows, cols = np.asarray(rows), np.asarray(cols)
        r0, r1 = int(rows.min()), int(rows.max()) + 1
        c0, c1 = int(cols.min()), int(cols.max()) + 1
        win = rasterio.windows.Window(c0, r0, c1 - c0, r1 - r0)
        tile = r.read(1, window=win).astype(np.float64)
        nod = r.nodata
        tif_info = dict(crs=str(r.crs), res=[float(r.res[0]), float(r.res[1])],
                        nodata=nod, dtype=str(r.dtypes[0]))
    tile[tile == nod] = np.nan
    v = tile[rows - r0, cols - c0].reshape(nx5, ny5)
    m = np.isfinite(v) & np.isfinite(zn5)
    d = zn5[m] - v[m]
    # slope from the GeoTIFF itself (cell-to-cell), to relate the bias to
    # roughness: bias should grow with local slope for a per-cell maximum
    gy_, gx_ = np.gradient(np.where(np.isfinite(v), v, np.nan), 0.5)
    slope = np.hypot(gx_, gy_)
    ms = m & np.isfinite(slope)
    q = np.nanpercentile(slope[ms], [25, 50, 75])
    by_slope = {}
    for lo, hi, name in ((0, q[0], "slope_q0_q25"), (q[0], q[1], "slope_q25_q50"),
                         (q[1], q[2], "slope_q50_q75"), (q[2], np.inf, "slope_q75_max")):
        sel = ms & (slope >= lo) & (slope < hi)
        by_slope[name] = dict(n=int(sel.sum()),
                              median_bias=float(np.median(zn5[sel] - v[sel])))
    # same at the station cells that run_roof uses
    dstat = zn5[ix[okc], iy[okc]] - v[ix[okc], iy[okc]]
    tif_cmp = dict(geotiff=tif_info, n_cells=int(m.sum()),
                   percell_max_minus_geotiff=stats(d),
                   by_local_slope=by_slope,
                   at_v6_stations=stats(dstat),
                   note="positive = our 0.5 m per-cell maximum sits above "
                        "Metashape's interpolated 6 cm DEM at the same cell")
    print("tif", json.dumps(tif_cmp, indent=1))

    out = dict(
        inputs=dict(old=str(OLD), new=str(NEW), footprint_halfwidth_m=FOOT_M,
                    centreline=str(REF / "centreline.csv")),
        density=density, nearest_neighbour=nn, horizontal_shift=shift,
        dem_difference_same_grid=dem_diff,
        dem_difference_at_v6_stations=station_diff,
        new_percell_max_vs_geotiff=tif_cmp,
    )
    (OUT / "datum_transfer.json").write_text(json.dumps(out, indent=2))

    # figure: dz map + histogram
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    sc = axes[0].scatter(xs, ys, c=dz - np.median(dz), s=1, cmap="RdBu_r",
                         vmin=-0.2, vmax=0.2, lw=0)
    axes[0].set_aspect("equal"); axes[0].set_title("NN dz minus median (m)")
    plt.colorbar(sc, ax=axes[0])
    axes[1].hist(dz, bins=200, color="0.3")
    axes[1].set_title(f"NN dz new-old: median {np.median(dz):.3f} std {dz.std():.3f}")
    dd = zn5 - zo5
    im = axes[2].imshow(dd.T - np.nanmedian(dd), origin="lower", cmap="RdBu_r",
                        vmin=-0.3, vmax=0.3,
                        extent=[xmin5, xmin5 + nx5 * 0.5, ymin5, ymin5 + ny5 * 0.5])
    axes[2].set_title("per-cell-max DEM new-old, same grid (m)")
    plt.colorbar(im, ax=axes[2])
    fig.tight_layout()
    fig.savefig(OUT / "fig_datum_transfer.png", dpi=130)
    print("wrote", OUT / "datum_transfer.json")


if __name__ == "__main__":
    main()
