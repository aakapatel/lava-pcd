"""Independent surface check: drone photogrammetry vs national elevation data.

Compares the project DEM (0.5 m, from Birgir's photogrammetric cloud, local
frame with LAVA_PCD_ORIGIN 479158/7089826, hypothesised WGS84 UTM 27N) against
ArcticDEM mosaic v4.1 (2 m, EPSG:3413, ellipsoidal heights) and IslandsDEM
v1.0 (10 m, EPSG:3057, orthometric), both independent of the photogrammetry.

Answers three review items with measurements instead of assumptions:
  - sigma_DEM: the shape agreement (std of the difference after removing a
    constant datum offset, which absorbs geoid and vertical-reference gaps);
  - doming/tilt: plane and quadratic fits to the difference surface over the
    strip (the SfM systematic that would masquerade as an along-tube trend);
  - the CRS hypothesis and any horizontal offset of the photogrammetric
    model, via a small grid search that maximises agreement.

Snow caveat: acquisition dates differ; snow present in either product appears
as a positive/negative bias patch, not removable here, so the std is an upper
bound on the photogrammetric error.

Output: analysis_out_v6/national_dem_check.json, fig_national_dem_check.png.
Run:  PYTHONPATH=src ANALYSIS_OUT=analysis_out_v6 \
      .venv/bin/python analysis/national_dem_check.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v6")))
ORIGIN = np.array([479158.0, 7089826.0])
SRC_EPSG = 32627          # hypothesis: WGS84 / UTM 27N


# DEM cache name follows run_roof's convention (dem_grid_<ROOF_DEM_PCD stem>);
# default unchanged (the v6 merged-cloud cache).
DEM_NPZ = "dem_grid_%s.npz" % Path(os.environ.get(
    "ROOF_DEM_PCD", "full_surface_and_subsurface_merged.pcd")).stem


def load_ours():
    d = np.load(OUT / DEM_NPZ)
    z = d["z"].astype(np.float64)          # x-major: z[ix, iy]
    return z, float(d["xmin"]), float(d["ymin"]), float(d["res"])


def sample_tif(path, xs, ys, epsg_dst):
    tr = Transformer.from_crs(SRC_EPSG, epsg_dst, always_xy=True)
    X, Y = tr.transform(xs, ys)
    with rasterio.open(path) as r:
        vals = np.array([v[0] for v in r.sample(np.c_[X, Y])], dtype=np.float64)
        nod = r.nodata
    if nod is not None:
        vals[vals == nod] = np.nan
    vals[np.abs(vals) > 1e4] = np.nan
    return vals


def analyse(name, path, epsg, z, xmin, ymin, res, step=4):
    ix, iy = np.meshgrid(np.arange(0, z.shape[0], step),
                         np.arange(0, z.shape[1], step), indexing="ij")
    zz = z[ix, iy]
    ok = np.isfinite(zz)
    xs = ORIGIN[0] + xmin + ix[ok] * res
    ys = ORIGIN[1] + ymin + iy[ok] * res
    zloc = zz[ok]

    # horizontal micro-search: shift our coordinates to maximise agreement
    best = None
    for dx in np.arange(-10, 10.1, 2.0):
        for dy in np.arange(-10, 10.1, 2.0):
            v = sample_tif(path, xs + dx, ys + dy, epsg)
            m = np.isfinite(v)
            if m.sum() < 500:
                continue
            d = zloc[m] - v[m]
            sd = float(np.std(d))
            if best is None or sd < best[0]:
                best = (sd, float(dx), float(dy))
    sd0, bdx, bdy = best
    v = sample_tif(path, xs + bdx, ys + bdy, epsg)
    m = np.isfinite(v)
    d = zloc[m] - v[m]
    xr, yr = xs[m] - xs[m].mean(), ys[m] - ys[m].mean()

    # plane (tilt) fit
    A1 = np.c_[xr, yr, np.ones(m.sum())]
    c1, *_ = np.linalg.lstsq(A1, d, rcond=None)
    tilt_m_per_km = float(np.hypot(c1[0], c1[1]) * 1000)
    r1 = d - A1 @ c1
    # quadratic (doming) fit on top
    A2 = np.c_[xr ** 2, yr ** 2, xr * yr, xr, yr, np.ones(m.sum())]
    c2, *_ = np.linalg.lstsq(A2, d, rcond=None)
    quad = A2 @ c2 - (A2[:, 3:] @ c2[3:])
    dome_amp = float(np.percentile(quad, 98) - np.percentile(quad, 2))
    r2 = d - A2 @ c2

    return dict(
        product=name, n=int(m.sum()),
        best_shift_m=[bdx, bdy],
        datum_offset_m=float(np.median(d)),
        std_raw_m=float(np.std(d)),
        tilt_m_per_km=tilt_m_per_km,
        std_after_plane_m=float(np.std(r1)),
        doming_amplitude_m=dome_amp,
        std_after_quadratic_m=float(np.std(r2))), (xs[m], ys[m], d)


def main() -> None:
    z, xmin, ymin, res = load_ours()
    res_a, cloud_a = analyse("ArcticDEM v4.1 2 m (ellipsoidal)",
                             ROOT / "maps/national_dem_site.tif", 3413,
                             z, xmin, ymin, res)
    res_i, cloud_i = analyse("IslandsDEM v1.0 10 m (orthometric)",
                             ROOT / "maps/islandsdem10m_site.tif", 3057,
                             z, xmin, ymin, res)
    out = dict(
        note="Independent surface check (review items E47/R2.9/R3.12): "
             "photogrammetric DEM vs national products over the survey "
             "footprint. Datum offset absorbs geoid/ellipsoid and reference "
             "differences; std and doming measure shape error, upper-bounded "
             "by snow/temporal change. CRS hypothesis for the local frame: "
             f"EPSG:{SRC_EPSG} with origin {ORIGIN.tolist()}.",
        arcticdem=res_a, islandsdem=res_i)
    (OUT / "national_dem_check.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, (xs, ys, d), t in zip(axes, [cloud_a, cloud_i],
                                  ["ArcticDEM 2 m", "IslandsDEM 10 m"]):
        dd = d - np.median(d)
        sc = ax.scatter(xs - ORIGIN[0], ys - ORIGIN[1], c=dd, s=2,
                        cmap="RdBu_r", vmin=-2, vmax=2, lw=0)
        ax.set_aspect("equal")
        ax.set_title(f"ours minus {t} (datum removed)")
        ax.set_xlabel("x (m, local)")
        plt.colorbar(sc, ax=ax, label="dz (m)")
    fig.tight_layout()
    fig.savefig(OUT / "fig_national_dem_check.png", dpi=140)


if __name__ == "__main__":
    main()
