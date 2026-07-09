"""T3 validation: is the merged cloud's surface in the same frame as the
existing aerial_crop DEM / centreline?

Builds a 0.5 m surface DEM from full_surface_and_subsurface_merged.pcd (per-cell
max z, same recipe as run_roof.build_dem) and compares z_DEM(merged) vs
z_DEM(aerial_crop) at the centreline stations that BOTH cover. If they agree to
< ~1 m, the merged surface is co-registered and can extend tau coverage to the
170 stations aerial_crop misses. If they disagree systematically, the merged
cloud is in a different frame and must not be sampled at the existing centreline
(need Aakash's separate surface cloud + transform).

Writes analysis_out/dem_grid_merged.npz and prints the comparison.
Run:  env -u PYTHONPATH .venv/bin/python analysis/merged_dem_validate.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
MERGED = ROOT / "maps/full_surface_and_subsurface_merged.pcd"
DEM_RES = 0.5


def build_dem(path: Path) -> dict:
    xmin = ymin = np.inf
    xmax = ymax = -np.inf
    with BinaryPcdReader(path) as r:
        for c in r.chunks():
            xmin = min(xmin, c[:, 0].min()); xmax = max(xmax, c[:, 0].max())
            ymin = min(ymin, c[:, 1].min()); ymax = max(ymax, c[:, 1].max())
    nx = int(np.ceil((xmax - xmin) / DEM_RES)) + 1
    ny = int(np.ceil((ymax - ymin) / DEM_RES)) + 1
    zmax = np.full((nx, ny), -np.inf, dtype=np.float32)
    with BinaryPcdReader(path) as r:
        for c in r.chunks():
            ix = ((c[:, 0] - xmin) / DEM_RES).astype(np.int64).clip(0, nx - 1)
            iy = ((c[:, 1] - ymin) / DEM_RES).astype(np.int64).clip(0, ny - 1)
            np.maximum.at(zmax, (ix, iy), c[:, 2])
    filled = zmax > -np.inf
    _, (ii, jj) = ndi.distance_transform_edt(~filled, return_indices=True)
    zfill = zmax[ii, jj]
    zmed = ndi.median_filter(zfill, size=3, mode="nearest")
    z = np.where(filled, zmed, np.nan).astype(np.float32)
    return dict(z=z, xmin=float(xmin), ymin=float(ymin), res=DEM_RES,
                filled=filled)


def dem_at(dem, x, y):
    ix = ((x - dem["xmin"]) / dem["res"]).astype(np.int64)
    iy = ((y - dem["ymin"]) / dem["res"]).astype(np.int64)
    z = dem["z"]
    ok = (ix >= 0) & (ix < z.shape[0]) & (iy >= 0) & (iy < z.shape[1])
    out = np.full(len(x), np.nan)
    out[ok] = z[ix[ok], iy[ok]]
    return out


def main():
    rows = list(csv.DictReader(open(OUT / "morphometry.csv")))
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    zc = np.array([float(r["z_ceil"]) if r["z_ceil"] else np.nan for r in rows])

    print("building merged 0.5 m DEM (streamed) ...")
    dm = build_dem(MERGED)
    # Same cache name run_roof.py uses when ROOF_DEM_PCD points at the merged
    # cloud, so the validation and the roof run share one DEM.
    np.savez_compressed(OUT / f"dem_grid_{MERGED.stem}.npz",
                        z=dm["z"], xmin=dm["xmin"], ymin=dm["ymin"],
                        res=dm["res"])

    # aerial_crop DEM (existing cache used by run_roof)
    d = np.load(OUT / "dem_grid.npz")
    da = dict(z=d["z"], xmin=float(d["xmin"]), ymin=float(d["ymin"]),
              res=float(d["res"]))

    zdem_m = dem_at(dm, x, y)
    zdem_a = dem_at(da, x, y)

    # Only compare where the merged cell was actually filled (not EDT-inpainted)
    ix = ((x - dm["xmin"]) / dm["res"]).astype(np.int64)
    iy = ((y - dm["ymin"]) / dm["res"]).astype(np.int64)
    inb = (ix >= 0) & (ix < dm["z"].shape[0]) & (iy >= 0) & (iy < dm["z"].shape[1])
    real = np.zeros(len(rows), bool)
    real[inb] = dm["filled"][ix[inb], iy[inb]]

    both = ~np.isnan(zdem_m) & ~np.isnan(zdem_a)
    diff = zdem_m[both] - zdem_a[both]
    tau_m = zdem_m - zc
    tau_a = zdem_a - zc

    def stats(a):
        a = a[~np.isnan(a)]
        return dict(n=int(len(a)), min=round(float(a.min()), 2),
                    max=round(float(a.max()), 2),
                    median=round(float(np.median(a)), 2),
                    mean=round(float(a.mean()), 2))

    res = dict(
        n_stations=len(rows),
        n_merged_dem=int((~np.isnan(zdem_m)).sum()),
        n_merged_dem_realcell=int(real.sum()),
        n_aerial_dem=int((~np.isnan(zdem_a)).sum()),
        overlap_n=int(both.sum()),
        dem_diff_merged_minus_aerial=dict(
            median=round(float(np.median(diff)), 2),
            mean=round(float(diff.mean()), 2),
            std=round(float(diff.std()), 2),
            p10=round(float(np.percentile(diff, 10)), 2),
            p90=round(float(np.percentile(diff, 90)), 2),
            absmedian=round(float(np.median(np.abs(diff))), 2)),
        tau_merged=stats(tau_m[real]),
        tau_aerial=stats(tau_a),
    )
    print(json.dumps(res, indent=2))
    (OUT / "merged_dem_validate.json").write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
