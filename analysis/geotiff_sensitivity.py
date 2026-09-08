"""v9 datum transfer, step D.7 (sensitivity row only, not the primary basis).

Recomputes the per-station overburden with Birgir's Metashape GeoTIFF DEM
(6.16 cm, EPSG:8088, ISH2004 heights) sampled at the station coordinates
instead of our 0.5 m per-cell-maximum rasterisation of the point cloud:

    tau_tif = z_tif(x, y) - z_ceil        vs     tau (roof_thickness.csv)

Stations are reprojected local -> EPSG:32627 -> EPSG:8088; the GeoTIFF is
read once as a window covering the stations. Reports the per-station
difference distribution (tau_tif minus tau) over the intact stations and the
effect on the headline statistics, appended to ANALYSIS_OUT/datum_transfer.json
under "geotiff_sensitivity". The primary basis is not changed.

Run:  ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
      .venv/bin/python analysis/geotiff_sensitivity.py
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
TIF = Path(os.environ.get(
    "ISN16_DEM_TIF",
    str(ROOT.parent / "Birgir_data_and_papers_08_sep"
        / "20250828_Raufarholshellir_DEM_ISN16.tif")))
ORIGIN = np.array([479158.0, 7089826.0])


def stats(v):
    v = np.asarray(v, dtype=np.float64)
    v = v[np.isfinite(v)]
    return dict(n=int(len(v)), mean=float(v.mean()), median=float(np.median(v)),
                std=float(v.std()), p5=float(np.percentile(v, 5)),
                p95=float(np.percentile(v, 95)), min=float(v.min()),
                max=float(v.max()))


def main() -> None:
    rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    g = lambda k: np.array([float(r[k]) if r[k] else np.nan for r in rows])
    s, x, y, zdem, zceil, tau = g("s"), g("x"), g("y"), g("z_dem"), g("z_ceil"), g("tau")
    klass = np.array([r["class"] for r in rows])
    tr = Transformer.from_crs(32627, 8088, always_xy=True)
    X, Y = tr.transform(x + ORIGIN[0], y + ORIGIN[1])
    with rasterio.open(TIF) as r:
        rr, cc = rasterio.transform.rowcol(r.transform, X, Y)
        rr, cc = np.asarray(rr), np.asarray(cc)
        r0, c0 = int(rr.min()) - 1, int(cc.min()) - 1
        win = rasterio.windows.Window(c0, r0, int(cc.max()) - c0 + 2,
                                      int(rr.max()) - r0 + 2)
        tile = r.read(1, window=win).astype(np.float64)
        tile[tile == r.nodata] = np.nan
        res = float(r.res[0])
    ztif = tile[rr - r0, cc - c0]
    tau_tif = ztif - zceil
    it = klass == "intact"
    d = tau_tif[it] - tau[it]
    dz = ztif[it] - zdem[it]
    with open(OUT / "roof_thickness_geotiff_sensitivity.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["s", "z_dem_percellmax", "z_dem_geotiff", "tau_percellmax",
                    "tau_geotiff", "class"])
        for i in range(len(rows)):
            w.writerow([f"{s[i]:.1f}", f"{zdem[i]:.3f}", f"{ztif[i]:.3f}",
                        f"{tau[i]:.3f}", f"{tau_tif[i]:.3f}", klass[i]])
    out = dict(
        geotiff=str(TIF), geotiff_res_m=res, n_intact=int(it.sum()),
        tau_geotiff_minus_tau_percellmax=stats(d),
        z_geotiff_minus_z_percellmax=stats(dz),
        headline=dict(
            median_tau_percellmax=float(np.nanmedian(tau[it])),
            median_tau_geotiff=float(np.nanmedian(tau_tif[it])),
            max_tau_percellmax=float(np.nanmax(tau[it])),
            max_tau_geotiff=float(np.nanmax(tau_tif[it])),
            q25_q75_percellmax=[float(v) for v in np.nanpercentile(tau[it], [25, 75])],
            q25_q75_geotiff=[float(v) for v in np.nanpercentile(tau_tif[it], [25, 75])],
            n_geotiff_nonpositive=int(np.sum(tau_tif[it] <= 0))),
        note="sensitivity only: Metashape's interpolated DEM sampled at the "
             "station against our 0.5 m per-cell maximum; negative difference "
             "= the per-cell maximum sits higher than the interpolated surface",
    )
    p = OUT / "datum_transfer.json"
    dt = json.loads(p.read_text())
    dt["geotiff_sensitivity"] = out
    p.write_text(json.dumps(dt, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
