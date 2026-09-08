"""v9 datum transfer, step D.2: bring Birgir's re-delivered orthometric point
cloud into the analysis local frame.

Input (2026-09-08 delivery): 20250828_Raufarholshellir_Pointcloud_ISN16.laz,
EPSG:8088 (ISN2016 / Lambert 2016), heights orthometric via the ISH2004
geoid, LAS 1.2 point format 2 (XYZ, intensity, RGB), 590,319,423 points.

X/Y are reprojected EPSG:8088 -> EPSG:32627 with pyproj (the datum step
"ISN2016 to WGS 84 (1)" is a null transformation; the pipeline pyproj
reports is recorded in the JSON), the repo's local origin (479158, 7089826)
is subtracted, and Z is left orthometric. The file is streamed in chunks,
never loaded whole. One pass produces:

    maps/aerial_isn16_ortho_crop.pcd   full resolution, bbox of aerial_crop.pcd
    maps/aerial_isn16_ortho_dem.pcd    full resolution, bbox of the v6 DEM
                                       source (full_surface_and_subsurface_merged)
    maps/aerial_isn16_ortho_10cm.pcd   10 cm voxel centroids over the union bbox

All carry the # LAVA_PCD_ORIGIN header of the repo's writer and the packed
rgb field. Header facts, the pyproj transformer description and the point
counts go to ANALYSIS_OUT/surface_products_isn16.json.

Run:  ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
      .venv/bin/python analysis/convert_isn16_surface.py
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import laspy
import numpy as np
import pyproj

from lava_pcd.io.pcd_writer import BinaryPcdWriter, pack_columns
from lava_pcd.voxel import VoxelDownsampler

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
LAZ = Path(os.environ.get(
    "ISN16_LAZ",
    str(ROOT.parent / "Birgir_data_and_papers_08_sep"
        / "20250828_Raufarholshellir_Pointcloud_ISN16.laz")))
ORIGIN = (479158.0, 7089826.0, 0.0)
TARGET_EPSG = "EPSG:32627"
CHUNK = 5_000_000
VOXEL = 0.10

# Local-frame bboxes (x0, x1, y0, y1). aerial_crop.pcd measured bounds and the
# v6 DEM cache extent (dem_grid_full_surface_and_subsurface_merged.npz).
BBOX_CROP = (1299.11, 1598.40, 409.98, 721.56)
BBOX_DEM = (1200.7193603515625, 1520.7193603515625,
            500.17303466796875, 881.6730346679688)
BBOX_UNION = (min(BBOX_CROP[0], BBOX_DEM[0]), max(BBOX_CROP[1], BBOX_DEM[1]),
              min(BBOX_CROP[2], BBOX_DEM[2]), max(BBOX_CROP[3], BBOX_DEM[3]))

COLS = ("x", "y", "z", "r", "g", "b")
FIELDS = ("x", "y", "z", "rgb")


def in_box(x, y, b):
    return (x >= b[0]) & (x <= b[1]) & (y >= b[2]) & (y <= b[3])


def main() -> None:
    OUT.mkdir(exist_ok=True)
    t0 = time.time()
    with laspy.open(LAZ, laz_backend=laspy.LazBackend.Lazrs) as f:
        h = f.header
        crs = h.parse_crs()
        tr = pyproj.Transformer.from_crs(crs, TARGET_EPSG, always_xy=True)
        info = dict(
            file=str(LAZ), las_version=str(h.version),
            point_format=int(h.point_format.id), point_count=int(h.point_count),
            scales=[float(v) for v in h.scales],
            offsets=[float(v) for v in h.offsets],
            mins=[float(v) for v in h.mins], maxs=[float(v) for v in h.maxs],
            crs_epsg=crs.to_epsg(), crs_name=crs.name,
            reprojection=dict(target=TARGET_EPSG,
                              pyproj_description=tr.description,
                              pyproj_definition=tr.definition,
                              pyproj_accuracy_m=tr.accuracy,
                              pyproj_version=pyproj.__version__,
                              proj_version=pyproj.proj_version_str),
            local_origin=list(ORIGIN),
            bbox_crop_local=list(BBOX_CROP), bbox_dem_local=list(BBOX_DEM),
            bbox_union_local=list(BBOX_UNION), voxel_m=VOXEL,
        )
        print(json.dumps(info, indent=2), flush=True)
        n_total = int(h.point_count)
        vox = VoxelDownsampler(VOXEL)
        n_crop = n_dem = n_union = 0
        zmin = np.inf
        zmax = -np.inf
        mins, maxs = h.mins, h.maxs
        w_crop = BinaryPcdWriter(ROOT / "maps/aerial_isn16_ortho_crop.pcd",
                                 max_points=n_total, fields=FIELDS, origin=ORIGIN)
        w_dem = BinaryPcdWriter(ROOT / "maps/aerial_isn16_ortho_dem.pcd",
                                max_points=n_total, fields=FIELDS, origin=ORIGIN)
        with w_crop, w_dem:
            done = 0
            for pts in f.chunk_iterator(CHUNK):
                gx = np.asarray(pts.x, dtype=np.float64)
                gy = np.asarray(pts.y, dtype=np.float64)
                gz = np.asarray(pts.z, dtype=np.float64)
                ok = ((gx >= mins[0] - 1e-3) & (gx <= maxs[0] + 1e-3)
                      & (gy >= mins[1] - 1e-3) & (gy <= maxs[1] + 1e-3)
                      & (gz >= mins[2] - 1e-3) & (gz <= maxs[2] + 1e-3))
                gx, gy, gz = gx[ok], gy[ok], gz[ok]
                ux, uy = tr.transform(gx, gy)
                lx, ly = ux - ORIGIN[0], uy - ORIGIN[1]
                mu = in_box(lx, ly, BBOX_UNION)
                done += len(pts)
                if mu.any():
                    idx = np.where(ok)[0][mu]
                    r = np.asarray(pts.red)[idx] / 257.0
                    g = np.asarray(pts.green)[idx] / 257.0
                    b = np.asarray(pts.blue)[idx] / 257.0
                    arr = np.column_stack([lx[mu], ly[mu], gz[mu], r, g, b])
                    zmin = min(zmin, float(gz[mu].min()))
                    zmax = max(zmax, float(gz[mu].max()))
                    n_union += int(mu.sum())
                    vox.add(arr)
                    mc = in_box(arr[:, 0], arr[:, 1], BBOX_CROP)
                    if mc.any():
                        w_crop.write_chunk(pack_columns(arr[mc], COLS))
                        n_crop += int(mc.sum())
                    md = in_box(arr[:, 0], arr[:, 1], BBOX_DEM)
                    if md.any():
                        w_dem.write_chunk(pack_columns(arr[md], COLS))
                        n_dem += int(md.sum())
                print(f"  {done/1e6:7.1f} M read, union {n_union/1e6:6.1f} M, "
                      f"crop {n_crop/1e6:5.1f} M, dem {n_dem/1e6:5.1f} M, "
                      f"{time.time()-t0:6.0f} s", flush=True)
        pts10 = vox.result()
        with BinaryPcdWriter(ROOT / "maps/aerial_isn16_ortho_10cm.pcd",
                             max_points=len(pts10), fields=FIELDS,
                             origin=ORIGIN) as w:
            w.write_chunk(pack_columns(pts10, COLS))
    info.update(dict(
        points_read=done, points_in_union_bbox=n_union,
        points_crop=n_crop, points_dem=n_dem, points_10cm=int(len(pts10)),
        z_range_union_local=[zmin, zmax], seconds=round(time.time() - t0, 1)))
    (OUT / "surface_products_isn16.json").write_text(json.dumps(info, indent=2))
    print(json.dumps({k: info[k] for k in ("points_read", "points_in_union_bbox",
                                            "points_crop", "points_dem",
                                            "points_10cm", "z_range_union_local",
                                            "seconds")}, indent=2))


if __name__ == "__main__":
    main()
