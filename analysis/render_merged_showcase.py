"""Session 2026-07-10: publication renders of the co-registered model built on
the FULL-EXTENT merged cloud (full_surface_and_subsurface_merged.pcd).

Improvements over render_hero.py:
  - surface = merged_surface_10cm.pcd (full extent: the conduit no longer
    leaves the surface footprint), true colour
  - tube = flf_10cm_aerial.pcd (the corrected FAST-LIO basis of the paper)
    coloured by DEPTH BELOW THE LOCAL SURFACE (perceptual colormap), which is
    the quantity the paper is about, instead of raw z with a rainbow
  - white background for print, higher point size, tighter framing

Outputs (analysis_out/renders/):
    showcase_oblique.png   full-extent cutaway hero
    showcase_topdown.png   orthographic plan with cutaway
    showcase_closeup.png   skylight rim, both data sets meeting
    showcase_side.png      long oblique down the tube axis

Run:  DISPLAY=:1 env -u PYTHONPATH ./.venv/bin/python analysis/render_merged_showcase.py
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
import open3d as o3d
import matplotlib.pyplot as _plt
from scipy.spatial import cKDTree

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
REN = OUT / "renders"
REN.mkdir(exist_ok=True)

SURFACE = ROOT / "maps/merged_surface_10cm.pcd"
TUBE = ROOT / "maps/flf_10cm_aerial.pcd"
W, H = 3840, 2160
CMAP = os.environ.get("SHOWCASE_CMAP", "inferno")
DEPTH_MAX = float(os.environ.get("SHOWCASE_DEPTH_MAX", "14"))  # m, colour scale

dem = np.load(OUT / "dem_grid_full_surface_and_subsurface_merged.npz")
DEMZ, XMIN, YMIN, RES = dem["z"], float(dem["xmin"]), float(dem["ymin"]), float(dem["res"])
NY, NX = DEMZ.shape
# fill DEM gaps by nearest valid cell so depth-below-surface never goes NaN
from scipy.ndimage import distance_transform_edt
_nan = ~np.isfinite(DEMZ)
if _nan.any():
    _, idx = distance_transform_edt(_nan, return_indices=True)
    DEMZ = DEMZ[idx[0], idx[1]]


def dem_at(xy: np.ndarray) -> np.ndarray:
    ix = np.clip(((xy[:, 0] - XMIN) / RES).astype(np.int64), 0, NX - 1)
    iy = np.clip(((xy[:, 1] - YMIN) / RES).astype(np.int64), 0, NY - 1)
    return DEMZ[iy, ix]


def load_xyz_rgb(path: Path) -> tuple[np.ndarray, np.ndarray | None]:
    xyz_parts, rgb_parts = [], []
    with BinaryPcdReader(path) as r:
        has_rgb = "rgb" in r.fields
        rgb_i = r.fields.index("rgb") if has_rgb else -1
        for chunk in r.chunks():
            xyz_parts.append(chunk[:, :3].astype(np.float64))
            if has_rgb:
                packed = chunk[:, rgb_i].astype(np.float32).view(np.uint32)
                rgb_parts.append(np.stack([(packed >> 16) & 255, (packed >> 8) & 255,
                                           packed & 255], 1).astype(np.float64) / 255.0)
    return np.vstack(xyz_parts), (np.vstack(rgb_parts) if rgb_parts else None)


def centreline_xy() -> np.ndarray:
    rows = list(csv.DictReader(open(OUT / "centreline.csv")))
    return np.array([[float(r["x"]), float(r["y"])] for r in rows])


def build_geoms(cutaway: float) -> tuple[list, np.ndarray]:
    sxyz, srgb = load_xyz_rgb(SURFACE)
    if cutaway > 0:
        d, _ = cKDTree(centreline_xy()).query(sxyz[:, :2], workers=-1)
        keep = d > cutaway
        sxyz, srgb = sxyz[keep], srgb[keep]
    # gentle brightness lift so the moss terrain prints lighter
    srgb = np.clip(srgb * 1.12 + 0.02, 0, 1)
    spc = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(sxyz))
    spc.colors = o3d.utility.Vector3dVector(srgb)
    print(f"surface: {len(sxyz):,} pts (cutaway {cutaway} m)")

    txyz, _ = load_xyz_rgb(TUBE)
    depth = dem_at(txyz[:, :2]) - txyz[:, 2]
    t = np.clip(depth / DEPTH_MAX, 0, 1)
    # shallow (thin roof) = bright orange, deep = dark purple; avoid the
    # near-white and near-black ends of the map so nothing washes out on print
    tcol = _plt.get_cmap(CMAP)(0.88 - 0.74 * t)[:, :3]
    tpc = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(txyz))
    tpc.colors = o3d.utility.Vector3dVector(tcol)
    print(f"tube: {len(txyz):,} pts, depth p1/p99 = "
          f"{np.percentile(depth,1):.1f}/{np.percentile(depth,99):.1f} m")
    return [spc, tpc], txyz


def capture(geoms, lookat, front, up, zoom, out_png, point_size=2.0):
    vis = o3d.visualization.Visualizer()
    vis.create_window(width=W, height=H, visible=False)
    for g in geoms:
        vis.add_geometry(g)
    opt = vis.get_render_option()
    opt.point_size = point_size
    opt.background_color = np.array([1.0, 1.0, 1.0])
    vc = vis.get_view_control()
    vc.set_lookat(lookat)
    vc.set_front(front)
    vc.set_up(up)
    vc.set_zoom(zoom)
    vis.poll_events()
    vis.update_renderer()
    vis.capture_screen_image(str(out_png), do_render=True)
    vis.destroy_window()
    print("wrote", out_png)


def main() -> None:
    sky = json.loads((OUT / "aerial_holes.json").read_text())
    centres = np.array([h["centroid"] for h in sky["holes"][:3]])

    geoms, txyz = build_geoms(cutaway=8.0)
    c = txyz.mean(0)
    capture(geoms, lookat=c, front=[0.55, -0.55, 0.63], up=[0, 0, 1],
            zoom=0.26, out_png=REN / "showcase_oblique.png")
    capture(geoms, lookat=c, front=[0.0, 0.0, 1.0], up=[0, 1, 0],
            zoom=0.40, out_png=REN / "showcase_topdown.png")
    capture(geoms, lookat=c, front=[0.92, -0.18, 0.35], up=[0, 0, 1],
            zoom=0.22, out_png=REN / "showcase_side.png")

    geoms, _ = build_geoms(cutaway=0.0)
    mid = centres.mean(0)
    capture(geoms, lookat=mid, front=[0.35, -0.35, 0.87], up=[0, 1, 0],
            zoom=0.05, out_png=REN / "showcase_closeup.png", point_size=2.5)


if __name__ == "__main__":
    main()
