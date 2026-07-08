"""S5: hero renders of the merged surface-subsurface model (Fig. 2 material).

Uses the open3d *legacy* Visualizer with an invisible window (the new
OffscreenRenderer needs EGL, unavailable on macOS). Scene: the aerial
photogrammetric surface in true RGB + the tube (10 cm working copy, aerial
frame) coloured by elevation, + the centreline as a tube of small spheres,
+ skylight markers.

Outputs (analysis_out/renders/):
    hero_oblique.png      the money shot: surface with the glowing tube below
    hero_topdown.png      orthographic top view
    hero_side.png         along-tube elevation view
    skylight_closeup.png  camera above S1 looking through the opening

Run:  .venv/bin/python analysis/render_hero.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import open3d as o3d
from matplotlib import cm

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
REN = OUT / "renders"
REN.mkdir(exist_ok=True)

AERIAL = ROOT / "maps/aerial_10cm.pcd"
TUBE = ROOT / "maps/tube_10cm_aerial.pcd"
W, H = 3840, 2160


def load_cloud(path: Path) -> tuple[np.ndarray, np.ndarray | None]:
    """Return xyz float64 and rgb in [0,1] if the cloud has a packed rgb field."""
    xyz_parts, rgb_parts = [], []
    with BinaryPcdReader(path) as r:
        has_rgb = "rgb" in r.fields
        rgb_i = r.fields.index("rgb") if has_rgb else -1
        for chunk in r.chunks():
            xyz_parts.append(chunk[:, :3].astype(np.float64))
            if has_rgb:
                packed = chunk[:, rgb_i].astype(np.float32).view(np.uint32)
                rgb = np.stack([(packed >> 16) & 255, (packed >> 8) & 255,
                                packed & 255], axis=1).astype(np.float64) / 255.0
                rgb_parts.append(rgb)
    xyz = np.vstack(xyz_parts)
    return xyz, (np.vstack(rgb_parts) if rgb_parts else None)


def make_scene() -> list[o3d.geometry.Geometry]:
    geoms = []
    axyz, argb = load_cloud(AERIAL)
    apc = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(axyz))
    if argb is not None:
        apc.colors = o3d.utility.Vector3dVector(argb)
    geoms.append(apc)
    print(f"aerial: {len(axyz):,} pts, rgb={'yes' if argb is not None else 'no'}")

    txyz, _ = load_cloud(TUBE)
    z = txyz[:, 2]
    zlo, zhi = np.percentile(z, 1), np.percentile(z, 99)
    tcol = cm.get_cmap("turbo")((z - zlo) / max(zhi - zlo, 1e-6))[:, :3]
    tpc = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(txyz))
    tpc.colors = o3d.utility.Vector3dVector(tcol)
    geoms.append(tpc)
    print(f"tube: {len(txyz):,} pts, z {zlo:.1f}..{zhi:.1f}")

    cl_path = OUT / "centreline.csv"
    if cl_path.exists():
        rows = list(csv.DictReader(open(cl_path)))
        pts = np.array([[float(r["x"]), float(r["y"]), float(r["z"])]
                        for r in rows])
        ls = o3d.geometry.LineSet(
            o3d.utility.Vector3dVector(pts),
            o3d.utility.Vector2iVector([[i, i + 1] for i in range(len(pts) - 1)]))
        ls.colors = o3d.utility.Vector3dVector(
            np.tile([[1.0, 0.1, 0.1]], (len(pts) - 1, 1)))
        geoms.append(ls)
    return geoms


def capture(geoms, lookat, front, up, zoom, out_png, point_size=1.5):
    vis = o3d.visualization.Visualizer()
    vis.create_window(width=W, height=H, visible=False)
    for g in geoms:
        vis.add_geometry(g)
    opt = vis.get_render_option()
    opt.point_size = point_size
    opt.background_color = np.array([0.02, 0.02, 0.04])
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
    geoms = make_scene()
    sky = json.loads((OUT / "aerial_holes.json").read_text())
    centres = np.array([h["centroid"] for h in sky["holes"][:3]])
    mid = centres.mean(0)

    # Scene centre: middle of the tube footprint.
    tube_xyz = np.asarray(geoms[1].points)
    c = tube_xyz.mean(0)

    capture(geoms, lookat=c, front=[0.55, -0.55, 0.63], up=[0, 0, 1],
            zoom=0.32, out_png=REN / "hero_oblique.png")
    capture(geoms, lookat=c, front=[0.0, 0.0, 1.0], up=[0, 1, 0],
            zoom=0.42, out_png=REN / "hero_topdown.png")
    capture(geoms, lookat=c, front=[0.9, -0.2, 0.38], up=[0, 0, 1],
            zoom=0.25, out_png=REN / "hero_side.png")
    capture(geoms, lookat=mid, front=[0.25, -0.25, 0.93], up=[0, 1, 0],
            zoom=0.06, out_png=REN / "skylight_closeup.png", point_size=2.0)


if __name__ == "__main__":
    main()
