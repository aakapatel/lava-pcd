"""Recover the camera of the archived showcase renders, so vector overlays
(skylight markers, along-tube ticks) can be drawn on them at the right place.

render_merged_showcase.py drives the Open3D legacy visualiser with
set_lookat / set_front / set_up / set_zoom. That is a deterministic pinhole
camera: Open3D puts the eye at lookat + front * (zoom * maxExtent) / tan(fov/2)
with a 60 degree vertical field of view, where maxExtent is the largest side of
the axis-aligned bounding box of the geometry that was added. This script
reconstructs that bounding box from the same two point clouds and the same
8 m cutaway, writes the camera to analysis_out_v10/showcase_camera.json, and
verifies it against the rendered image: the projected tube points must cover
the tube-coloured pixels of the render.

Nothing here re-renders anything; the archived PNGs in
analysis_out_v9/renders/ are the only images used.

Run (repo root):
    ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
        .venv/bin/python analysis/showcase_camera.py
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import binary_dilation
from scipy.spatial import cKDTree

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
OUT10 = ROOT / "analysis_out_v10"
SURFACE = Path(os.environ.get("SHOWCASE_SURFACE_PCD",
                              str(ROOT / "maps/aerial_isn16_ortho_surface_10cm.pcd")))
TUBE = Path(os.environ.get("SHOWCASE_TUBE_PCD",
                           str(ROOT / "maps/flf_10cm_ed_ortho.pcd")))
W, H = 3840, 2160          # render_merged_showcase.py window size
FOV = 60.0                 # Open3D ViewControl FIELD_OF_VIEW_DEFAULT
CUTAWAY = 8.0              # metres, as in render_merged_showcase.build_geoms
VIEWS = {                  # name: (front, up, zoom) from render_merged_showcase
    "oblique":  ([0.55, -0.55, 0.63], [0.0, 0.0, 1.0], 0.26),
    "topdown":  ([0.0, 0.0, 1.0], [0.0, 1.0, 0.0], 0.40),
    "side":     ([0.92, -0.18, 0.35], [0.0, 0.0, 1.0], 0.22),
}
Image.MAX_IMAGE_PIXELS = None


def load_xyz(path: Path) -> np.ndarray:
    parts = []
    with BinaryPcdReader(path) as r:
        for c in r.chunks():
            parts.append(c[:, :3].astype(np.float64))
    return np.vstack(parts)


def camera(bbox_min, bbox_max, lookat, front, up, zoom) -> dict:
    front = np.asarray(front, float)
    front = front / np.linalg.norm(front)
    right = np.cross(np.asarray(up, float), front)
    right /= np.linalg.norm(right)
    up2 = np.cross(front, right)
    up2 /= np.linalg.norm(up2)
    extent = float(np.max(np.asarray(bbox_max) - np.asarray(bbox_min)))
    distance = zoom * extent / np.tan(np.radians(FOV * 0.5))
    eye = np.asarray(lookat, float) + front * distance
    return dict(front=front.tolist(), right=right.tolist(), up=up2.tolist(),
                eye=eye.tolist(), distance=float(distance), extent=extent,
                width=W, height=H, fov_deg=FOV)


def project(cam: dict, P) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """World points -> pixel coordinates of the full-size render."""
    P = np.atleast_2d(np.asarray(P, float))
    d = P - np.asarray(cam["eye"])
    depth = -(d @ np.asarray(cam["front"]))
    t = np.tan(np.radians(cam["fov_deg"] * 0.5))
    aspect = cam["width"] / cam["height"]
    px = ((d @ np.asarray(cam["right"])) / (depth * aspect * t) + 1) * 0.5 * cam["width"]
    py = (1 - (d @ np.asarray(cam["up"])) / (depth * t)) * 0.5 * cam["height"]
    return px, py, depth


def tube_mask(img: np.ndarray) -> np.ndarray:
    """Pixels rendered with the inferno tube ramp (bright orange to dark
    purple); the terrain is desaturated green/tan and the ground is white."""
    a = img.astype(float) / 255.0
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    mx, mn = a.max(-1), a.min(-1)
    sat = np.where(mx > 0, (mx - mn) / np.maximum(mx, 1e-6), 0.0)
    warm = (r > g + 0.18) & (sat > 0.45)
    dark = (mx < 0.45) & (b > g) & (r > g)
    return warm | dark


def main() -> None:
    sxyz = load_xyz(SURFACE)
    txyz = load_xyz(TUBE)
    rows = list(csv.DictReader(open(OUT / "centreline.csv")))
    cxy = np.array([[float(r["x"]), float(r["y"])] for r in rows])
    d, _ = cKDTree(cxy).query(sxyz[:, :2], workers=-1)
    sxyz = sxyz[d > CUTAWAY]
    bmin = np.minimum(sxyz.min(0), txyz.min(0))
    bmax = np.maximum(sxyz.max(0), txyz.max(0))
    lookat = txyz.mean(0)
    print(f"surface {len(sxyz):,} pts after {CUTAWAY:g} m cutaway, "
          f"tube {len(txyz):,} pts")
    print(f"bbox {bmin.round(2).tolist()} .. {bmax.round(2).tolist()}, "
          f"max extent {float(np.max(bmax - bmin)):.2f} m")

    out = dict(source="analysis/showcase_camera.py", surface=str(SURFACE),
               tube=str(TUBE), cutaway_m=CUTAWAY,
               bbox_min=bmin.tolist(), bbox_max=bmax.tolist(),
               lookat=lookat.tolist(), views={})
    sub = txyz[::17]
    for name, (front, up, zoom) in VIEWS.items():
        cam = camera(bmin, bmax, lookat, front, up, zoom)
        png = OUT / "renders" / f"showcase_{name}.png"
        if not png.exists():
            print(f"  {name}: no render, camera written without a check")
            out["views"][name] = cam
            continue
        img = np.asarray(Image.open(png).convert("RGB"))
        assert img.shape[:2] == (H, W), (name, img.shape)
        ref = tube_mask(img)
        px, py, dep = project(cam, sub)
        grid = np.zeros((H, W), bool)
        ok = (dep > 0) & (px >= 0) & (px < W) & (py >= 0) & (py < H)
        grid[py[ok].astype(int), px[ok].astype(int)] = True
        grid = binary_dilation(grid, iterations=4)
        cover = float((grid & ref).sum() / max(ref.sum(), 1))
        cam["check_render"] = png.name
        cam["check_cover_fraction"] = cover
        out["views"][name] = cam
        print(f"  {name}: {cover * 100:.1f}% of the tube-coloured pixels of "
              f"{png.name} are covered by the projected tube cloud")
        if cover < 0.9:
            raise SystemExit(f"{name}: camera model does not reproduce the render")

    OUT10.mkdir(parents=True, exist_ok=True)
    dst = OUT10 / "showcase_camera.json"
    dst.write_text(json.dumps(out, indent=2))
    print("wrote", dst)


if __name__ == "__main__":
    main()
