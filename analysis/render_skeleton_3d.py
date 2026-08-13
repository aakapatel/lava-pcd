"""Session 2026-07-10: 3D view of the L1-medial centreline and the section
"ribs" it samples, over the interior shell (extends ED Fig 7).

Scene: FAST-LIO tube shell in faint grey; every 20th station highlighted as a
rib (points in a 0.5 m slab about the section plane) coloured by the section
area; the centreline drawn as a dark red polyline of small spheres.

Output: analysis_out/renders/skeleton_3d.png (+ copy in manuscript figures/)

Run:  DISPLAY=:1 env -u PYTHONPATH ./.venv/bin/python analysis/render_skeleton_3d.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import open3d as o3d
import matplotlib.pyplot as plt

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
import os
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out")))
TUBE_PCD = Path(os.environ.get("SKEL3D_TUBE_PCD", str(ROOT / "maps/flf_10cm_aerial.pcd")))
REN = OUT / "renders"
FIGS = (ROOT.parent / "_Nature__Autonomous_aerial_reconnaissance_of_a_basaltic_"
        "lava_tube_reveals_interior_morphology_and_roof_thickness_distribution_"
        "for_planetary_subsurface" / "figures")
W, H = 3840, 1700
RIB_EVERY = 20  # stations
RIB_HALF = 0.5  # m slab half-width

rows = list(csv.DictReader(open(OUT / "morphometry.csv")))
S = np.array([float(r["s"]) for r in rows])
P = np.array([[float(r["x"]), float(r["y"]), float(r["z"])] for r in rows])
T = np.array([[float(r["tx"]), float(r["ty"]), float(r["tz"])] for r in rows])
A = np.array([float(r["area"]) for r in rows])

xyz_parts = []
with BinaryPcdReader(TUBE_PCD) as r:
    for chunk in r.chunks():
        xyz_parts.append(chunk[:, :3].astype(np.float64))
xyz = np.vstack(xyz_parts)

# cut away the shell above the local centreline so the view looks INTO the
# conduit (otherwise the ceiling hides the ribs and the centreline)
from scipy.spatial import cKDTree
_, near = cKDTree(P[:, :2]).query(xyz[:, :2], workers=-1)
xyz = xyz[xyz[:, 2] < P[near, 2] + 1.2]

# rib assignment: for each highlighted station, points within the slab and
# within 1.6x the station equivalent radius laterally
col = np.full((len(xyz), 3), 0.82)
rib_ids = list(range(3, len(rows) - 3, RIB_EVERY))
amin, amax = np.percentile(A, 2), np.percentile(A, 98)
cmap = plt.get_cmap("viridis")
for i in rib_ids:
    d = xyz - P[i]
    along = d @ T[i]
    lat = np.linalg.norm(d - np.outer(along, T[i]), axis=1)
    rmax = 1.9 * np.sqrt(max(A[i], 1.0) / np.pi) + 2.5
    m = (np.abs(along) < RIB_HALF) & (lat < rmax)
    c = cmap(np.clip((A[i] - amin) / (amax - amin), 0, 1))[:3]
    col[m] = c

pc = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(xyz))
pc.colors = o3d.utility.Vector3dVector(col)

# centreline as a chain of small spheres (visible at any zoom)
cl = []
for i in range(0, len(P), 2):
    sp = o3d.geometry.TriangleMesh.create_sphere(radius=0.35, resolution=8)
    sp.translate(P[i])
    sp.paint_uniform_color([0.75, 0.05, 0.05])
    sp.compute_vertex_normals()
    cl.append(sp)


def capture(front, up, zoom, out_png, point_size=2.2):
    vis = o3d.visualization.Visualizer()
    vis.create_window(width=W, height=H, visible=False)
    vis.add_geometry(pc)
    for g in cl:
        vis.add_geometry(g)
    opt = vis.get_render_option()
    opt.point_size = point_size
    opt.background_color = np.array([1.0, 1.0, 1.0])
    vc = vis.get_view_control()
    vc.set_lookat((P.min(0) + P.max(0)) / 2)
    vc.set_front(front)
    vc.set_up(up)
    vc.set_zoom(zoom)
    vis.poll_events()
    vis.update_renderer()
    vis.capture_screen_image(str(out_png), do_render=True)
    vis.destroy_window()
    print("wrote", out_png)


capture(front=[0.42, -0.50, 0.76], up=[0, 0, 1], zoom=0.46,
        out_png=REN / "skeleton_3d.png")

# trim white margins before handing to the manuscript
from PIL import Image
im = np.asarray(Image.open(REN / "skeleton_3d.png"))
nz_r = np.where((im[..., :3].min(-1) < 245).any(1))[0]
nz_c = np.where((im[..., :3].min(-1) < 245).any(0))[0]
pad = 20
crop = im[max(nz_r[0] - pad, 0):nz_r[-1] + pad, max(nz_c[0] - pad, 0):nz_c[-1] + pad]
Image.fromarray(crop).save(FIGS / "skeleton_3d.png")
print("copied trimmed copy to manuscript figures/")
