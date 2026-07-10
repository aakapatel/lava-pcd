"""Split full_surface_and_subsurface_merged.pcd into working copies for the
visualization push (session 2026-07-10).

The merged cloud (186M pts, aerial local frame, rgb) contains the full-extent
photogrammetric surface plus a shallow subsurface component. The paper's tau
basis is the FAST-LIO tube (flf_10cm_aerial.pcd) + the DEM rasterised from this
merged cloud, so for rendering we extract:

  maps/merged_surface_10cm.pcd   surface-only, true colour, 10 cm voxel
  maps/merged_surface_25cm.pcd   coarser copy for interactive/oblique scenes
  maps/merged_sub_10cm.pcd       the merged file's own subsurface component
                                 (diagnostic; the paper tube remains flf)

Surface membership: z >= DEM(x, y) - BAND, with the DEM the cached 0.5 m
per-cell-max grid built from this same cloud (T3). Points in cells without DEM
coverage are kept as surface (edge cells).

Run:  env -u PYTHONPATH ./.venv/bin/python analysis/extract_merged_components.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import open3d as o3d

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
MERGED = ROOT / "maps/full_surface_and_subsurface_merged.pcd"
BAND = 1.2  # m below the per-cell surface still counted as surface

dem = np.load(OUT / "dem_grid_full_surface_and_subsurface_merged.npz")
Z, xmin, ymin, res = dem["z"], float(dem["xmin"]), float(dem["ymin"]), float(dem["res"])
ny, nx = Z.shape


def classify(xyz: np.ndarray) -> np.ndarray:
    """True = surface."""
    ix = ((xyz[:, 0] - xmin) / res).astype(np.int64)
    iy = ((xyz[:, 1] - ymin) / res).astype(np.int64)
    inside = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    out = np.ones(len(xyz), dtype=bool)
    zdem = Z[iy[inside], ix[inside]]
    ok = np.isfinite(zdem)
    sub = inside.copy()
    sub[inside] = ok & (xyz[inside, 2] < zdem - BAND)
    out[sub] = False
    return out


def main() -> None:
    surf_xyz, surf_rgb, sub_xyz, sub_rgb = [], [], [], []
    n_done = 0
    with BinaryPcdReader(MERGED) as r:
        rgb_i = r.fields.index("rgb")
        for chunk in r.chunks():
            xyz = chunk[:, :3]
            packed = chunk[:, rgb_i].astype(np.float32).view(np.uint32)
            rgb = np.stack([(packed >> 16) & 255, (packed >> 8) & 255,
                            packed & 255], axis=1).astype(np.uint8)
            is_surf = classify(xyz.astype(np.float64))
            surf_xyz.append(xyz[is_surf]); surf_rgb.append(rgb[is_surf])
            sub_xyz.append(xyz[~is_surf]); sub_rgb.append(rgb[~is_surf])
            n_done += len(chunk)
            print(f"\r{n_done/1e6:.0f}M pts", end="", flush=True)
    print()

    for name, xs, cs, voxels in (
        ("surface", surf_xyz, surf_rgb, (0.10, 0.25)),
        ("sub", sub_xyz, sub_rgb, (0.10,)),
    ):
        xyz = np.vstack(xs).astype(np.float64)
        rgb = np.vstack(cs).astype(np.float64) / 255.0
        print(f"{name}: {len(xyz):,} pts, z {xyz[:,2].min():.1f}..{xyz[:,2].max():.1f}")
        pc = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(xyz))
        pc.colors = o3d.utility.Vector3dVector(rgb)
        for v in voxels:
            ds = pc.voxel_down_sample(v)
            out = ROOT / f"maps/merged_{name}_{int(v*100)}cm.pcd"
            o3d.io.write_point_cloud(str(out), ds, compressed=False)
            print(f"  wrote {out.name}: {len(ds.points):,} pts @ {v} m")


if __name__ == "__main__":
    main()
