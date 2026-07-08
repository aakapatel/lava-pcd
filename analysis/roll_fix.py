"""S1.6: fix the weakly constrained roll about the near-collinear skylight axis.

Three skylights whose centroids are nearly collinear leave one rotation
unconstrained in the 6-DOF constellation solve: roll about the line through
them (the solver warned about exactly this). Any roll keeps the three rim
centroids fixed, so rim residuals cannot detect it, but it swings the rest of
the conduit above or below the terrain.

Two physical facts close the gap, both measured by the data already in hand:
  1. The conduit must lie beneath the surface model everywhere (no part of the
     interior may protrude above the terrain).
  2. The terrain patches that the onboard lidar mapped *through* the three
     apertures coincide with the same terrain in the surface model.

We rotate the registered tube cloud about the best-fit skylight axis by theta,
and minimise  cost(theta) = protrusion(theta) + patch_misfit(theta), where
protrusion is the mean height of tube points above DEM + margin, and
patch_misfit is the trimmed median |z - DEM| of tube points inside the
skylight columns near surface level. The optimum is written as
transform_final.json and the tube working copies are re-exported.

Run:  .venv/bin/python analysis/roll_fix.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lava_pcd.geometry import transform_points
from lava_pcd.holes import HoleSet
from lava_pcd.io.pcd_reader import BinaryPcdReader
from lava_pcd.register import Transform, vertical_residuals

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
TUBE = ROOT / "maps/tube_10cm.pcd"      # SLAM frame (untransformed)
MARGIN = 0.5                             # protrusion allowance (m)
PATCH_Z = 3.0                            # |z - DEM| window for aperture patches
THETAS = np.deg2rad(np.arange(-15.0, 15.01, 0.25))


def load_dem() -> dict:
    d = np.load(OUT / "dem_grid.npz")
    return dict(z=d["z"], xmin=float(d["xmin"]), ymin=float(d["ymin"]),
                res=float(d["res"]))


def dem_at(dem: dict, xy: np.ndarray) -> np.ndarray:
    ix = ((xy[:, 0] - dem["xmin"]) / dem["res"]).astype(np.int64)
    iy = ((xy[:, 1] - dem["ymin"]) / dem["res"]).astype(np.int64)
    z = dem["z"]
    ok = (ix >= 0) & (ix < z.shape[0]) & (iy >= 0) & (iy < z.shape[1])
    out = np.full(len(xy), np.nan)
    out[ok] = z[ix[ok], iy[ok]]
    return out


def rot_about_axis(axis: np.ndarray, theta: float) -> np.ndarray:
    a = axis / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)


def main() -> None:
    tf = Transform.from_json(OUT / "transform_landmark.json")
    M0 = tf.array
    dem = load_dem()
    sky = HoleSet.from_json(OUT / "aerial_holes.json")
    matched_aerial = sorted({a for _, a in tf.inliers})
    anchors = np.asarray(tf.anchors, float)

    # Axis: best-fit line through the three anchors (aerial frame).
    centre = anchors.mean(0)
    _, _, Vt = np.linalg.svd(anchors - centre)
    axis = Vt[0]
    lin_resid = np.linalg.norm(
        (anchors - centre) - np.outer((anchors - centre) @ axis, axis), axis=1)
    print(f"skylight axis {np.round(axis, 3)}, collinearity residuals "
          f"{np.round(lin_resid, 2)} m")

    # Tube cloud in the aerial frame (subsample for the search).
    parts = []
    with BinaryPcdReader(TUBE) as r:
        for chunk in r.chunks():
            parts.append(chunk[:, :3].astype(np.float64))
    P = transform_points(np.vstack(parts), M0)
    rng = np.random.default_rng(3)
    if len(P) > 800_000:
        P = P[rng.choice(len(P), 800_000, replace=False)]

    # Aperture columns for the patch term.
    cols = []
    for a_idx in matched_aerial:
        h = sky.holes[a_idx]
        cols.append((h.centroid[0], h.centroid[1], max(h.semi_major, 1.5) * 1.3))

    def costs(theta: float) -> tuple[float, float]:
        R = rot_about_axis(axis, theta)
        Q = (P - centre) @ R.T + centre
        zd = dem_at(dem, Q[:, :2])
        ok = np.isfinite(zd)
        prot = float(np.mean(np.clip(Q[ok, 2] - zd[ok] - MARGIN, 0, None)))
        in_col = np.zeros(len(Q), bool)
        for cx, cy, rad in cols:
            in_col |= np.hypot(Q[:, 0] - cx, Q[:, 1] - cy) <= rad
        sel = in_col & ok & (np.abs(Q[:, 2] - zd) <= PATCH_Z)
        patch = (float(np.median(np.abs(Q[sel, 2] - zd[sel])))
                 if sel.sum() > 50 else PATCH_Z)
        return prot, patch

    prot = np.empty(len(THETAS))
    patch = np.empty(len(THETAS))
    for i, th in enumerate(THETAS):
        prot[i], patch[i] = costs(th)
    cost = prot + patch
    i_star = int(np.argmin(cost))
    th_star = float(THETAS[i_star])
    print(f"theta* = {np.degrees(th_star):+.2f} deg  "
          f"protrusion {prot[i_star]:.3f} m  patch misfit {patch[i_star]:.3f} m "
          f"(initial: prot {prot[len(THETAS)//2]:.3f}, "
          f"patch {patch[len(THETAS)//2]:.3f})")

    # Compose the final transform: M = T(centre) R T(-centre) M0.
    R = rot_about_axis(axis, th_star)
    A = np.eye(4); A[:3, :3] = R
    A[:3, 3] = centre - R @ centre
    M = A @ M0
    final = Transform(
        matrix=[list(map(float, row)) for row in M],
        inliers=tf.inliers, rms=tf.rms, mode=tf.mode + "+roll",
        tube_up=tf.tube_up, aerial_up=tf.aerial_up,
        anchors=tf.anchors, anchor_axes=tf.anchor_axes,
        warnings=list(tf.warnings) + [
            f"roll about the skylight axis fixed at {np.degrees(th_star):+.2f} "
            f"deg by the terrain-consistency search (roll_fix.py)"],
    )
    final.to_json(OUT / "transform_final.json")

    a_hs = HoleSet.from_json(OUT / "aerial_holes.json")
    t_hs = HoleSet.from_json(OUT / "tube_holes.json")
    vres = vertical_residuals(a_hs, t_hs, final)
    print("vertical residuals after roll fix:", np.round(vres, 3))

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(np.degrees(THETAS), prot, label="protrusion above DEM (m)")
    ax.plot(np.degrees(THETAS), patch, label="aperture patch |dz| (m)")
    ax.plot(np.degrees(THETAS), cost, "k--", lw=1, label="total")
    ax.axvline(np.degrees(th_star), color="r", ls=":", lw=1)
    ax.set_xlabel("roll about skylight axis (deg)")
    ax.set_ylabel("cost (m)")
    ax.legend()
    fig.savefig(OUT / "fig_roll_search.png", dpi=140, bbox_inches="tight")

    diag = dict(theta_deg=round(np.degrees(th_star), 2),
                protrusion_m=round(prot[i_star], 4),
                patch_misfit_m=round(patch[i_star], 4),
                protrusion_before_m=round(float(prot[len(THETAS) // 2]), 4),
                patch_before_m=round(float(patch[len(THETAS) // 2]), 4),
                vertical_residuals_m=[round(float(v), 3) for v in vres],
                axis=[round(float(v), 4) for v in axis],
                collinearity_residuals_m=[round(float(v), 3) for v in lin_resid])
    (OUT / "roll_fix.json").write_text(json.dumps(diag, indent=2))


if __name__ == "__main__":
    main()
