"""T5: register the day2 second_long_flight FAST-LIO map to the skylight/aerial
frame (non-destructive; writes slf-prefixed outputs so the committed
first_long_flight S1 registration is untouched).

Reuses the constellation matcher. The FAST-LIO map is z-up (unlike the Mac
first-flight SLAM frame), so ceiling-hole detection uses up=(0,0,1). GICP refine
is skipped (rejected by its own safety net in S1; the constellation solve is the
applied method).

Inputs:  maps/slf_10cm.pcd, analysis_out/aerial_holes.json
Outputs: analysis_out/slf_tube_holes.json, slf_transform_landmark.json,
         slf_registration_report.json, fig_slf_tube_holes.png, fig_slf_match.png

Run:  env -u PYTHONPATH .venv/bin/python analysis/run_registration_slf.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lava_pcd.geometry import transform_points
from lava_pcd.holes import HoleSet, build_occupancy, detect_holes
from lava_pcd.register import match_constellations, vertical_residuals

ROOT = Path(__file__).resolve().parents[1]
MAPS = ROOT / "maps"
OUT = ROOT / "analysis_out"

AERIAL = MAPS / "aerial_crop.pcd"
# Tube cloud selectable so the same z-up ceiling registration can be run on any
# FAST-LIO map (second_long_flight default; pass a path + tag to reuse for the
# first_long_flight FAST-LIO map). Outputs are prefixed by <tag>_.
TUBE = Path(sys.argv[1]) if len(sys.argv) > 1 else (MAPS / "slf_10cm.pcd")
TAG = sys.argv[2] if len(sys.argv) > 2 else "slf"

TUBE_PARAMS = dict(
    mode="ceiling", up=(0.0, 0.0, 1.0), resolution=1.0,
    min_area=4.0, max_area=400.0, edge_margin=3.0, ceiling_jump=2.0,
    merge_factor=1.3,
)
TOLERANCE = 3.0   # slightly looser than the first-flight 2.5 to allow for the
                  # different SLAM frame / rim sampling


def hole_summary(hs):
    return [dict(id=h.id, centroid=[round(c, 3) for c in h.centroid],
                 area=round(h.area, 1), semi_major=round(h.semi_major, 2),
                 semi_minor=round(h.semi_minor, 2)) for h in hs.holes]


def render_holes(pcd, params, holes, out_png):
    grid = build_occupancy(pcd, mode=params["mode"], up=params.get("up"),
                           resolution=params["resolution"])
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.imshow(np.log1p(grid.counts.T), origin="lower", extent=grid.extent,
              cmap="viridis")
    R = np.asarray(grid.R)
    for h in holes.holes:
        c = R @ np.asarray(h.centroid)
        ax.plot(c[0], c[1], "r+", ms=10)
        ax.add_patch(plt.matplotlib.patches.Ellipse(
            (c[0], c[1]), 2 * h.semi_major, 2 * h.semi_minor,
            angle=np.degrees(h.orientation), fill=False, color="r", lw=1.2))
        ax.annotate(str(h.id), (c[0], c[1]), color="w", fontsize=9,
                    xytext=(4, 4), textcoords="offset points")
    ax.set_title(f"{pcd.name}: {len(holes.holes)} ceiling holes")
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)


def render_match(aerial_hs, tube_hs, tf, out_png):
    up = np.asarray(tf.aerial_up, float); up /= np.linalg.norm(up)
    e1 = np.cross(up, [0, 0, 1.0])
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.array([1.0, 0, 0])
    e1 /= np.linalg.norm(e1); e2 = np.cross(up, e1)
    P = lambda X: np.c_[X @ e1, X @ e2]
    A, B = aerial_hs.centroids(), tube_hs.centroids()
    Bt = transform_points(B, np.asarray(tf.matrix))
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, Bx, title in ((axes[0], B, "before"), (axes[1], Bt, "after")):
        pa, pb = P(A), P(Bx)
        ax.scatter(pa[:, 0], pa[:, 1], c="tab:blue", label="aerial", s=40)
        ax.scatter(pb[:, 0], pb[:, 1], c="tab:orange", label="tube", s=40, marker="x")
        for t_i, a_i in tf.inliers:
            ax.plot([pa[a_i, 0], pb[t_i, 0]], [pa[a_i, 1], pb[t_i, 1]],
                    c="tab:green", lw=0.8)
        ax.set_title(f"{title} (RMS {tf.rms:.2f} m)" if title == "after" else title)
        ax.set_aspect("equal"); ax.legend(fontsize=8)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    aerial_hs = HoleSet.from_json(OUT / "aerial_holes.json")
    print(f"[1/3] aerial holes (cached): {len(aerial_hs.holes)}")

    print(f"[2/3] ceiling hole detection on {TUBE.name} ...")
    tube_hs = detect_holes(TUBE, **TUBE_PARAMS)
    tube_hs.to_json(OUT / f"{TAG}_tube_holes.json")
    render_holes(TUBE, TUBE_PARAMS, tube_hs, OUT / f"fig_{TAG}_tube_holes.png")
    print(f"      {len(tube_hs.holes)} holes")

    print("[3/3] constellation match ...")
    tf = match_constellations(aerial_hs, tube_hs, mode="auto", tolerance=TOLERANCE)
    tf.to_json(OUT / f"{TAG}_transform_landmark.json")
    vres = vertical_residuals(aerial_hs, tube_hs, tf)
    A, B = aerial_hs.centroids(), tube_hs.centroids()
    ruler = []
    pairs = list(tf.inliers)
    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            (t1, a1), (t2, a2) = pairs[i], pairs[j]
            da = float(np.linalg.norm(A[a1] - A[a2]))
            dt = float(np.linalg.norm(B[t1] - B[t2]))
            ruler.append(dict(pair=f"a{a1}-a{a2}", aerial_m=round(da, 3),
                              tube_m=round(dt, 3), diff_m=round(dt - da, 3)))
    report = dict(
        inputs=dict(aerial=AERIAL.name, tube=TUBE.name, tag=TAG,
                    tube_source=f"{TAG} FAST-LIO scans.pcd"),
        tube_params={k: (list(v) if isinstance(v, tuple) else v)
                     for k, v in TUBE_PARAMS.items()},
        tolerance=TOLERANCE,
        aerial_holes=hole_summary(aerial_hs),
        tube_holes=hole_summary(tube_hs),
        landmark=dict(mode=tf.mode, matches=len(tf.inliers),
                      rms=round(tf.rms, 3), margin=tf.margin,
                      shape_score=round(tf.shape_score, 3),
                      inliers=[list(p) for p in tf.inliers],
                      vertical_residuals_m=[round(v, 3) for v in vres],
                      warnings=list(tf.warnings)),
        scale_ruler=ruler,
    )
    (OUT / f"{TAG}_registration_report.json").write_text(json.dumps(report, indent=2))
    print(f"      mode={tf.mode} matches={len(tf.inliers)} RMS={tf.rms:.3f} "
          f"margin={tf.margin} vres={[f'{v:+.2f}' for v in vres]}")
    render_match(aerial_hs, tube_hs, tf, OUT / f"fig_{TAG}_match.png")
    print(f"wrote {TAG}_registration_report.json")


if __name__ == "__main__":
    main()
