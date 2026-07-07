"""S1: headless skylight registration + Z-axis analysis (plan/02 stage S1).

Reproduces the full holes -> register -> refine pipeline without any GUI and
writes every intermediate product and metric needed by the manuscript:

    analysis_out/
        aerial_holes.json / tube_holes.json     landmark sets
        transform_landmark.json                 constellation-only (4dof)
        transform_refined_d4.json               + rim GICP, 4-DOF   (default)
        transform_refined_d6.json               + rim GICP, 6-DOF
        registration_report.json                all metrics, all variants
        fig_aerial_holes.png / fig_tube_holes.png   occupancy + detections
        fig_match.png                           constellation before/after

Run:  .venv/bin/python analysis/run_registration.py
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lava_pcd.geometry import transform_points
from lava_pcd.holes import HoleSet, build_occupancy, detect_holes
from lava_pcd.merge import icp_refine
from lava_pcd.register import match_constellations, vertical_residuals

ROOT = Path(__file__).resolve().parents[1]
MAPS = ROOT / "maps"
OUT = ROOT / "analysis_out"
OUT.mkdir(exist_ok=True)

AERIAL = MAPS / "aerial_crop.pcd"
TUBE_FULL = MAPS / "first_long_flight_full.pcd"
TUBE = MAPS / "tube_10cm.pcd"

# Detection parameters. The aerial photogrammetry sees the tube floor *through*
# each skylight, so the hole is not empty in a plain occupancy image; the
# relative-density mode catches the drop against the local ground density.
AERIAL_PARAMS = dict(
    mode="aerial", resolution=1.0, relative=0.5, smooth=0.5,
    min_area=2.0, max_area=400.0,
)
TUBE_PARAMS = dict(
    mode="ceiling", up=(0.0, 0.2, 0.98), resolution=1.0,
    min_area=4.0, max_area=400.0, edge_margin=3.0, ceiling_jump=2.0,
    merge_factor=1.3,
)
TOLERANCE = 2.5
RIM_HEIGHT = 2.0
RIM_INFLATE = 1.5


def hole_summary(hs: HoleSet) -> list[dict]:
    return [
        dict(id=h.id, centroid=[round(c, 3) for c in h.centroid],
             radius=round(h.radius, 2), area=round(h.area, 1),
             semi_major=round(h.semi_major, 2), semi_minor=round(h.semi_minor, 2))
        for h in hs.holes
    ]


def render_holes(pcd: Path, params: dict, holes: HoleSet, out_png: Path) -> None:
    grid = build_occupancy(pcd, mode=params["mode"], up=params.get("up"),
                           resolution=params["resolution"])
    fig, ax = plt.subplots(figsize=(10, 7))
    img = np.log1p(grid.counts.T)
    ax.imshow(img, origin="lower", extent=grid.extent, cmap="viridis")
    R = np.asarray(grid.R)
    for h in holes.holes:
        c = R @ np.asarray(h.centroid)
        ax.plot(c[0], c[1], "r+", ms=10)
        ax.add_patch(plt.matplotlib.patches.Ellipse(
            (c[0], c[1]), 2 * h.semi_major, 2 * h.semi_minor,
            angle=np.degrees(h.orientation), fill=False, color="r", lw=1.2))
        ax.annotate(str(h.id), (c[0], c[1]), color="w", fontsize=9,
                    xytext=(4, 4), textcoords="offset points")
    ax.set_title(f"{pcd.name}: {len(holes.holes)} holes")
    ax.set_xlabel("up-plane a (m)")
    ax.set_ylabel("up-plane b (m)")
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)


def render_match(aerial_hs: HoleSet, tube_hs: HoleSet, tf, out_png: Path) -> None:
    """Before/after constellation plot, projected on the aerial up-plane."""
    up = np.asarray(tf.aerial_up, dtype=float)
    up /= np.linalg.norm(up)
    e1 = np.cross(up, [0.0, 0.0, 1.0])
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.array([1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(up, e1)
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
        ax.set_title(f"{title}  (RMS {tf.rms:.2f} m)" if title == "after" else title)
        ax.set_aspect("equal")
        ax.legend(fontsize=8)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)


def refine_variant(tf, dof: int, tag: str, report: dict, aerial_hs, tube_hs):
    """Run rim GICP at the given DOF and record shift/metrics."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ref = icp_refine(AERIAL, TUBE, tf, dof=dof,
                         rim_height=RIM_HEIGHT, rim_inflate=RIM_INFLATE)
        notes = [str(w.message) for w in caught]
    d = np.asarray(ref.matrix)[:3, 3] - np.asarray(tf.matrix)[:3, 3]
    up = np.asarray(tf.aerial_up, dtype=float)
    up /= np.linalg.norm(up)
    report[f"refined_{tag}"] = dict(
        rms=round(ref.rms, 3),
        translation_shift_m=round(float(np.linalg.norm(d)), 3),
        vertical_shift_m=round(float(d @ up), 3),
        vertical_residuals_m=[round(v, 3) for v in
                              vertical_residuals(aerial_hs, tube_hs, ref)],
        warnings=notes,
    )
    ref.to_json(OUT / f"transform_refined_{tag}.json")
    return ref


def main() -> None:
    report: dict = dict(
        inputs=dict(aerial=str(AERIAL.name), tube=str(TUBE.name),
                    tube_source=str(TUBE_FULL.name)),
        aerial_params=AERIAL_PARAMS, tube_params={
            k: (list(v) if isinstance(v, tuple) else v) for k, v in TUBE_PARAMS.items()},
        tolerance=TOLERANCE, rim_height=RIM_HEIGHT, rim_inflate=RIM_INFLATE,
    )

    print("[1/5] aerial hole detection ...")
    aerial_hs = detect_holes(AERIAL, **AERIAL_PARAMS)
    aerial_hs.to_json(OUT / "aerial_holes.json")
    report["aerial_holes"] = hole_summary(aerial_hs)
    render_holes(AERIAL, AERIAL_PARAMS, aerial_hs, OUT / "fig_aerial_holes.png")
    print(f"      {len(aerial_hs.holes)} holes")

    print("[2/5] tube ceiling hole detection ...")
    tube_hs = detect_holes(TUBE, **TUBE_PARAMS)
    tube_hs.to_json(OUT / "tube_holes.json")
    report["tube_holes"] = hole_summary(tube_hs)
    render_holes(TUBE, TUBE_PARAMS, tube_hs, OUT / "fig_tube_holes.png")
    print(f"      {len(tube_hs.holes)} holes")

    print("[3/5] constellation match ...")
    tf = match_constellations(aerial_hs, tube_hs, mode="auto", tolerance=TOLERANCE)
    tf.to_json(OUT / "transform_landmark.json")
    vres = vertical_residuals(aerial_hs, tube_hs, tf)
    report["landmark"] = dict(
        mode=tf.mode, matches=len(tf.inliers), rms=round(tf.rms, 3),
        margin=tf.margin, shape_score=round(tf.shape_score, 3),
        inliers=[list(p) for p in tf.inliers],
        vertical_residuals_m=[round(v, 3) for v in vres],
        warnings=list(tf.warnings),
    )
    print(f"      mode={tf.mode} matches={len(tf.inliers)} RMS={tf.rms:.3f} "
          f"margin={tf.margin} vres={[f'{v:+.2f}' for v in vres]}")
    render_match(aerial_hs, tube_hs, tf, OUT / "fig_match.png")

    print("[4/5] rim GICP refinement (4-DOF and 6-DOF) ...")
    ref4 = refine_variant(tf, 4, "d4", report, aerial_hs, tube_hs)
    refine_variant(tf, 6, "d6", report, aerial_hs, tube_hs)

    print("[5/5] scale/drift ruler (inter-skylight distances) ...")
    A = aerial_hs.centroids()
    B = tube_hs.centroids()
    pairs, ruler = list(tf.inliers), []
    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            (t1, a1), (t2, a2) = pairs[i], pairs[j]
            da = float(np.linalg.norm(A[a1] - A[a2]))
            dt = float(np.linalg.norm(B[t1] - B[t2]))
            ruler.append(dict(pair=f"a{a1}-a{a2}", aerial_m=round(da, 3),
                              tube_m=round(dt, 3), diff_m=round(dt - da, 3),
                              diff_pct=round(100 * (dt - da) / da, 3)))
    report["scale_ruler"] = ruler

    (OUT / "registration_report.json").write_text(json.dumps(report, indent=2))
    print(f"wrote {OUT/'registration_report.json'}")


if __name__ == "__main__":
    sys.exit(main())
