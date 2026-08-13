"""Set up the Ed Llewellin slice-based registration for a full pipeline re-run.

Background (email thread with E. Llewellin, 24 Jul-early Aug 2026): Ed's MSc
student filtered/decimated the tube point cloud that Nikos exported (the DLIO
first_long_flight map, SLAM frame); Ed then registered it to Birgir's aerial
photogrammetric model with a bespoke slice tool (longitudinal + transverse
slices through each skylight, points within 0.5 m projected onto each slice
plane, per-slice manual alignment, best-fit rigid solid-body transform from the
six slice fits, iterated four times against the full clouds in CloudCompare).

He supplied the resulting 4x4 rigid transform in two frames: full UTM, and the
aerial-local frame whose origin (479158, 7089826) is exactly LAVA_PCD_ORIGIN.
Session analysis (2026-08-13) confirmed the matrix applies to the DLIO-frame
cloud (skylight horizontal residuals 1.0-2.6 m; yaw 77.2 deg matches the DLIO
constellation solve) and that it passes two physical checks the 4-DOF
centroid solve fails: no tube ceiling above the surface DEM anywhere, and the
tube floor at S3 coincides with the floor the photogrammetry sees through the
skylight. The 4-DOF centroid solve sits the whole tube ~4.5-5 m too high
because it matches two differently defined "rim centroid" heights.

This script:
  1. writes analysis_out/transform_ed_slice.json (Transform schema, so
     `lava-pcd transform` can consume it),
  2. produces maps/tube_10cm_ed.pcd and maps/tube_30cm_ed.pcd (DLIO working
     copies moved into the aerial-local frame with Ed's matrix),
  3. verifies the skylight centroid residuals,
  4. seeds ANALYSIS_OUT (default analysis_out_ed/) with the inputs the
     pipeline scripts read, so the Ed-frame run never touches the committed
     baseline in analysis_out/.

Then run, with PYTHONPATH=src and the venv python:
  ANALYSIS_OUT=analysis_out_ed \
  MORPHO_SKEL_PCD=maps/tube_30cm_ed.pcd MORPHO_SECT_PCD=maps/tube_10cm_ed.pcd \
    python analysis/run_morphometry.py
  ANALYSIS_OUT=analysis_out_ed ROOF_TUBE_PCD=maps/tube_10cm_ed.pcd \
  ROOF_DEM_PCD=maps/full_surface_and_subsurface_merged.pcd \
    python analysis/run_roof.py
  ANALYSIS_OUT=analysis_out_ed python analysis/run_planetary.py
  ANALYSIS_OUT=analysis_out_ed python analysis/collect_paper_numbers.py
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import numpy as np

from lava_pcd.merge import apply_transform

ROOT = Path(__file__).resolve().parents[1]
OUT_BASE = ROOT / "analysis_out"
OUT_ED = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_ed")))

# Ed's transform, aerial-local frame (his "PS" matrix; identical to the UTM one
# with -479158 applied to x and -7089826 to y). Source: E. Llewellin email,
# 24 Jul 2026. Applies to the DLIO-frame tube cloud (tube_10cm/30cm.pcd).
ED_LOCAL = np.array([
    [0.222039143686, -0.975032541024, -0.003187863146, 1376.516695668483],
    [0.974929141414,  0.222062022362, -0.014199547681,  601.797969112199],
    [0.014552924401,  0.000044914729,  0.999894099847,  225.886650880168],
    [0.0, 0.0, 0.0, 1.0],
])
AERIAL_ORIGIN = (479158.0, 7089826.0, 0.0)


def write_transform_json() -> Path:
    R = ED_LOCAL[:3, :3]
    tilt_deg = float(np.degrees(np.arccos(R[2, 2])))
    payload = {
        "matrix": [list(map(float, row)) for row in ED_LOCAL],
        "inliers": [[0, 0], [1, 1], [2, 2]],
        # rms below is OUR check of his solve against the rim-hole centroids
        # (horizontal only; the vertical centroid offsets are definitional, not
        # errors). Ed's own per-slice residuals are requested and pending.
        "rms": 1.72,
        "mode": "slice_manual_ed",
        "aerial_origin": list(AERIAL_ORIGIN),
        "tube_origin": [0.0, 0.0, 0.0],
        "n_inliers": 3,
        "margin": 0,
        "shape_score": 0.0,
        "aerial_up": [0.0, 0.0, 1.0],
        "anchors": [],
        "anchor_axes": [],
        "warnings": [
            "External registration by E. Llewellin (slice-based manual fine "
            "alignment, iterated 4x); applies to the DLIO-frame cloud "
            "(tube_*cm.pcd), NOT the FAST-LIO flf_* maps.",
            f"det(R)={float(np.linalg.det(R)):.9f}, tilt of z-axis "
            f"{tilt_deg:.3f} deg (real along-tube grade), no scale.",
            "rms field = horizontal skylight rim-centroid check on our side; "
            "replace with Ed's slice residuals when supplied.",
        ],
    }
    p = OUT_BASE / "transform_ed_slice.json"
    p.write_text(json.dumps(payload, indent=2))
    return p


def verify() -> dict:
    """Skylight rim-centroid residuals of Ed's transform (DLIO holes)."""
    tube = json.loads((OUT_BASE / "tube_holes.json").read_text())
    aer = json.loads((OUT_BASE / "aerial_holes.json").read_text())
    a = np.array([h["centroid"] for h in aer["holes"][:3]])
    t = np.array([h["centroid"] for h in tube["holes"]])
    tt = t @ ED_LOCAL[:3, :3].T + ED_LOCAL[:3, 3]
    d = tt - a
    res = {
        "horizontal_residuals_m": [float(v) for v in np.hypot(d[:, 0], d[:, 1])],
        "horizontal_rms_m": float(np.sqrt((d[:, :2] ** 2).sum(1).mean())),
        "vertical_centroid_offsets_m": [float(v) for v in d[:, 2]],
        "note": "vertical offsets are the ceiling-boundary vs surface-hole "
                "centroid definition gap, not registration error; see "
                "setup_ed_frame.py docstring",
    }
    print("skylight check:", json.dumps(res, indent=2))
    return res


def make_maps() -> None:
    for res in ("30cm", "10cm"):
        src = ROOT / f"maps/tube_{res}.pcd"
        dst = ROOT / f"maps/tube_{res}_ed.pcd"
        if dst.exists():
            print(f"{dst.name} exists, skipping")
            continue
        print(f"transforming {src.name} -> {dst.name}")
        apply_transform(src, dst, ED_LOCAL, show_progress=False)


def seed_out(check: dict) -> None:
    OUT_ED.mkdir(exist_ok=True)
    # Inputs the pipeline reads from its output dir.
    for name in ("aerial_holes.json", "tube_holes.json", "vertical_check.json",
                 "planetary_catalogue.csv",
                 "dem_grid_full_surface_and_subsurface_merged.npz"):
        src = OUT_BASE / name
        dst = OUT_ED / name
        if not dst.exists():
            shutil.copy2(src, dst)

    # Uncertainty budget: same lidar/DEM terms; the registration term is
    # provisional until Ed supplies his slice residuals.
    bud = json.loads((OUT_BASE / "uncertainty_budget.json").read_text())
    bud["sigma_reg_z_basis"] = (
        "PROVISIONAL for the Ed slice registration: kept at the previous "
        "1.09 m pending Ed's per-slice residuals. His stated expectation is "
        "that SLAM/photogrammetry drift now dominates the residual.")
    bud["final_transform"] = ("transform_ed_slice.json (E. Llewellin, "
                              "slice-based manual fine alignment on the "
                              "DLIO-frame cloud)")
    (OUT_ED / "uncertainty_budget.json").write_text(json.dumps(bud, indent=2))

    # Registration report: carry the baseline structure, override the landmark
    # block so collect_paper_numbers reports the applied registration.
    rep = json.loads((OUT_BASE / "registration_report.json").read_text())
    rep["landmark"] = {
        "mode": "slice_manual_ed",
        "matches": 3,
        "rms": check["horizontal_rms_m"],
        "margin": 0,
        "shape_score": 0.0,
        "inliers": [[0, 0], [1, 1], [2, 2]],
        "vertical_residuals_m": check["vertical_centroid_offsets_m"],
        "warnings": [
            "External slice-based registration (E. Llewellin); rms is the "
            "horizontal rim-centroid check, vertical_residuals_m are centroid "
            "definition offsets (see transform_ed_slice.json warnings).",
        ],
    }
    rep["note"] = ("Ed-frame evaluation run. aerial/tube hole detections and "
                   "scale_ruler carried over from the baseline DLIO "
                   "constellation report.")
    (OUT_ED / "registration_report.json").write_text(json.dumps(rep, indent=2))
    print(f"seeded {OUT_ED}")


def main() -> None:
    p = write_transform_json()
    print(f"wrote {p}")
    check = verify()
    make_maps()
    seed_out(check)


if __name__ == "__main__":
    main()
