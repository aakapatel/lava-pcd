"""Validate the Llewellin slice-based registration and set its uncertainty.

The 4-DOF constellation solve was validated against the skylight rim
centroids it was fitted to, which is circular and, as it turned out, biased:
the tube ceiling-boundary centroid and the aerial surface-hole centroid are
different physical heights, and matching them lifted the whole map ~4.5 m.
This script validates Ed's slice registration against measurements it was
NOT fitted to on our side, and derives the vertical registration uncertainty
entering the roof thickness budget from those residuals.

Checks:
  1. Floor-through-skylight residual (S3). The aerial photogrammetry images
     the tube floor through the largest skylight; the tube lidar mapped the
     same floor. For every aerial below-rim point near S3, the signed vertical
     offset to the local tube-floor surface is computed. This is a direct,
     fit-independent probe of the vertical datum at the anchor.
  2. Rim-centroid horizontal residuals (from setup_ed_frame.py's check).
  3. Two-SLAM consistency: per-zone ceiling agreement between the Ed-frame
     DLIO cloud and the FAST-LIO cloud chained into the same frame
     (chain_flf_to_ed.py). Bounds SLAM-drift error away from the anchors.
  4. Physical non-negativity: stations with tau < 0 over intact ground.

sigma_reg_z is set to the quadrature of the S3 floor residual RMS (datum at
the anchors) and the worst-zone two-SLAM std (drift away from them); it is
provisional until Ed supplies his per-slice residuals. The budget and the
roof outputs in ANALYSIS_OUT (default analysis_out_ed) are then refreshed by
re-running run_roof.py and collect_paper_numbers.py.

Run:  PYTHONPATH=src .venv/bin/python analysis/validate_ed_registration.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT_BASE = ROOT / "analysis_out"
OUT_ED = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_ed")))


def load(p: Path, keep=None) -> np.ndarray:
    parts = []
    with BinaryPcdReader(p) as r:
        for c in r.chunks():
            q = c[:, :3].astype(np.float64)
            if keep is not None:
                q = q[keep(q)]
            parts.append(q)
    return np.vstack(parts)


def main() -> None:
    aer_holes = json.loads((OUT_BASE / "aerial_holes.json").read_text())
    s3 = np.array(aer_holes["holes"][2]["centroid"])  # largest skylight
    r_disc = 5.0

    def near_s3(q):
        return np.hypot(q[:, 0] - s3[0], q[:, 1] - s3[1]) < r_disc

    aer = load(ROOT / "maps/aerial_10cm.pcd", keep=near_s3)
    tube = load(ROOT / "maps/tube_10cm_ed.pcd", keep=near_s3)

    # Aerial below-rim (through-skylight) points: below the local surface.
    surf = np.percentile(aer[:, 2], 90)
    a_in = aer[aer[:, 2] < surf - 3.0]

    # Tube floor surface: per-column lowest band of the interior cloud.
    t_floor = []
    txy = cKDTree(tube[:, :2])
    dz = []
    for p in a_in:
        idx = txy.query_ball_point(p[:2], 0.5)
        if len(idx) < 5:
            continue
        z = tube[idx, 2]
        floor = np.median(z[z < np.percentile(z, 30)])  # local floor band
        dz.append(p[2] - floor)
    dz = np.array(dz)
    floor_stats = dict(
        n=int(len(dz)),
        median_m=float(np.median(dz)),
        rms_m=float(np.sqrt((dz ** 2).mean())),
        p10_m=float(np.percentile(dz, 10)),
        p90_m=float(np.percentile(dz, 90)),
    )
    print("S3 floor-through-skylight residual (aerial minus tube floor):")
    print(json.dumps(floor_stats, indent=2))

    chain = json.loads((OUT_BASE / "transform_flf_ed_chain.json").read_text())
    zones = chain["ceiling_offset_validation"]
    worst_std = max(v["std"] for v in zones.values())
    worst_mean = max(abs(v["mean"]) for v in zones.values())
    print(f"two-SLAM ceiling agreement: worst zone |mean| {worst_mean:.2f} m, "
          f"worst std {worst_std:.2f} m")

    rep = json.loads((OUT_ED / "registration_report.json").read_text())
    rim_h = rep["landmark"]["rms"]

    rows = [l.split(",") for l in
            (OUT_ED / "roof_thickness.csv").read_text().splitlines()[1:]]
    # columns: s,x,y,z_dem,z_ceil,tau,sigma_tau,span,coverage,ghost_frac,class
    neg = [r for r in rows if r[-1] == "inconsistent"]

    sigma_anchor = floor_stats["rms_m"]
    sigma_drift = float(np.sqrt(worst_mean ** 2 + worst_std ** 2))
    sigma_reg = float(np.sqrt(sigma_anchor ** 2 + sigma_drift ** 2))

    validation = dict(
        floor_through_skylight_s3=floor_stats,
        rim_horizontal_rms_m=rim_h,
        two_slam_zones=zones,
        n_negative_tau_stations=len(neg),
        sigma_reg_components=dict(anchor_floor_rms_m=sigma_anchor,
                                  drift_rss_m=sigma_drift),
        sigma_reg_z_m=round(sigma_reg, 2),
        note="sigma_reg = quadrature(anchor floor RMS at S3, worst-zone "
             "two-SLAM RSS). PROVISIONAL until Ed's per-slice residuals "
             "arrive; validated against measurements the registration was "
             "not fitted to on our side.",
    )
    (OUT_ED / "registration_validation.json").write_text(
        json.dumps(validation, indent=2))
    print(f"sigma_reg_z = {sigma_reg:.2f} m "
          f"(anchor {sigma_anchor:.2f} + drift {sigma_drift:.2f} in quadrature)")

    bud = json.loads((OUT_ED / "uncertainty_budget.json").read_text())
    bud["sigma_reg_z_m"] = round(sigma_reg, 2)
    bud["sigma_reg_z_basis"] = validation["note"]
    bud["sigma_tau_m"] = round(float(np.sqrt(
        sigma_reg ** 2 + bud["sigma_lidar_z_m"] ** 2 + bud["sigma_dem_z_m"] ** 2)), 2)
    bud["horizontal_rim_rms_m"] = rim_h
    (OUT_ED / "uncertainty_budget.json").write_text(json.dumps(bud, indent=2))
    print(f"updated budget: sigma_tau = {bud['sigma_tau_m']} m")


if __name__ == "__main__":
    main()
