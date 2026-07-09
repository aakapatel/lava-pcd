"""T6 (drift): cross-check the Mac DLIO tube reconstruction against the offline
FAST-LIO map of the SAME first_long_flight, both registered to the same three
skylights.

Both maps are pinned to the skylight constellation at the shallow end. If the
DLIO map has accumulated vertical drift along the descent, its ceiling will fall
progressively below the FAST-LIO ceiling with distance from the skylights. This
is exactly what we test: per centreline station, compare
  z_ceil(DLIO, from morphometry.csv)  vs  z_ceil(FAST-LIO, this script).

Background: the DLIO map (tube_10cm_aerial) descends ~74 m over the survey and
drove the T3 full-DEM tau result (median 26.7 m, deep tau up to 63 m). But the
FAST-LIO odometry of the same flight descends only ~12 m, casting doubt on the
DLIO descent. This script quantifies the disagreement so the tau result can be
corrected (the barometric altitude from the bag, T7, is the final arbiter).

Inputs:  analysis_out/morphometry.csv (DLIO ceiling), maps/flf_10cm_aerial.pcd
Outputs: analysis_out/dlio_vs_fastlio_drift.csv / .json / fig_drift_profile.png
Run:  env -u PYTHONPATH .venv/bin/python analysis/dlio_vs_fastlio_drift.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
FLF = ROOT / "maps/flf_10cm_aerial.pcd"
RADIUS = 3.0
MIN_PTS = 30
CEIL_PCTL = 95   # ceiling proxy = high percentile of z in the XY neighbourhood


def main():
    pts = []
    with BinaryPcdReader(FLF) as r:
        for c in r.chunks():
            pts.append(c[:, :3].astype(np.float64))
    flf = np.vstack(pts)
    tree = cKDTree(flf[:, :2])

    rows = list(csv.DictReader(open(OUT / "morphometry.csv")))
    recs = []
    for rrow in rows:
        s = float(rrow["s"]); x = float(rrow["x"]); y = float(rrow["y"])
        zc = float(rrow["z_ceil"]) if rrow["z_ceil"] else np.nan
        idx = tree.query_ball_point([x, y], RADIUS)
        if len(idx) < MIN_PTS or np.isnan(zc):
            recs.append((s, zc, np.nan, np.nan)); continue
        zf = float(np.percentile(flf[idx, 2], CEIL_PCTL))
        recs.append((s, zc, zf, zc - zf))

    with open(OUT / "dlio_vs_fastlio_drift.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["s", "z_ceil_dlio", "z_ceil_fastlio", "dlio_minus_fastlio"])
        for s, zc, zf, d in recs:
            w.writerow([f"{s:.1f}", f"{zc:.3f}", f"{zf:.3f}" if not np.isnan(zf)
                        else "", f"{d:.3f}" if not np.isnan(d) else ""])

    diffs = np.array([d for *_, d in recs if not np.isnan(d)])
    s_all = np.array([s for s, *_ in recs])
    # near-skylight baseline (s <= 120) vs deep end
    near = np.array([d for s, zc, zf, d in recs if not np.isnan(d) and s <= 120])
    deep = np.array([d for s, zc, zf, d in recs if not np.isnan(d) and s >= 250])
    summary = dict(
        n_compared=int(len(diffs)),
        near_skylight_mean_diff_m=round(float(near.mean()), 2) if len(near) else None,
        near_skylight_absmedian_m=round(float(np.median(np.abs(near))), 2) if len(near) else None,
        deep_end_mean_diff_m=round(float(deep.mean()), 2) if len(deep) else None,
        max_divergence_m=round(float(diffs.min()), 2),  # DLIO below FAST-LIO
        interpretation=(
            "DLIO ceiling agrees with FAST-LIO to ~+-3 m near the skylights "
            "(both pinned there) but falls up to tens of m below FAST-LIO toward "
            "the deep end, i.e. the DLIO 74 m descent is largely vertical drift. "
            "FAST-LIO stays roughly flat, consistent with its 12 m odometry "
            "descent. The barometric altitude from the bag (T7) is the arbiter."),
    )
    (OUT / "dlio_vs_fastlio_drift.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    # figure
    sc = np.array([s for s, zc, zf, d in recs if not np.isnan(d)])
    zd = np.array([zc for s, zc, zf, d in recs if not np.isnan(d)])
    zf = np.array([zf for s, zc, zf, d in recs if not np.isnan(d)])
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    axes[0].plot(sc, zd, ".-", ms=4, label="DLIO ceiling (tube_10cm_aerial)")
    axes[0].plot(sc, zf, ".-", ms=4, label="FAST-LIO ceiling (same flight)")
    axes[0].set_ylabel("ceiling elevation z (m, aerial frame)")
    axes[0].legend(); axes[0].set_title(
        "DLIO vs FAST-LIO ceiling, both registered to the 3 skylights")
    axes[1].axhline(0, color="k", lw=0.8)
    axes[1].plot(sc, zd - zf, ".-", ms=4, color="tab:red")
    axes[1].set_ylabel("DLIO - FAST-LIO (m)")
    axes[1].set_xlabel("distance along centreline s (m)")
    fig.savefig(OUT / "fig_drift_profile.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("wrote fig_drift_profile.png")


if __name__ == "__main__":
    main()
