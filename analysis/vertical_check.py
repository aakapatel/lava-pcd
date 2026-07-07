"""S1.2d: per-skylight vertical agreement between the two clouds.

For each matched skylight, gather the rim band in both clouds (points inside the
skylight's inflated ellipse, within +-band of the opening level, tube points
mapped through the final transform) and split it into the *lip* percentile bands.
The overlap surface both sensors observe is the upper rim lip, so we compare
upper-quantile elevations rather than means (the aerial sees the surface plus
the snow-cone floor through the hole; the tube lidar sees the underside).

Output: analysis_out/vertical_check.json with, per skylight, the aerial and
tube elevation quantiles inside the collar and their differences. The spread of
the q90 differences across the three skylights is the vertical registration
uncertainty entering the roof-thickness quadrature (sigma_reg).

Run: .venv/bin/python analysis/vertical_check.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from lava_pcd.geometry import transform_points
from lava_pcd.io.pcd_reader import BinaryPcdReader
from lava_pcd.register import Transform

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
BAND = 3.0        # vertical band around the opening level (m)
INFLATE = 1.6     # ellipse inflation to catch the lip
QS = (50, 75, 90, 95)


def collar_mask(xy: np.ndarray, cx: float, cy: float, a: float, b: float,
                ex: float, ey: float) -> np.ndarray:
    """Points inside the inflated ellipse (axes a,b along (ex,ey) / perp)."""
    dx, dy = xy[:, 0] - cx, xy[:, 1] - cy
    u = dx * ex + dy * ey
    v = -dx * ey + dy * ex
    return (u / (INFLATE * a)) ** 2 + (v / (INFLATE * b)) ** 2 <= 1.0


def gather(path: Path, anchors: np.ndarray, axes: np.ndarray,
           matrix: np.ndarray | None) -> list[np.ndarray]:
    per_hole: list[list[np.ndarray]] = [[] for _ in anchors]
    with BinaryPcdReader(path) as reader:
        for chunk in reader.chunks():
            pts = chunk[:, :3].astype(np.float64)
            if matrix is not None:
                pts = transform_points(pts, matrix)
            for j, (anc, ax) in enumerate(zip(anchors, axes)):
                a, b, ex, ey = ax[0], max(ax[1], 0.5), ax[2], ax[3]
                near = np.abs(pts[:, 2] - anc[2]) <= BAND
                if not near.any():
                    continue
                sel = pts[near]
                m = collar_mask(sel[:, :2], anc[0], anc[1], a, b, ex, ey)
                if m.any():
                    per_hole[j].append(sel[m][:, 2])
    return [np.concatenate(z) if z else np.empty(0) for z in per_hole]


def main() -> None:
    tf = Transform.from_json(OUT / "transform_landmark.json")
    anchors = np.asarray(tf.anchors, float)
    axes = np.asarray(tf.anchor_axes, float)
    za = gather(ROOT / "maps/aerial_crop.pcd", anchors, axes, None)
    zt = gather(ROOT / "maps/tube_10cm.pcd", anchors, axes, tf.array)

    report = []
    for j, (a_z, t_z) in enumerate(zip(za, zt)):
        entry = dict(skylight=j, n_aerial=int(len(a_z)), n_tube=int(len(t_z)))
        for q in QS:
            qa = float(np.percentile(a_z, q)) if len(a_z) else None
            qt = float(np.percentile(t_z, q)) if len(t_z) else None
            entry[f"q{q}_aerial"] = round(qa, 3) if qa is not None else None
            entry[f"q{q}_tube"] = round(qt, 3) if qt is not None else None
            entry[f"q{q}_diff"] = (round(qt - qa, 3)
                                   if qa is not None and qt is not None else None)
        report.append(entry)
        print(entry)

    diffs = [e["q90_diff"] for e in report if e["q90_diff"] is not None]
    summary = dict(
        per_skylight=report,
        q90_diff_mean=round(float(np.mean(diffs)), 3),
        q90_diff_std=round(float(np.std(diffs)), 3),
        q90_diff_max_abs=round(float(np.max(np.abs(diffs))), 3),
        band_m=BAND, inflate=INFLATE,
        note="tube-minus-aerial rim-lip elevation; overlap surface is the lip, "
             "so upper quantiles are the comparable statistic",
    )
    (OUT / "vertical_check.json").write_text(json.dumps(summary, indent=2))
    print("q90 diffs:", diffs)


if __name__ == "__main__":
    main()
