"""Chain the FAST-LIO first-flight map into Ed's registration frame.

Ed's slice-based transform (analysis_out/transform_ed_slice.json, see
setup_ed_frame.py) applies to the DLIO-frame cloud only. The FAST-LIO map of
the same flight is our cleaner reconstruction (no multipass ghost near the
skylights), so to use Ed's registration with it we need the rigid transform
between the two SLAM frames:

    T_aerial<-flf  =  ED_LOCAL  @  T_dlio<-flf

T_dlio<-flf is estimated by GICP between the 30 cm working copies of the two
maps (same rock, two SLAM solutions), initialised from the two skylight
registrations (ED_LOCAL^-1 @ FLF_4DOF) and refined in two stages. The GICP
residual also measures the relative drift between the two SLAM solutions,
which bounds how far "one rigid transform" can be trusted.

Outputs:
    analysis_out/transform_flf_ed_chain.json   (Transform schema + GICP stats)
    maps/flf_10cm_ed.pcd, maps/flf_30cm_ed.pcd
    printed validation: per-zone ceiling offset between the Ed-frame DLIO and
    Ed-frame FAST-LIO clouds along the baseline centreline (expect ~0 if the
    chain is consistent; structure in s = relative SLAM drift).

Run:  PYTHONPATH=src .venv/bin/python analysis/chain_flf_to_ed.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import small_gicp
from scipy.spatial import cKDTree

from lava_pcd.io.pcd_reader import BinaryPcdReader
from lava_pcd.merge import apply_transform

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"


def load(p: Path) -> np.ndarray:
    parts = []
    with BinaryPcdReader(p) as r:
        for c in r.chunks():
            parts.append(c[:, :3].astype(np.float64))
    return np.vstack(parts)


def main() -> None:
    ED = np.array(json.loads((OUT / "transform_ed_slice.json").read_text())["matrix"])
    FLF = np.array(json.loads((OUT / "flf_transform_landmark.json").read_text())["matrix"])

    dlio = load(ROOT / "maps/tube_30cm.pcd")     # target (DLIO SLAM frame)
    flf = load(ROOT / "maps/flf_30cm.pcd")       # source (FAST-LIO SLAM frame)
    print(f"target (DLIO) {len(dlio):,} pts, source (FAST-LIO) {len(flf):,} pts")

    # Init: through the two aerial registrations. The ~4.5 m vertical bias of
    # the 4-DOF solve is inside GICP's coarse capture range.
    T0 = np.linalg.inv(ED) @ FLF
    stats = []
    T = T0
    for max_dist, it in ((6.0, 30), (1.5, 30)):
        res = small_gicp.align(dlio, flf, init_T_target_source=T,
                               registration_type="GICP",
                               downsampling_resolution=0.3,
                               max_correspondence_distance=max_dist,
                               max_iterations=it, num_threads=8)
        T = np.asarray(res.T_target_source)
        stats.append(dict(max_dist=max_dist, converged=bool(res.converged),
                          iterations=int(res.iterations),
                          error=float(res.error),
                          num_inliers=int(res.num_inliers)))
        print(f"GICP stage max_dist={max_dist}: converged={res.converged} "
              f"iters={res.iterations} inliers={res.num_inliers}")

    d = T @ np.linalg.inv(T0)
    print(f"correction w.r.t. init: dt={np.round(d[:3,3],3)} m, "
          f"rot {np.degrees(np.arccos(np.clip((np.trace(d[:3,:3])-1)/2,-1,1))):.3f} deg")

    CHAIN = ED @ T   # aerial <- flf
    Rz = CHAIN[:3, :3]
    payload = {
        "matrix": [list(map(float, r)) for r in CHAIN],
        "inliers": [], "rms": 0.0, "mode": "ed_chain_gicp",
        "aerial_origin": [479158.0, 7089826.0, 0.0],
        "tube_origin": [0.0, 0.0, 0.0],
        "n_inliers": 0, "margin": 0, "shape_score": 0.0,
        "aerial_up": [0.0, 0.0, 1.0], "anchors": [], "anchor_axes": [],
        "warnings": [
            "T_aerial<-flf = ED_LOCAL @ T_dlio<-flf(GICP). Carries the "
            "FAST-LIO map into Ed's slice-registered frame.",
            f"tilt of z-axis {float(np.degrees(np.arccos(Rz[2,2]))):.3f} deg",
        ],
        "gicp_stages": stats,
        "T_dlio_from_flf": [list(map(float, r)) for r in T],
    }
    p = OUT / "transform_flf_ed_chain.json"
    p.write_text(json.dumps(payload, indent=2))
    print(f"wrote {p}")

    for res_tag in ("30cm", "10cm"):
        src = ROOT / f"maps/flf_{res_tag}.pcd"
        dst = ROOT / f"maps/flf_{res_tag}_ed.pcd"
        print(f"transforming {src.name} -> {dst.name}")
        apply_transform(src, dst, CHAIN, show_progress=False)

    # Validation: per-station p99 ceiling of the two Ed-frame clouds along the
    # baseline centreline. Offset ~0 => chain consistent; trend in s =>
    # relative SLAM drift between the DLIO and FAST-LIO solutions.
    import csv
    rows = list(csv.DictReader(open(OUT / "centreline.csv")))
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    s = np.array([float(r["s"]) for r in rows])
    de = dlio @ ED[:3, :3].T + ED[:3, 3]
    fe = flf @ CHAIN[:3, :3].T + CHAIN[:3, 3]
    td, tf = cKDTree(de[:, :2]), cKDTree(fe[:, :2])
    off = np.full(len(x), np.nan)
    for i in range(len(x)):
        a = td.query_ball_point([x[i], y[i]], 2.5)
        b = tf.query_ball_point([x[i], y[i]], 2.5)
        if len(a) > 30 and len(b) > 30:
            off[i] = np.percentile(de[a, 2], 99) - np.percentile(fe[b, 2], 99)
    print("ceiling offset (DLIO_ed minus FLF_ed) by zone:")
    zones = {}
    for lo, hi in ((0, 80), (80, 160), (160, 240), (240, 340)):
        m = (s >= lo) & (s < hi) & ~np.isnan(off)
        if m.sum():
            zones[f"s{lo}-{hi}"] = dict(mean=float(np.nanmean(off[m])),
                                        std=float(np.nanstd(off[m])),
                                        n=int(m.sum()))
            print(f"  s {lo:3d}-{hi:3d}: mean {np.nanmean(off[m]):+6.2f} m "
                  f"std {np.nanstd(off[m]):.2f} (n={m.sum()})")
    payload["ceiling_offset_validation"] = zones
    p.write_text(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
