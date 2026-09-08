"""Chain the SECOND long flight into the registered frame and compare roofs.

second_long_flight could never be registered through the skylights alone
(near-collinear roll ambiguity, see slf_registration_finding.json). Here it
is carried into the slice-registered frame the same way as the first flight:
a map-to-map GICP against the DLIO first-flight cloud (the two flights
overlap over the entrance and skylight reach), composed with Ed's transform.
The per-column ceiling difference against the first flight in the common
frame is the closest thing the campaign has to a flight-to-flight
repeatability measurement.

Outputs: analysis_out/transform_slf_ed_chain.json, maps/slf_10cm_ed.pcd,
         analysis_out_v6/flight2_repeatability.json

Run:  PYTHONPATH=src .venv/bin/python analysis/chain_slf_to_ed.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import small_gicp
from scipy.spatial import cKDTree

from lava_pcd.io.pcd_reader import BinaryPcdReader
from lava_pcd.merge import apply_transform

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
# repeatability output dir env-selectable (default unchanged).
OUT6 = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v6")))
# v9 datum-transfer rerun: SLF_CHAIN_REUSE=1 skips the GICP fit and the map
# rewrite, reuses the committed transform_slf_ed_chain.json, and recomputes the
# repeatability from the (env-selectable) chained maps along OUT6's centreline.
REUSE = os.environ.get("SLF_CHAIN_REUSE", "") == "1"
F1_PCD = Path(os.environ.get("SLF_F1_PCD", str(ROOT / "maps/flf_30cm_ed.pcd")))
F2_PCD = os.environ.get("SLF_F2_PCD")   # only used with SLF_CHAIN_REUSE=1


def load(p: Path) -> np.ndarray:
    parts = []
    with BinaryPcdReader(p) as r:
        for c in r.chunks():
            parts.append(c[:, :3].astype(np.float64))
    return np.vstack(parts)


def main() -> None:
    ED = np.array(json.loads((OUT / "transform_ed_slice.json").read_text())["matrix"])
    SLF4 = np.array(json.loads((OUT / "slf_transform_landmark.json").read_text())["matrix"])

    if REUSE:
        CHAIN = np.array(json.loads(
            (OUT / "transform_slf_ed_chain.json").read_text())["matrix"])
        print("SLF_CHAIN_REUSE=1: reusing transform_slf_ed_chain.json, no GICP, "
              "no map rewrite")
        repeatability(CHAIN, None)
        return

    dlio = load(ROOT / "maps/tube_30cm.pcd")     # target (first flight, DLIO frame)
    slf = load(ROOT / "maps/slf_30cm.pcd")       # source (second flight, FAST-LIO frame)
    print(f"target {len(dlio):,} pts, source {len(slf):,} pts")

    T = np.linalg.inv(ED) @ SLF4                 # init through the two registrations
    stats = []
    for max_dist, it in ((8.0, 40), (2.0, 40)):
        res = small_gicp.align(dlio, slf, init_T_target_source=T,
                               registration_type="GICP",
                               downsampling_resolution=0.3,
                               max_correspondence_distance=max_dist,
                               max_iterations=it, num_threads=8)
        T = np.asarray(res.T_target_source)
        stats.append(dict(max_dist=max_dist, converged=bool(res.converged),
                          iterations=int(res.iterations),
                          num_inliers=int(res.num_inliers)))
        print(f"GICP max_dist={max_dist}: converged={res.converged} "
              f"inliers={res.num_inliers}")

    CHAIN = ED @ T
    tilt = float(np.degrees(np.arccos(np.clip(CHAIN[2, 2], -1, 1))))
    payload = dict(matrix=[list(map(float, r)) for r in CHAIN],
                   inliers=[], rms=0.0, mode="ed_chain_gicp_slf",
                   aerial_origin=[479158.0, 7089826.0, 0.0],
                   tube_origin=[0.0, 0.0, 0.0], n_inliers=0, margin=0,
                   shape_score=0.0, aerial_up=[0.0, 0.0, 1.0], anchors=[],
                   anchor_axes=[], gicp_stages=stats,
                   warnings=[f"tilt {tilt:.2f} deg",
                             "second_long_flight -> registered frame via "
                             "map-to-map GICP against the first-flight cloud"])
    (OUT / "transform_slf_ed_chain.json").write_text(json.dumps(payload, indent=2))

    for res_tag in ("30cm", "10cm"):
        apply_transform(ROOT / f"maps/slf_{res_tag}.pcd",
                        ROOT / f"maps/slf_{res_tag}_ed.pcd", CHAIN,
                        show_progress=False)
    repeatability(CHAIN, slf)


def repeatability(CHAIN: np.ndarray, slf: np.ndarray | None) -> None:
    # flight-to-flight ceiling repeatability along the v6 centreline
    import csv
    rows = list(csv.DictReader(open(OUT6 / "centreline.csv")))
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    s = np.array([float(r["s"]) for r in rows])
    f1 = load(F1_PCD)                            # flight 1, chained
    if slf is None:                              # reuse mode: chained map on disk
        f2 = load(Path(F2_PCD) if F2_PCD else ROOT / "maps/slf_30cm_ed.pcd")
    else:
        f2 = slf @ CHAIN[:3, :3].T + CHAIN[:3, 3]    # flight 2, chained
    t1, t2 = cKDTree(f1[:, :2]), cKDTree(f2[:, :2])
    off = np.full(len(x), np.nan)
    for i in range(len(x)):
        a = t1.query_ball_point([x[i], y[i]], 2.5)
        b = t2.query_ball_point([x[i], y[i]], 2.5)
        if len(a) > 30 and len(b) > 30:
            off[i] = np.percentile(f1[a, 2], 99) - np.percentile(f2[b, 2], 99)
    ok = ~np.isnan(off)
    rep = dict(n_stations=int(ok.sum()),
               s_range=[float(s[ok].min()), float(s[ok].max())] if ok.any() else None,
               median=float(np.nanmedian(off)),
               rms=float(np.sqrt(np.nanmean(off[ok] ** 2))),
               p5=float(np.nanpercentile(off, 5)),
               p95=float(np.nanpercentile(off, 95)),
               note="flight1 minus flight2 p99 ceiling in 2.5 m columns "
                    "along the v6 centreline, both flights in the registered "
                    "frame; coverage limited to the reach flight 2 mapped")
    (OUT6 / "flight2_repeatability.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
