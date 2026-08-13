"""Full-pipeline repeatability: the same flight processed end to end through
two independent lidar-inertial estimators (onboard DLIO in analysis_out_ed,
offline FAST-LIO2 in analysis_out_v6), each with its own L1-medial centreline,
cross-sections and roof profile, in the common registered frame.

Reports per-station differences in tau, ceiling and DEM at horizontally
matched stations, plus morphometry aggregates, to analysis_out_v6/
basis_comparison.json. This is the pipeline-level counterpart of the
map-level chain check in chain_flf_to_ed.py.

Run:  PYTHONPATH=src .venv/bin/python analysis/compare_bases.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "analysis_out_ed"    # DLIO basis
B = ROOT / "analysis_out_v6"    # FAST-LIO chained basis (primary)


def load(p):
    rows = list(csv.DictReader(open(p / "roof_thickness.csv")))
    g = lambda k: np.array([float(r[k]) if r[k] else np.nan for r in rows])
    return dict(s=g("s"), x=g("x"), y=g("y"), tau=g("tau"),
                zceil=g("z_ceil"), zdem=g("z_dem"),
                klass=np.array([r["class"] for r in rows]))


def main() -> None:
    a, b = load(A), load(B)
    t = cKDTree(np.c_[a["x"], a["y"]])
    dist, idx = t.query(np.c_[b["x"], b["y"]])
    ok = (dist < 2.0) & (a["klass"][idx] == "intact") & (b["klass"] == "intact")
    dtau = b["tau"][ok] - a["tau"][idx[ok]]
    dceil = b["zceil"][ok] - a["zceil"][idx[ok]]
    sb = b["s"][ok]
    zones = {}
    for lo, hi in ((0, 80), (80, 160), (160, 240), (240, 360)):
        m = (sb >= lo) & (sb < hi)
        if m.sum():
            zones[f"s{lo}-{hi}"] = dict(mean=float(np.mean(dtau[m])),
                                        std=float(np.std(dtau[m])),
                                        n=int(m.sum()))
    ma = json.loads((A / "morphometry_summary.json").read_text())
    mb = json.loads((B / "morphometry_summary.json").read_text())
    out = dict(
        n_matched=int(ok.sum()),
        median_horiz_offset_m=float(np.median(dist[ok])),
        dtau=dict(median=float(np.median(dtau)),
                  rms=float(np.sqrt(np.mean(dtau ** 2))),
                  p5=float(np.percentile(dtau, 5)),
                  p95=float(np.percentile(dtau, 95))),
        dceil=dict(median=float(np.median(dceil)),
                   rms=float(np.sqrt(np.mean(dceil ** 2)))),
        dtau_zones=zones,
        morphometry=dict(
            centreline_m=[ma["total_centreline_m"], mb["total_centreline_m"]],
            area_median=[ma["area_m2"]["median"], mb["area_m2"]["median"]],
            width_median=[ma["width_m"]["median"], mb["width_m"]["median"]],
            height_median=[ma["height_m"]["median"], mb["height_m"]["median"]]),
        note="A = onboard DLIO basis (analysis_out_ed); B = offline FAST-LIO2 "
             "chained basis (analysis_out_v6, primary). Both share the "
             "slice-based registration datum; differences measure estimator "
             "sensitivity of the full pipeline, not the datum.")
    (B / "basis_comparison.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
