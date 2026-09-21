"""Per-station ceiling difference between the two long flights (Fig. 3c).

analysis_out_v9/flight2_repeatability.json stores only the summary of the
flight-to-flight comparison. The new Figure 3c needs the per-station series,
so this script repeats the SAME computation as
analysis/chain_slf_to_ed.py::repeatability (p99 ceiling in a 2.5 m radius
column around every one-metre centreline station, flight 1 minus flight 2,
both maps already in the registered frame) and writes the station table.

It never touches analysis_out_v9: the output goes to
analysis_out_v10/flight2_station_offsets.csv, and the summary recomputed from
the station table is asserted against the archived JSON so the figure cannot
drift from the published numbers.

Run (repo root):
    ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
        .venv/bin/python analysis/flight2_station_offsets.py
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
OUT10 = ROOT / "analysis_out_v10"
F1_PCD = Path(os.environ.get("SLF_F1_PCD", str(ROOT / "maps/flf_30cm_ed_ortho.pcd")))
F2_PCD = Path(os.environ.get("SLF_F2_PCD", str(ROOT / "maps/slf_30cm_ed_ortho.pcd")))
RADIUS = 2.5          # m, column radius, as in chain_slf_to_ed.py
MIN_PTS = 30          # minimum returns in a column, as in chain_slf_to_ed.py


def load(p: Path) -> np.ndarray:
    parts = []
    with BinaryPcdReader(p) as r:
        for c in r.chunks():
            parts.append(c[:, :3].astype(np.float64))
    return np.vstack(parts)


def main() -> None:
    rows = list(csv.DictReader(open(OUT / "centreline.csv")))
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    s = np.array([float(r["s"]) for r in rows])
    f1, f2 = load(F1_PCD), load(F2_PCD)
    print(f"flight 1 {len(f1):,} pts ({F1_PCD.name}), "
          f"flight 2 {len(f2):,} pts ({F2_PCD.name})")
    t1, t2 = cKDTree(f1[:, :2]), cKDTree(f2[:, :2])
    off = np.full(len(x), np.nan)
    z1 = np.full(len(x), np.nan)
    z2 = np.full(len(x), np.nan)
    for i in range(len(x)):
        a = t1.query_ball_point([x[i], y[i]], RADIUS)
        b = t2.query_ball_point([x[i], y[i]], RADIUS)
        if len(a) > MIN_PTS and len(b) > MIN_PTS:
            z1[i] = np.percentile(f1[a, 2], 99)
            z2[i] = np.percentile(f2[b, 2], 99)
            off[i] = z1[i] - z2[i]

    OUT10.mkdir(parents=True, exist_ok=True)
    dst = OUT10 / "flight2_station_offsets.csv"
    with open(dst, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["s", "x", "y", "z_ceil_flight1", "z_ceil_flight2", "dz"])
        for i in range(len(x)):
            w.writerow([f"{s[i]:.4f}", f"{x[i]:.4f}", f"{y[i]:.4f}",
                        "" if np.isnan(z1[i]) else f"{z1[i]:.4f}",
                        "" if np.isnan(z2[i]) else f"{z2[i]:.4f}",
                        "" if np.isnan(off[i]) else f"{off[i]:.4f}"])
    print("wrote", dst)

    ok = ~np.isnan(off)
    got = dict(n_stations=int(ok.sum()),
               median=float(np.nanmedian(off)),
               rms=float(np.sqrt(np.nanmean(off[ok] ** 2))),
               p5=float(np.nanpercentile(off, 5)),
               p95=float(np.nanpercentile(off, 95)))
    ref = json.loads((OUT / "flight2_repeatability.json").read_text())
    print("recomputed", json.dumps(got, indent=2))
    assert got["n_stations"] == ref["n_stations"], (got, ref)
    for k in ("median", "rms", "p5", "p95"):
        assert abs(got[k] - ref[k]) < 1e-6, (k, got[k], ref[k])
    print("summary matches analysis_out_v9/flight2_repeatability.json")


if __name__ == "__main__":
    main()
