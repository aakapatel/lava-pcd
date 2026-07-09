"""T4 acceptance test: is the +5-8 m double-ceiling ghost present in a FAST-LIO
scans.pcd?

The Mac single-flight map had the skylight stretch duplicated ~5-8 m vertically
(outbound vs return pass, no loop closure). The accumulated FAST-LIO map carries
no per-point time, so we test the map geometry directly: grid the XY plane, build
a per-cell vertical histogram, isolate the ceiling band (top of each column), and
flag cells whose ceiling band is split into two peaks separated by > GAP metres.
A clean (loop-closed / drift-free) map shows a single ceiling mode almost
everywhere; a ghosted map shows a broad population of doubled-ceiling cells.

LIMITATION: this is a purely geometric pre-check. It only tests columns tall
enough (span > 12 m) to separate a doubled ceiling from an ordinary floor+ceiling
pair, and in these tube maps only ~1-4% of populated columns qualify, so a low
count here is encouraging but NOT a final acceptance. The definitive ghost test
is DEM-based (does the mapped ceiling rise above the co-registered surface DEM?
run_roof.py multipass flag) and requires registering the map to the aerial frame
first (T5). An earlier, cruder version of this script that searched the top 8 m
of EVERY column reported ~44% 'doubled' cells; inspection showed those were
floor+ceiling pairs of short passages at hover points, not the ghost, hence the
tall-column + floor-below refinement here.

Assumes +z is up (lidar/FAST-LIO convention). Streams the 2.4 GB cloud twice.

Run:  env -u PYTHONPATH .venv/bin/python analysis/fastlio_ghost_diag.py <scans.pcd>
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
XY_RES = 2.0        # horizontal cell size (m)
Z_RES = 0.25        # vertical histogram bin (m)
GAP = 3.0           # min separation to call two ceiling modes distinct (m)
CEIL_BAND = 8.0     # search the top this-many-m of each column for a doubled ceiling
MIN_PTS = 60        # ignore sparse cells


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path.home() / "rosbags/LavaTube_Iceland/April2025/Iceland_Apr_25/day2/"
        "fastlio_maps/second_long_flight/scans.pcd")
    tag = path.parent.name

    # pass 1: bbox
    xmin = ymin = zmin = np.inf
    xmax = ymax = zmax = -np.inf
    n = 0
    with BinaryPcdReader(path) as r:
        for c in r.chunks():
            xmin = min(xmin, c[:, 0].min()); xmax = max(xmax, c[:, 0].max())
            ymin = min(ymin, c[:, 1].min()); ymax = max(ymax, c[:, 1].max())
            zmin = min(zmin, c[:, 2].min()); zmax = max(zmax, c[:, 2].max())
            n += len(c)
    nx = int(np.ceil((xmax - xmin) / XY_RES)) + 1
    ny = int(np.ceil((ymax - ymin) / XY_RES)) + 1
    nz = int(np.ceil((zmax - zmin) / Z_RES)) + 1
    print(f"{tag}: {n:,} pts  x[{xmin:.1f},{xmax:.1f}] y[{ymin:.1f},{ymax:.1f}] "
          f"z[{zmin:.1f},{zmax:.1f}]  grid {nx}x{ny}x{nz}")

    # pass 2: per-cell vertical histogram (uint16 counts)
    hist = np.zeros((nx * ny, nz), dtype=np.uint16)
    with BinaryPcdReader(path) as r:
        for c in r.chunks():
            ix = ((c[:, 0] - xmin) / XY_RES).astype(np.int64).clip(0, nx - 1)
            iy = ((c[:, 1] - ymin) / XY_RES).astype(np.int64).clip(0, ny - 1)
            iz = ((c[:, 2] - zmin) / Z_RES).astype(np.int64).clip(0, nz - 1)
            flat = ix * ny + iy
            np.add.at(hist, (flat, iz), 1)

    zc = zmin + (np.arange(nz) + 0.5) * Z_RES
    tot = hist.sum(1)
    cells = np.where(tot >= MIN_PTS)[0]

    # A doubled CEILING (the ghost) is two peaks near the TOP of a column that
    # ALSO has a distinct floor band well below them; a plain floor+ceiling pair
    # in a short (thin) passage is NOT a ghost. So we require the column to be
    # tall enough (full span > MIN_SPAN) that the true floor is clearly separated
    # from the two ceiling peaks, and require real point mass below the lower of
    # the two ceiling peaks (the floor).
    MIN_SPAN = 12.0        # only tall columns can exhibit a *doubled ceiling*
    FLOOR_GAP = 5.0        # floor must sit at least this far below the lower peak
    n_double = 0
    n_tall = 0
    gaps = []
    dbl_xy = []
    for cell in cells:
        h = hist[cell].astype(float)
        occ = np.where(h > 0)[0]
        if len(occ) == 0:
            continue
        full_span = zc[occ.max()] - zc[occ.min()]
        if full_span < MIN_SPAN:
            continue                      # thin passage: floor+ceiling only
        n_tall += 1
        top_z = zc[occ.max()]
        band = (zc >= top_z - CEIL_BAND) & (zc <= top_z + 0.01)
        hb = h.copy(); hb[~band] = 0
        if hb.sum() < MIN_PTS * 0.3:
            continue
        k = np.array([0.25, 0.5, 0.25])
        hs = np.convolve(hb, k, mode="same")
        thr = 0.15 * hs.max()
        peaks, _ = find_peaks(hs, height=thr, distance=int(GAP / Z_RES))
        if len(peaks) >= 2:
            pk_z = zc[peaks]
            sep = pk_z.max() - pk_z.min()
            lower_peak = pk_z.min()
            # floor = point mass at least FLOOR_GAP below the lower ceiling peak
            floor_mass = h[zc <= lower_peak - FLOOR_GAP].sum()
            if sep >= GAP and floor_mass >= MIN_PTS * 0.3:
                n_double += 1
                gaps.append(sep)
                ci = cell // ny; cj = cell % ny
                dbl_xy.append((xmin + ci * XY_RES, ymin + cj * XY_RES))

    frac = n_double / max(n_tall, 1)
    print(f"cells with >= {MIN_PTS} pts: {len(cells)}; "
          f"tall enough (span > {MIN_SPAN} m) to test: {n_tall}")
    print(f"doubled-ceiling cells (>= {GAP} m split, floor below): {n_double} "
          f"({100 * frac:.1f}% of tall columns)")
    if gaps:
        print(f"  ceiling-split separation: median {np.median(gaps):.1f} m, "
              f"p90 {np.percentile(gaps, 90):.1f} m, max {np.max(gaps):.1f} m")

    # QA: top-down map of doubled-ceiling cells
    fig, ax = plt.subplots(figsize=(11, 9))
    occ_cells = [( (cl // ny) * XY_RES + xmin, (cl % ny) * XY_RES + ymin)
                 for cl in cells]
    ox, oy = zip(*occ_cells)
    ax.scatter(ox, oy, s=3, c="0.8", label="occupied")
    if dbl_xy:
        dx, dy = zip(*dbl_xy)
        ax.scatter(dx, dy, s=8, c="red", label="doubled ceiling")
    ax.set_aspect("equal"); ax.legend(); ax.set_title(
        f"{tag}: doubled-ceiling cells = {n_double} of {len(cells)} "
        f"({100*frac:.1f}%)")
    fig.savefig(OUT / f"fig_ghost_{tag}.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote fig_ghost_{tag}.png")


if __name__ == "__main__":
    main()
