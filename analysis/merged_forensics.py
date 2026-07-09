"""T3 forensics: characterize maps/full_surface_and_subsurface_merged.pcd.

Questions (plan/06 T3):
  1. Extent: does the merged cloud cover more of the tube than aerial_crop?
  2. Coverage: how many centreline stations gain a DEM value under the merged
     cloud vs aerial_crop (the v1 173-no-DEM wall)?
  3. Ghost: is the +5-8 m vertical duplicate of the skylight stretch still
     present in the merged cloud's subsurface component?

Streams the 3 GB cloud in chunks; nothing is held in full. Writes
analysis_out/merged_forensics.json and two QA PNGs.

Run:  env -u PYTHONPATH .venv/bin/python analysis/merged_forensics.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
MERGED = ROOT / "maps/full_surface_and_subsurface_merged.pcd"
AERIAL = ROOT / "maps/aerial_crop.pcd"
RES = 1.0   # coarse raster for extent/coverage maps


def bbox_and_raster(path: Path):
    """Two streaming passes: bbox, then per-cell max-z / min-z / count."""
    xmin = ymin = zmin = np.inf
    xmax = ymax = zmax = -np.inf
    n = 0
    with BinaryPcdReader(path) as r:
        for c in r.chunks():
            xmin = min(xmin, c[:, 0].min()); xmax = max(xmax, c[:, 0].max())
            ymin = min(ymin, c[:, 1].min()); ymax = max(ymax, c[:, 1].max())
            zmin = min(zmin, c[:, 2].min()); zmax = max(zmax, c[:, 2].max())
            n += len(c)
    nx = int(np.ceil((xmax - xmin) / RES)) + 1
    ny = int(np.ceil((ymax - ymin) / RES)) + 1
    zhi = np.full((nx, ny), -np.inf, dtype=np.float32)
    zlo = np.full((nx, ny), np.inf, dtype=np.float32)
    cnt = np.zeros((nx, ny), dtype=np.int64)
    with BinaryPcdReader(path) as r:
        for c in r.chunks():
            ix = ((c[:, 0] - xmin) / RES).astype(np.int64).clip(0, nx - 1)
            iy = ((c[:, 1] - ymin) / RES).astype(np.int64).clip(0, ny - 1)
            np.maximum.at(zhi, (ix, iy), c[:, 2])
            np.minimum.at(zlo, (ix, iy), c[:, 2])
            np.add.at(cnt, (ix, iy), 1)
    return dict(n=n, xmin=float(xmin), ymin=float(ymin), xmax=float(xmax),
                ymax=float(ymax), zmin=float(zmin), zmax=float(zmax),
                nx=nx, ny=ny, zhi=zhi, zlo=zlo, cnt=cnt)


def coverage_at_stations(rast, stations):
    ix = ((stations[:, 0] - rast["xmin"]) / RES).astype(np.int64)
    iy = ((stations[:, 1] - rast["ymin"]) / RES).astype(np.int64)
    ok = (ix >= 0) & (ix < rast["nx"]) & (iy >= 0) & (iy < rast["ny"])
    hit = np.zeros(len(stations), dtype=bool)
    hit[ok] = rast["cnt"][ix[ok], iy[ok]] > 0
    return hit


def main():
    st = []
    with open(OUT / "centreline.csv") as f:
        for row in csv.DictReader(f):
            st.append([float(row["x"]), float(row["y"]), float(row["z"])])
    stations = np.array(st)

    print("[1/2] rasterizing merged cloud (3 GB, streamed) ...")
    m = bbox_and_raster(MERGED)
    print("[2/2] rasterizing aerial_crop ...")
    a = bbox_and_raster(AERIAL)

    hit_m = coverage_at_stations(m, stations)
    hit_a = coverage_at_stations(a, stations)

    # Ghost test: over the tube footprint (cells the aerial_crop covers AND the
    # merged cloud has a tall vertical span), how many cells show a z-range
    # far larger than a bare photogrammetric surface would (surface alone is
    # ~2.5D, range < a few m per cell). A large per-cell z-range means the
    # subsurface tube column is present; a *bimodal* column with mass well
    # above the local surface is the multipass ghost.
    span = m["zhi"] - m["zlo"]
    tall = (m["cnt"] > 0) & (span > 5.0)
    result = dict(
        merged=dict(n_points=m["n"],
                    x=[round(m["xmin"], 1), round(m["xmax"], 1)],
                    y=[round(m["ymin"], 1), round(m["ymax"], 1)],
                    z=[round(m["zmin"], 1), round(m["zmax"], 1)],
                    extent_x_m=round(m["xmax"] - m["xmin"], 1),
                    extent_y_m=round(m["ymax"] - m["ymin"], 1)),
        aerial_crop=dict(n_points=a["n"],
                         x=[round(a["xmin"], 1), round(a["xmax"], 1)],
                         y=[round(a["ymin"], 1), round(a["ymax"], 1)],
                         extent_x_m=round(a["xmax"] - a["xmin"], 1),
                         extent_y_m=round(a["ymax"] - a["ymin"], 1)),
        stations_total=len(stations),
        stations_with_merged_coverage=int(hit_m.sum()),
        stations_with_aerial_crop_coverage=int(hit_a.sum()),
        stations_gained=int((hit_m & ~hit_a).sum()),
        n_cells_tall_span_gt5m=int(tall.sum()),
    )
    (OUT / "merged_forensics.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))

    # QA figure: top-down coverage, aerial_crop vs merged, centreline overlaid.
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for ax, rast, title in ((axes[0], a, "aerial_crop cell coverage"),
                            (axes[1], m, "merged cloud cell coverage")):
        occ = (rast["cnt"] > 0).astype(float)
        ax.imshow(occ.T, origin="lower",
                  extent=[rast["xmin"], rast["xmax"], rast["ymin"], rast["ymax"]],
                  cmap="Greys", vmax=1.4, aspect="equal")
        ax.plot(stations[:, 0], stations[:, 1], "r-", lw=1.2, label="centreline")
        ax.set_title(title)
        ax.legend(loc="upper right")
    fig.savefig(OUT / "fig_merged_coverage.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # QA figure: per-cell z-range over merged cloud (ghost/subsurface map).
    fig, ax = plt.subplots(figsize=(12, 7))
    im = ax.imshow(np.where(m["cnt"].T > 0, span.T, np.nan), origin="lower",
                   extent=[m["xmin"], m["xmax"], m["ymin"], m["ymax"]],
                   cmap="magma", aspect="equal", vmax=40)
    ax.plot(stations[:, 0], stations[:, 1], "c-", lw=1.0)
    plt.colorbar(im, label="per-cell z-range (m)")
    ax.set_title("merged cloud vertical span per 1 m cell "
                 "(tall = subsurface column present)")
    fig.savefig(OUT / "fig_merged_zspan.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote merged_forensics.json, fig_merged_coverage.png, "
          "fig_merged_zspan.png")


if __name__ == "__main__":
    main()
