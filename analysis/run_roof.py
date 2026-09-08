"""S3: roof-thickness profile, stability ratio, and shielding mass.

tau_k = z_DEM(x_k, y_k) - z_ceil(k): the surface DEM above each centreline
station minus the highest interior ceiling point of the station's cross-section
(both in the aerial frame after the S1 registration). Stations inside a skylight
footprint are classed 'skylight' (tau -> 0 by definition there); everything else
is 'intact'. A station inside the inflated aperture screen is an aperture
station whether or not the surface DEM has a value there: a surface-only DEM
(the v9 orthometric crop) has no cell inside the openings, whereas the earlier
merged cloud filled them with interior points, so from v9 on 'skylight' takes
precedence over 'no_dem' (decision 2026-09-08; no effect on earlier runs, where
every aperture cell had a value). The empirical stability statement is conservative: every intact
surveyed span satisfies tau/L >= kappa_env (the observed minimum), while the
skylight (failed) stations sit at tau ~ 0.

Shielding is reported as areal mass rho * tau (kg/m^2 and g/cm^2) with basalt
rho = 3000 kg/m^3, compared against the terrestrial atmospheric column
(~1033 g/cm^2); no dose-transport modelling is attempted here.

Inputs: analysis_out/centreline.csv, morphometry.csv, aerial_holes.json,
        uncertainty_budget.json, maps/aerial_crop.pcd
Outputs: analysis_out/roof_thickness.csv, roof_summary.json,
         fig_roof_profile.png, fig_roof_footprint.png, dem_grid.npz

Run:  .venv/bin/python analysis/run_roof.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi

from lava_pcd.holes import HoleSet
from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
# Output directory selectable so an alternative registration can be evaluated
# side by side without clobbering the committed baseline (default analysis_out).
import os
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out")))
# DEM source. Default is the cropped aerial surface (aerial_crop.pcd), which
# only covers ~half the surveyed centreline. Set ROOF_DEM_PCD to the
# full-extent surface (e.g. full_surface_and_subsurface_merged.pcd, whose
# surface component is co-registered with aerial_crop to <0.05 m; see
# analysis/merged_dem_validate.py) to extend tau to every station. A distinct
# DEM cache is used per source so the two never clobber each other.
import os
AERIAL = Path(os.environ.get("ROOF_DEM_PCD", str(ROOT / "maps/aerial_crop.pcd")))
DEM_CACHE = OUT / (f"dem_grid_{AERIAL.stem}.npz"
                   if os.environ.get("ROOF_DEM_PCD") else "dem_grid.npz")

DEM_RES = 0.5
RHO_BASALT = 2600.0          # kg/m^3 (assumed bulk, vesicular crust+soil; range 2200-3000)
ATMOSPHERE_G_CM2 = 1033.0    # terrestrial atmospheric column
SKYLIGHT_INFLATE = 1.6
MIN_COVERAGE = 0.5


def build_dem() -> dict:
    """Streamed 0.5 m DEM raster (per-cell max, 3x3 median-filtered)."""
    xmin = ymin = np.inf
    xmax = ymax = -np.inf
    with BinaryPcdReader(AERIAL) as r:
        for c in r.chunks():
            xmin = min(xmin, c[:, 0].min()); xmax = max(xmax, c[:, 0].max())
            ymin = min(ymin, c[:, 1].min()); ymax = max(ymax, c[:, 1].max())
    nx = int(np.ceil((xmax - xmin) / DEM_RES)) + 1
    ny = int(np.ceil((ymax - ymin) / DEM_RES)) + 1
    zmax = np.full((nx, ny), -np.inf, dtype=np.float32)
    with BinaryPcdReader(AERIAL) as r:
        for c in r.chunks():
            ix = ((c[:, 0] - xmin) / DEM_RES).astype(np.int64).clip(0, nx - 1)
            iy = ((c[:, 1] - ymin) / DEM_RES).astype(np.int64).clip(0, ny - 1)
            np.maximum.at(zmax, (ix, iy), c[:, 2])
    filled = zmax > -np.inf
    # Nearest-fill the empty cells, then a 3x3 median to knock out spikes;
    # cells that were empty are restored to NaN afterwards.
    _, (ii, jj) = ndi.distance_transform_edt(~filled, return_indices=True)
    zfill = zmax[ii, jj]
    zmed = ndi.median_filter(zfill, size=3, mode="nearest")
    z = np.where(filled, zmed, np.nan).astype(np.float32)
    return dict(z=z, xmin=float(xmin), ymin=float(ymin), res=DEM_RES)


def dem_at(dem: dict, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    ix = ((x - dem["xmin"]) / dem["res"]).astype(np.int64)
    iy = ((y - dem["ymin"]) / dem["res"]).astype(np.int64)
    z = dem["z"]
    ok = (ix >= 0) & (ix < z.shape[0]) & (iy >= 0) & (iy < z.shape[1])
    out = np.full(len(x), np.nan)
    out[ok] = z[ix[ok], iy[ok]]
    return out


def main() -> None:
    dem_npz = DEM_CACHE
    if dem_npz.exists():
        d = np.load(dem_npz)
        dem = dict(z=d["z"], xmin=float(d["xmin"]), ymin=float(d["ymin"]),
                   res=float(d["res"]))
        print("DEM loaded from cache")
    else:
        print("building DEM ...")
        dem = build_dem()
        np.savez_compressed(dem_npz, **dem)

    rows = list(csv.DictReader(open(OUT / "morphometry.csv")))
    sky = HoleSet.from_json(OUT / "aerial_holes.json")
    budget = json.loads((OUT / "uncertainty_budget.json").read_text())
    sigma_tau = float(budget["sigma_tau_m"])

    # Multipass (ghost) flag: the raw single-flight map contains a vertically
    # displaced duplicate of the skylight stretch (outbound vs return pass
    # without loop closure). Stations whose 4 m column holds a significant
    # fraction of points ABOVE the local surface are internally inconsistent
    # and are excluded from the roof statistics until the loop-closed SLAM
    # re-run (office PC, stage S6.2) replaces this map.
    from scipy.spatial import cKDTree
    from lava_pcd.io.pcd_reader import BinaryPcdReader as _R
    tube_pcd = Path(os.environ.get("ROOF_TUBE_PCD",
                                   str(ROOT / "maps/tube_10cm_aerial.pcd")))
    parts = []
    with _R(tube_pcd) as r:
        for c in r.chunks():
            parts.append(c[:, :3].astype(np.float64))
    tube_pts = np.vstack(parts)
    txy = cKDTree(tube_pts[:, :2])

    s = np.array([float(r["s"]) for r in rows])
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    zc = np.array([float(r["z_ceil"]) if r["z_ceil"] else np.nan for r in rows])
    L = np.array([float(r["width"]) if r["width"] else np.nan for r in rows])
    cov = np.array([float(r["coverage"]) if r["coverage"] else 0.0 for r in rows])

    zdem = dem_at(dem, x, y)
    tau = zdem - zc

    # Skylight class: station horizontally inside an inflated skylight ellipse.
    klass = np.array(["intact"] * len(rows), dtype=object)
    for h in sky.holes:
        cx, cy = h.centroid[0], h.centroid[1]
        a = max(h.semi_major, 1.0) * SKYLIGHT_INFLATE
        d2 = np.hypot(x - cx, y - cy)
        klass[d2 <= a] = "skylight"
    # aperture stations keep their class even where the surface DEM has no
    # cell (surface-only DEM inside an opening); see the docstring
    klass[np.isnan(tau) & (klass != "skylight")] = "no_dem"
    klass[(cov < MIN_COVERAGE) & (klass == "intact")] = "low_coverage"

    ghost_frac = np.zeros(len(rows))
    for i in range(len(rows)):
        if np.isnan(zdem[i]):
            ghost_frac[i] = np.nan   # no surface cell: no column statistic
            continue
        idx = txy.query_ball_point([x[i], y[i]], 4.0)
        if len(idx) < 50:
            continue
        zcol = tube_pts[idx, 2]
        ghost_frac[i] = float(np.mean(zcol > zdem[i] + 1.0))
    klass[(ghost_frac > 0.05) & (klass == "intact")] = "multipass"

    # Physical-consistency screen: a non-positive roof over intact (non-skylight)
    # ground means the mapped ceiling sits at or above the surface DEM, which is
    # impossible for real roof; it signals residual registration/DEM error or an
    # unmapped opening. Exclude such stations from the intact statistics and from
    # kappa_env (roof thickness cannot be negative). Near the skylights the roof
    # genuinely pinches out (tau -> 0); those stations are already 'skylight' or
    # 'multipass'.
    klass[(tau <= 0) & (klass == "intact")] = "inconsistent"

    with open(OUT / "roof_thickness.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["s", "x", "y", "z_dem", "z_ceil", "tau", "sigma_tau",
                    "span", "coverage", "ghost_frac", "class"])
        for i in range(len(rows)):
            w.writerow([f"{s[i]:.1f}", f"{x[i]:.3f}", f"{y[i]:.3f}",
                        f"{zdem[i]:.3f}", f"{zc[i]:.3f}", f"{tau[i]:.3f}",
                        f"{sigma_tau:.3f}", f"{L[i]:.2f}", f"{cov[i]:.2f}",
                        f"{ghost_frac[i]:.3f}", klass[i]])

    it = (klass == "intact")
    ti, Li = tau[it], L[it]
    ratio = ti / Li
    kappa_env = float(np.nanmin(ratio))
    k_arg = int(np.nanargmin(ratio))
    rng = np.random.default_rng(11)
    boot = [np.nanmedian(rng.choice(ti, len(ti))) for _ in range(10_000)]
    shield = ti * RHO_BASALT / 10.0  # g/cm^2 (kg/m^2 -> /10)
    summary = dict(
        n_stations=len(rows), n_intact=int(it.sum()),
        n_skylight=int((klass == "skylight").sum()),
        n_no_dem=int((klass == "no_dem").sum()),
        n_low_coverage=int((klass == "low_coverage").sum()),
        n_multipass=int((klass == "multipass").sum()),
        n_inconsistent=int((klass == "inconsistent").sum()),
        tau_m=dict(min=round(float(np.nanmin(ti)), 2),
                   max=round(float(np.nanmax(ti)), 2),
                   median=round(float(np.nanmedian(ti)), 2),
                   q25=round(float(np.nanpercentile(ti, 25)), 2),
                   q75=round(float(np.nanpercentile(ti, 75)), 2),
                   median_ci95=[round(float(np.percentile(boot, 2.5)), 2),
                                round(float(np.percentile(boot, 97.5)), 2)],
                   sigma=sigma_tau),
        kappa_env=round(kappa_env, 4),
        kappa_env_station=dict(s=float(s[it][k_arg]),
                               tau=round(float(ti[k_arg]), 2),
                               span=round(float(Li[k_arg]), 2)),
        shielding_g_cm2=dict(median=round(float(np.nanmedian(shield)), 0),
                             min=round(float(np.nanmin(shield)), 0),
                             atmosphere_ratio_median=round(
                                 float(np.nanmedian(shield)) / ATMOSPHERE_G_CM2, 2)),
        rho_basalt=RHO_BASALT,
    )
    (OUT / "roof_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    # Profile figure (QA; publication styling in S5).
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    axes[0].plot(s[it], ti, ".", ms=3)
    axes[0].fill_between(s[it], ti - sigma_tau, ti + sigma_tau, alpha=0.2)
    axes[0].set_ylabel(r"roof thickness $\tau$ (m)")
    for h in sky.holes:
        d2 = np.hypot(x - h.centroid[0], y - h.centroid[1])
        if (d2 < 30).any():
            axes[0].axvline(s[np.argmin(d2)], color="orange", ls="--", lw=0.8)
    axes[1].plot(s, L, ".", ms=3, color="tab:green")
    axes[1].set_ylabel("span L (m)")
    axes[2].plot(s[it], ratio, ".", ms=3, color="tab:red")
    axes[2].axhline(kappa_env, color="k", ls=":", lw=1)
    axes[2].set_ylabel(r"$\tau/L$")
    axes[2].set_xlabel("distance along tube (m)")
    fig.savefig(OUT / "fig_roof_profile.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 8))
    sc = ax.scatter(x[it], y[it], c=ti, s=8, cmap="magma")
    ax.scatter(x[~it], y[~it], c="cyan", s=4, label="skylight/no-DEM")
    plt.colorbar(sc, label=r"$\tau$ (m)")
    ax.set_aspect("equal")
    ax.set_title("roof thickness on the surface footprint")
    fig.savefig(OUT / "fig_roof_footprint.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
