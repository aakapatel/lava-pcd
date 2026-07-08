"""S2: centreline (L1-medial skeleton) and along-tube morphometry.

Implements the L1-median skeleton of Huang et al. 2013 (SIGGRAPH), as applied
to lava tubes by Yang et al. 2024 (Int. J. Appl. Earth Obs.): sample points are
contracted onto the local L1-median of the shell cloud under a repulsion term,
with a growing support radius, then ordered into a polyline by the diameter
path of a minimum-spanning-tree, spline-smoothed and resampled at 1 m.

At each station the perpendicular cross-section is measured from the 10 cm
cloud by angular-binned radial profiles (72 bins, per-bin median radius), which
tolerate the incomplete shell (lidar shadowing) and report their own coverage.

Inputs (aerial = gravity-referenced frame, made by `lava-pcd transform`):
    maps/tube_30cm_aerial.pcd   skeleton input
    maps/tube_10cm_aerial.pcd   cross-section input
Outputs:
    analysis_out/centreline.csv          s, x, y, z, tx, ty, tz
    analysis_out/morphometry.csv         per-station section metrics
    analysis_out/morphometry_summary.json  aggregates for the paper
    analysis_out/fig_skeleton_topdown.png / fig_skeleton_side.png
    analysis_out/fig_sections_gallery.png  QA: 12 random sections

Run:  .venv/bin/python analysis/run_morphometry.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.sparse.csgraph import minimum_spanning_tree, shortest_path
from scipy.spatial import cKDTree

from lava_pcd.io.pcd_reader import BinaryPcdReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_out"
RNG = np.random.default_rng(7)

SKEL_INPUT = ROOT / "maps/tube_30cm_aerial.pcd"
SECT_INPUT = ROOT / "maps/tube_10cm_aerial.pcd"

# L1-skeleton parameters
N_Q = 60_000          # shell points used as the attraction set
X_SPACING = 2.0       # initial sample spacing (m)
H_LEVELS = (8.0, 12.0, 16.0)   # growing support radii (m); tube spans ~4-25 m
MU = 0.35             # repulsion weight (Huang et al. use 0.35)
ITERS = 40            # iterations per level
STEP_TOL = 1e-3       # convergence: max move (m)

# Cross-section parameters
SLAB = 0.75           # half-thickness of the station slab (m)
R_MAX = 30.0          # max radius searched around the centreline (m)
N_BINS = 72           # angular bins (5 degrees)
MIN_COVERAGE = 0.5    # report area stats only where >=50% of bins are hit
STATION = 1.0         # station spacing (m)


def load_xyz(path: Path, keep: int | None = None) -> np.ndarray:
    parts = []
    with BinaryPcdReader(path) as r:
        for chunk in r.chunks():
            parts.append(chunk[:, :3].astype(np.float64))
    pts = np.vstack(parts)
    if keep is not None and len(pts) > keep:
        pts = pts[RNG.choice(len(pts), keep, replace=False)]
    return pts


def l1_skeleton(Q: np.ndarray) -> np.ndarray:
    """Contract samples X onto the local L1-median axis of Q (Huang 2013)."""
    # Initial samples: grid-decimate Q at X_SPACING.
    key = np.floor(Q / X_SPACING).astype(np.int64)
    _, idx = np.unique(key, axis=0, return_index=True)
    X = Q[np.sort(idx)].copy()
    print(f"  skeleton: |Q|={len(Q)}, |X|={len(X)}")

    qtree = cKDTree(Q)
    for h in H_LEVELS:
        theta_den = (h / 2.0) ** 2
        for it in range(ITERS):
            xtree = cKDTree(X)
            # Attraction: neighbours in Q within h.
            nbrs = xtree.query_ball_tree(qtree, h)
            # Repulsion: neighbours in X within h.
            xnbrs = xtree.query_ball_tree(xtree, h)
            X_new = X.copy()
            for i in range(len(X)):
                qi = Q[nbrs[i]]
                if len(qi) == 0:
                    continue
                d = qi - X[i]
                r = np.linalg.norm(d, axis=1)
                r = np.maximum(r, 1e-6)
                a = np.exp(-(r ** 2) / theta_den) / r
                mean_term = (qi * a[:, None]).sum(0) / a.sum()

                xi = X[np.setdiff1d(xnbrs[i], [i])]
                if len(xi):
                    dxr = X[i] - xi
                    rr = np.maximum(np.linalg.norm(dxr, axis=1), 1e-6)
                    b = np.exp(-(rr ** 2) / theta_den) / (rr ** 2)
                    # directionality (sigma) from the X-neighbourhood covariance
                    C = np.cov((xi - X[i]).T) if len(xi) > 2 else np.eye(3)
                    w = np.linalg.eigvalsh(C)
                    sigma = w[-1] / max(w.sum(), 1e-9)
                    rep = (dxr * b[:, None]).sum(0) / b.sum()
                    X_new[i] = mean_term + MU * sigma * rep
                else:
                    X_new[i] = mean_term
            move = float(np.max(np.linalg.norm(X_new - X, axis=1)))
            X = X_new
            if move < STEP_TOL:
                break
        print(f"  h={h}: {it + 1} iters, last max move {move:.4f} m")
    return X


def order_polyline(X: np.ndarray) -> np.ndarray:
    """Order contracted points along the tube: MST diameter path."""
    # Collapse near-duplicates (contracted clusters) at 1 m.
    key = np.floor(X / 1.0).astype(np.int64)
    _, idx = np.unique(key, axis=0, return_index=True)
    P = X[np.sort(idx)]
    tree = cKDTree(P)
    dist = tree.sparse_distance_matrix(tree, max_distance=8.0)
    mst = minimum_spanning_tree(dist.tocsr())
    # Tree diameter: farthest node from node 0, then farthest from that.
    d0 = shortest_path(mst, directed=False, indices=0)
    a = int(np.nanargmax(np.where(np.isinf(d0), np.nan, d0)))
    da, pred = shortest_path(mst, directed=False, indices=a, return_predecessors=True)
    b = int(np.nanargmax(np.where(np.isinf(da), np.nan, da)))
    path = [b]
    while path[-1] != a and pred[path[-1]] >= 0:
        path.append(int(pred[path[-1]]))
    print(f"  polyline: {len(P)} nodes, diameter path {len(path)} nodes, "
          f"length {da[b]:.1f} m")
    return P[path]


def smooth_resample(poly: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Cubic smoothing spline through the polyline, resampled at 1 m."""
    seg = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    u = np.concatenate([[0.0], np.cumsum(seg)])
    tck, _ = interpolate.splprep(poly.T, u=u, s=len(poly) * 0.5, k=3)
    total = u[-1]
    s = np.arange(0.0, total, STATION)
    pts = np.array(interpolate.splev(s, tck)).T
    tang = np.array(interpolate.splev(s, tck, der=1)).T
    tang /= np.linalg.norm(tang, axis=1, keepdims=True)
    return pts, tang


def section_metrics(P: np.ndarray, station: np.ndarray, T: np.ndarray,
                    tree: cKDTree, cloud: np.ndarray) -> dict | None:
    """Angular-binned radial section in the plane normal to T at `station`."""
    idx = tree.query_ball_point(station, R_MAX)
    if len(idx) < 30:
        return None
    pts = cloud[idx]
    d = pts - station
    along = d @ T
    m = np.abs(along) <= SLAB
    if m.sum() < 30:
        return None
    sec = d[m]
    # In-plane basis: u horizontal, v completes (approximately vertical).
    uvec = np.cross(T, [0.0, 0.0, 1.0])
    nu = np.linalg.norm(uvec)
    if nu < 1e-3:
        uvec = np.array([1.0, 0.0, 0.0])
    else:
        uvec = uvec / nu
    vvec = np.cross(T, uvec)
    if vvec[2] < 0:
        vvec = -vvec
    su, sv = sec @ uvec, sec @ vvec
    # Recentre on the section median (the centreline may be slightly off-axis).
    cu, cv = np.median(su), np.median(sv)
    du, dv = su - cu, sv - cv
    r = np.hypot(du, dv)
    ang = np.arctan2(dv, du)
    bins = ((ang + np.pi) / (2 * np.pi) * N_BINS).astype(int).clip(0, N_BINS - 1)
    # Per bin, the section boundary is the *innermost dense cluster* of radii,
    # not the median of all returns: through a skylight the lidar also mapped
    # patches of the outside terrain, which appear as a second, farther layer
    # along upward rays and would otherwise inflate the ceiling and the area.
    rmed = np.full(N_BINS, np.nan)
    zmax_in = np.full(N_BINS, np.nan)   # top of the innermost layer (world v)
    for b in range(N_BINS):
        mb = bins == b
        rb = r[mb]
        if len(rb) < 2:
            continue
        order = np.argsort(rb)
        rs = rb[order]
        gaps = np.where(np.diff(rs) > 1.5)[0]
        stop = len(rs) if len(gaps) == 0 else None
        if stop is None:
            for g in gaps:
                if g + 1 >= 3:          # innermost cluster with >=3 points
                    stop = g + 1
                    break
            if stop is None:
                stop = len(rs)
        rmed[b] = np.median(rs[:stop])
        zmax_in[b] = float(np.max(dv[mb][order][:stop]))
    covered = int(np.isfinite(rmed).sum())
    coverage = covered / N_BINS
    area = float(np.nansum(0.5 * rmed ** 2 * (2 * np.pi / N_BINS)))
    # Robust extents: 1st-99th percentiles reject stray returns (e.g. dust
    # points metres below the floor).
    width = float(np.percentile(du, 99) - np.percentile(du, 1))
    height = float(np.percentile(dv, 99) - np.percentile(dv, 1))
    # Aspect ratio from the section's second moments.
    M = np.cov(np.vstack([du, dv]))
    ev = np.linalg.eigvalsh(M)
    eta = float(np.sqrt(max(ev[1], 1e-9) / max(ev[0], 1e-9)))
    # Interior ceiling: top of the innermost layer over near-vertical bins
    # (bin direction within ~45 degrees of up). z of the section frame is
    # station z + v (vvec is unit and near-vertical for a subhorizontal tube).
    bin_dirs = (np.arange(N_BINS) + 0.5) / N_BINS * 2 * np.pi - np.pi
    up_bins = np.sin(bin_dirs) > 0.7
    vceil = np.nanmax(np.where(up_bins, zmax_in, np.nan))
    vfloor = np.nanmin(zmax_in)  # conservative floor from innermost layer
    zpts = pts[m][:, 2]
    z_ceil = float(station[2] + cv + vceil) if np.isfinite(vceil) else float(zpts.max())
    return dict(area=area, coverage=coverage, width=width, height=height,
                eta=eta, d_eq=float(2 * np.sqrt(area / np.pi)),
                z_floor=float(np.percentile(zpts, 1)), z_ceil=z_ceil,
                n_pts=int(m.sum()))


def main() -> None:
    print("[1/4] loading clouds ...")
    Q = load_xyz(SKEL_INPUT, keep=N_Q)
    cloud = load_xyz(SECT_INPUT)
    tree = cKDTree(cloud)
    print(f"  sections cloud: {len(cloud)} pts")

    print("[2/4] L1-medial skeleton ...")
    cache = OUT / "skeleton_poly.npz"
    if cache.exists():
        poly = np.load(cache)["poly"]
        print(f"  polyline loaded from cache ({len(poly)} nodes)")
    else:
        X = l1_skeleton(Q)
        poly = order_polyline(X)
        np.savez_compressed(cache, poly=poly)
    stations, tangents = smooth_resample(poly)
    # Orient s=0 at the entrance end: the end nearest the skylights (the three
    # openings sit in the starting section of the surveyed tube).
    sky_y = 640.0  # mean skylight northing in the aerial local frame
    if abs(stations[0][1] - sky_y) > abs(stations[-1][1] - sky_y):
        stations, tangents = stations[::-1].copy(), -tangents[::-1].copy()
        print("  centreline flipped: s=0 at the entrance/skylight end")
    total_len = STATION * (len(stations) - 1)
    print(f"  centreline: {len(stations)} stations, {total_len:.1f} m")

    print("[3/4] cross-sections ...")
    rows = []
    for k, (sk, tk) in enumerate(zip(stations, tangents)):
        met = section_metrics(poly, sk, tk, tree, cloud)
        row = dict(s=k * STATION, x=sk[0], y=sk[1], z=sk[2],
                   tx=tk[0], ty=tk[1], tz=tk[2])
        if met:
            row.update(met)
        rows.append(row)

    import csv
    keys = ["s", "x", "y", "z", "tx", "ty", "tz", "area", "coverage", "width",
            "height", "eta", "d_eq", "z_floor", "z_ceil", "n_pts"]
    with open(OUT / "morphometry.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow({k: (f"{row[k]:.4f}" if isinstance(row.get(k), float)
                            else row.get(k, "")) for k in keys})
    with open(OUT / "centreline.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["s", "x", "y", "z", "tx", "ty", "tz"])
        for k, (sk, tk) in enumerate(zip(stations, tangents)):
            w.writerow([k * STATION, *[f"{v:.4f}" for v in (*sk, *tk)]])

    # Sinuosity over 50 m windows.
    W = 50
    sin_vals = []
    for k in range(len(stations)):
        lo, hi = max(0, k - W // 2), min(len(stations) - 1, k + W // 2)
        if hi - lo < W // 2:
            sin_vals.append(np.nan)
            continue
        pathlen = STATION * (hi - lo)
        chord = float(np.linalg.norm(stations[hi] - stations[lo]))
        sin_vals.append(pathlen / max(chord, 1e-6))
    sin_vals = np.array(sin_vals)

    # Stats exclude the first/last 3 stations (truncated rings at the map ends).
    ok = [r for r in rows[3:-3] if r.get("coverage", 0) >= MIN_COVERAGE]
    areas = np.array([r["area"] for r in ok])
    heights = np.array([r["height"] for r in ok])
    etas = np.array([r["eta"] for r in ok])
    boot = RNG.choice(areas, (10_000, len(areas))).mean(1) if len(areas) else []
    summary = dict(
        total_centreline_m=round(total_len, 1),
        n_stations=len(stations),
        n_valid_sections=len(ok),
        min_coverage=MIN_COVERAGE,
        area_m2=dict(min=round(float(areas.min()), 1),
                     max=round(float(areas.max()), 1),
                     median=round(float(np.median(areas)), 1),
                     q25=round(float(np.percentile(areas, 25)), 1),
                     q75=round(float(np.percentile(areas, 75)), 1),
                     mean_boot_ci95=[round(float(np.percentile(boot, 2.5)), 1),
                                     round(float(np.percentile(boot, 97.5)), 1)]),
        height_m=dict(min=round(float(heights.min()), 2),
                      median=round(float(np.median(heights)), 2),
                      max=round(float(heights.max()), 2)),
        eta=dict(median=round(float(np.median(etas)), 2),
                 q25=round(float(np.percentile(etas, 25)), 2),
                 q75=round(float(np.percentile(etas, 75)), 2)),
        width_m=dict(min=round(float(min(r["width"] for r in ok)), 2),
                     median=round(float(np.median([r["width"] for r in ok])), 2),
                     max=round(float(max(r["width"] for r in ok)), 2)),
        width_min_m=round(float(min(r["width"] for r in ok)), 2),
        sinuosity_50m=dict(mean=round(float(np.nanmean(sin_vals)), 3),
                           max=round(float(np.nanmax(sin_vals)), 3)),
    )
    (OUT / "morphometry_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    print("[4/4] renders ...")
    sub = cloud[RNG.choice(len(cloud), min(400_000, len(cloud)), replace=False)]
    for view, (i, j), name in ((("top"), (0, 1), "fig_skeleton_topdown.png"),
                               (("side"), (1, 2), "fig_skeleton_side.png")):
        fig, ax = plt.subplots(figsize=(14, 7))
        ax.scatter(sub[:, i], sub[:, j], s=0.3, c=sub[:, 2], cmap="viridis",
                   alpha=0.35, linewidths=0)
        ax.plot(stations[:, i], stations[:, j], "r-", lw=1.6, label="centreline")
        ax.set_aspect("equal")
        ax.legend()
        ax.set_title(f"L1-medial centreline, {view} view")
        fig.savefig(OUT / name, dpi=140, bbox_inches="tight")
        plt.close(fig)

    picks = [r for r in rows if r.get("coverage", 0) >= MIN_COVERAGE]
    picks = picks[:: max(1, len(picks) // 12)][:12]
    fig, axes = plt.subplots(3, 4, figsize=(16, 10))
    for ax, row in zip(axes.ravel(), picks):
        k = int(row["s"] / STATION)
        sk, tk = stations[k], tangents[k]
        idx = tree.query_ball_point(sk, R_MAX)
        pts = cloud[idx] - sk
        m = np.abs(pts @ tk) <= SLAB
        uvec = np.cross(tk, [0, 0, 1.0])
        uvec = uvec / max(np.linalg.norm(uvec), 1e-6)
        vvec = np.cross(tk, uvec)
        if vvec[2] < 0:
            vvec = -vvec
        ax.scatter(pts[m] @ uvec, pts[m] @ vvec, s=1.5, linewidths=0)
        ax.set_title(f"s={row['s']:.0f} m  A={row['area']:.0f} m$^2$ "
                     f"cov={row['coverage']:.2f}", fontsize=8)
        ax.set_aspect("equal")
    fig.savefig(OUT / "fig_sections_gallery.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("done")


if __name__ == "__main__":
    main()
