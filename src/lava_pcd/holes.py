"""Skylight ("hole") detection in aerial and lava-tube ``.pcd`` clouds.

A skylight is a collapse/opening that is a *void* in the aerial ground and a
*hole* in the tube ceiling. Detecting the set of skylights in each cloud yields a
sparse landmark constellation that :mod:`lava_pcd.register` matches to recover the
rigid transform between the two clouds (their only shared geometry).

Two detection modes, both reduced to "find enclosed anomalies in a 2-D grid taken
looking along the up-axis":

* ``aerial`` -- up is ``+Z`` (gravity-aligned survey). A skylight is an *enclosed
  empty region* of a top-down occupancy grid (laser passes through, no return).
* ``ceiling`` -- the tube is in an arbitrary SLAM frame, so an up-axis is
  estimated (or given via ``up``) and the cloud is rotated so up -> ``+Z``. The
  ceiling-height field (per-cell max Z) is built; a skylight is an enclosed cell
  whose ceiling height deviates strongly from the local ceiling (the roof is
  missing there -- the cell sees the sky above or the floor below) or has no
  return at all.

Detection is the fragile step, so results are meant to be reviewed (``review_holes``)
or placed by hand (``pick_holes``) before registration. Grids are accumulated by
streaming the cloud in chunks, so memory stays bounded by the *grid* size, not the
point count.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from lava_pcd.geometry import normalize, rotation_align
from lava_pcd.io.pcd_reader import BinaryPcdReader

DEFAULT_RESOLUTION = 1.0          # grid cell size in coordinate units (metres)
DEFAULT_MIN_AREA = 4.0            # ignore voids smaller than this (m^2)
DEFAULT_CEILING_JUMP = 1.0        # ceiling-height deviation flagged as "open" (m)
_MAX_GRID_CELLS = 50_000_000      # guard against an absurdly fine grid
_UP_SAMPLE = 200_000             # points sampled for up-axis estimation


@dataclass
class Hole:
    """A single detected skylight, in the cloud's local frame."""

    id: int
    centroid: tuple[float, float, float]
    radius: float
    area: float
    n_cells: int


@dataclass
class HoleSet:
    """All skylights found in one cloud, plus the frame they were found in."""

    holes: list[Hole]
    origin: tuple[float, float, float]
    up: tuple[float, float, float]
    mode: str
    resolution: float
    source_path: str = ""

    def centroids(self) -> np.ndarray:
        """``(M, 3)`` array of hole centroids (empty -> shape ``(0, 3)``)."""
        if not self.holes:
            return np.empty((0, 3), dtype=np.float64)
        return np.array([h.centroid for h in self.holes], dtype=np.float64)

    def radii(self) -> np.ndarray:
        """``(M,)`` array of equivalent radii."""
        return np.array([h.radius for h in self.holes], dtype=np.float64)

    def to_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "origin": list(self.origin),
            "up": list(self.up),
            "mode": self.mode,
            "resolution": self.resolution,
            "source_path": self.source_path,
            "holes": [asdict(h) for h in self.holes],
        }
        path.write_text(json.dumps(payload, indent=2))
        return path

    @classmethod
    def from_json(cls, path: str | Path) -> "HoleSet":
        data = json.loads(Path(path).read_text())
        holes = [
            Hole(
                id=int(h["id"]),
                centroid=tuple(h["centroid"]),
                radius=float(h["radius"]),
                area=float(h["area"]),
                n_cells=int(h["n_cells"]),
            )
            for h in data["holes"]
        ]
        return cls(
            holes=holes,
            origin=tuple(data["origin"]),
            up=tuple(data["up"]),
            mode=data["mode"],
            resolution=float(data["resolution"]),
            source_path=data.get("source_path", ""),
        )


# --------------------------------------------------------------------------- #
# up-axis estimation
# --------------------------------------------------------------------------- #
def _load_sample(input_path: Path, max_points: int, seed: int = 0) -> np.ndarray:
    """Stream the cloud and return a random XYZ subsample (``(k, 3)``)."""
    rng = np.random.default_rng(seed)
    with BinaryPcdReader(input_path) as reader:
        total = max(reader.point_count, 1)
        keep = min(1.0, max_points / total)
        parts: list[np.ndarray] = []
        for chunk in reader.chunks():
            if keep < 1.0:
                chunk = chunk[rng.random(len(chunk)) < keep]
            if len(chunk):
                parts.append(np.ascontiguousarray(chunk[:, :3], dtype=np.float64))
    return np.concatenate(parts) if parts else np.empty((0, 3), dtype=np.float64)


def estimate_up(xyz: np.ndarray, k: int = 30) -> tuple[float, float, float]:
    """Estimate the dominant surface-normal direction (approx. "up").

    Per-point normals are the smallest-eigenvector of each local neighbourhood;
    the dominant orientation is the top eigenvector of the normal scatter matrix
    ``sum(n n^T)``. For a tube the floor/ceiling are the dominant surfaces, so
    this approximates the vertical. The sign is chosen to point toward ``+Z``.

    This is a heuristic for an arbitrary SLAM frame -- prefer an explicit ``--up``
    when the geometry is known.
    """
    xyz = np.asarray(xyz, dtype=np.float64)
    n = len(xyz)
    if n < k + 1:
        return (0.0, 0.0, 1.0)
    tree = cKDTree(xyz)
    # Use a capped seed set for speed; normals from a subset are plenty.
    seeds = xyz if n <= 20_000 else xyz[np.random.default_rng(0).choice(n, 20_000, replace=False)]
    _, idx = tree.query(seeds, k=k + 1, workers=-1)
    scatter = np.zeros((3, 3))
    for nbr in idx:
        pts = xyz[nbr]
        cov = np.cov((pts - pts.mean(axis=0)).T)
        w, v = np.linalg.eigh(cov)
        nrm = v[:, 0]  # smallest eigenvalue -> surface normal
        scatter += np.outer(nrm, nrm)
    w, v = np.linalg.eigh(scatter)
    up = v[:, -1]  # dominant normal direction
    if up[2] < 0:
        up = -up
    return tuple(normalize(up))


# --------------------------------------------------------------------------- #
# grid accumulation (streamed)
# --------------------------------------------------------------------------- #
@dataclass
class _Grid:
    counts: np.ndarray   # (nx, ny) int
    zsum: np.ndarray     # (nx, ny) float, sum of z (working frame)
    zmax: np.ndarray     # (nx, ny) float, max z (-inf where empty)
    xmin: float
    ymin: float
    res: float

    @property
    def meanz(self) -> np.ndarray:
        out = np.full(self.counts.shape, np.nan)
        occ = self.counts > 0
        out[occ] = self.zsum[occ] / self.counts[occ]
        return out


def _accumulate_grid(input_path: Path, R: np.ndarray, res: float) -> _Grid:
    """Two streamed passes (extents, then grids) over the up-rotated cloud."""
    Rt = R.T
    # Pass 1: extents in the working (up -> +Z) frame.
    xmin = ymin = np.inf
    xmax = ymax = -np.inf
    with BinaryPcdReader(input_path) as reader:
        for chunk in reader.chunks():
            w = chunk[:, :3].astype(np.float64) @ Rt  # rotate: (R @ x^T)^T = x @ R^T
            xmin, xmax = min(xmin, w[:, 0].min()), max(xmax, w[:, 0].max())
            ymin, ymax = min(ymin, w[:, 1].min()), max(ymax, w[:, 1].max())
    if not np.isfinite(xmin):
        raise ValueError(f"{input_path.name}: no points to grid")

    nx = max(1, int(np.ceil((xmax - xmin) / res)) + 1)
    ny = max(1, int(np.ceil((ymax - ymin) / res)) + 1)
    if nx * ny > _MAX_GRID_CELLS:
        raise ValueError(
            f"grid would be {nx}x{ny} cells; increase --res (resolution) "
            f"(extent {xmax - xmin:.0f}x{ymax - ymin:.0f} at res {res})"
        )

    counts = np.zeros((nx, ny), dtype=np.int64)
    zsum = np.zeros((nx, ny), dtype=np.float64)
    zmax = np.full((nx, ny), -np.inf, dtype=np.float64)
    # Pass 2: accumulate.
    with BinaryPcdReader(input_path) as reader:
        for chunk in reader.chunks():
            w = chunk[:, :3].astype(np.float64) @ Rt
            ix = np.clip(((w[:, 0] - xmin) / res).astype(np.int64), 0, nx - 1)
            iy = np.clip(((w[:, 1] - ymin) / res).astype(np.int64), 0, ny - 1)
            np.add.at(counts, (ix, iy), 1)
            np.add.at(zsum, (ix, iy), w[:, 2])
            np.maximum.at(zmax, (ix, iy), w[:, 2])
    return _Grid(counts, zsum, zmax, float(xmin), float(ymin), float(res))


# --------------------------------------------------------------------------- #
# detection
# --------------------------------------------------------------------------- #
def _hole_cells_aerial(grid: _Grid) -> np.ndarray:
    """Enclosed empty cells of the occupancy grid (the void skylights)."""
    occ = grid.counts > 0
    filled = ndi.binary_fill_holes(occ)
    return filled & ~occ


def _hole_cells_ceiling(grid: _Grid, jump: float, win: int = 7) -> np.ndarray:
    """Enclosed cells where the ceiling is missing or deviates from its neighbours."""
    occ = grid.counts > 0
    H = grid.zmax.copy()
    # Fill empty cells with the global ceiling median so the local filter is stable.
    H_filled = np.where(occ, H, np.median(H[occ]) if occ.any() else 0.0)
    local = ndi.median_filter(H_filled, size=win, mode="nearest")
    anomaly = occ & (np.abs(H - local) > jump)
    empty_enclosed = ndi.binary_fill_holes(occ) & ~occ
    return anomaly | empty_enclosed


def _components_to_holes(
    hole_cells: np.ndarray, grid: _Grid, R: np.ndarray, min_area: float, max_area: float | None
) -> list[Hole]:
    """Label connected hole cells and turn each into a :class:`Hole` (local frame)."""
    labels, n = ndi.label(hole_cells)
    occ = grid.counts > 0
    meanz = grid.meanz
    cell_area = grid.res * grid.res
    holes: list[Hole] = []
    for lab in range(1, n + 1):
        cells = labels == lab
        n_cells = int(cells.sum())
        area = n_cells * cell_area
        if area < min_area or (max_area is not None and area > max_area):
            continue
        coords = np.argwhere(cells)  # (k, 2) of (ix, iy)
        cix, ciy = coords[:, 0].mean(), coords[:, 1].mean()
        wx = grid.xmin + (cix + 0.5) * grid.res
        wy = grid.ymin + (ciy + 0.5) * grid.res
        # Surrounding surface level: median z of occupied cells bordering the hole.
        border = ndi.binary_dilation(cells) & occ
        zvals = meanz[border]
        zvals = zvals[np.isfinite(zvals)]
        wz = float(np.median(zvals)) if zvals.size else float(np.nanmedian(meanz))
        local = R.T @ np.array([wx, wy, wz])  # un-rotate back to the cloud frame
        radius = float(np.sqrt(area / np.pi))
        holes.append(Hole(
            id=len(holes),
            centroid=(float(local[0]), float(local[1]), float(local[2])),
            radius=radius,
            area=float(area),
            n_cells=n_cells,
        ))
    return holes


def detect_holes(
    input_path: str | Path,
    mode: str = "aerial",
    up: tuple[float, float, float] | None = None,
    resolution: float = DEFAULT_RESOLUTION,
    min_area: float = DEFAULT_MIN_AREA,
    max_area: float | None = None,
    ceiling_jump: float = DEFAULT_CEILING_JUMP,
) -> HoleSet:
    """Detect skylight holes in ``input_path`` and return a :class:`HoleSet`.

    ``mode`` is ``"aerial"`` (enclosed voids in a top-down occupancy grid, up =
    ``+Z``) or ``"ceiling"`` (the tube ceiling, viewed along ``up``). For
    ``"ceiling"`` the ``up`` axis is used if given, else estimated. Centroids are
    returned in the cloud's **local** frame, alongside the ``up`` and ``origin``
    so :mod:`lava_pcd.register` can reconcile the two clouds.
    """
    input_path = Path(input_path)
    if input_path.suffix.lower() != ".pcd":
        raise ValueError("hole detection works on a .pcd input")
    if mode not in ("aerial", "ceiling"):
        raise ValueError(f"mode must be 'aerial' or 'ceiling', got {mode!r}")

    with BinaryPcdReader(input_path) as reader:
        origin = reader.origin

    if up is not None:
        up_vec = tuple(normalize(np.asarray(up, dtype=np.float64)))
    elif mode == "ceiling":
        up_vec = estimate_up(_load_sample(input_path, _UP_SAMPLE))
    else:
        up_vec = (0.0, 0.0, 1.0)

    R = rotation_align(np.asarray(up_vec), np.array([0.0, 0.0, 1.0]))
    grid = _accumulate_grid(input_path, R, resolution)

    if mode == "aerial":
        hole_cells = _hole_cells_aerial(grid)
    else:
        hole_cells = _hole_cells_ceiling(grid, ceiling_jump)

    holes = _components_to_holes(hole_cells, grid, R, min_area, max_area)
    # Renumber sequentially after filtering.
    for i, h in enumerate(holes):
        h.id = i
    return HoleSet(
        holes=holes,
        origin=origin,
        up=up_vec,
        mode=mode,
        resolution=resolution,
        source_path=str(input_path),
    )


# --------------------------------------------------------------------------- #
# interactive review / manual placement
# --------------------------------------------------------------------------- #
def _rotated_sample(input_path: Path, up_vec: tuple[float, float, float], max_points: int):
    sample = _load_sample(input_path, max_points)
    R = rotation_align(np.asarray(up_vec), np.array([0.0, 0.0, 1.0]))
    return sample @ R.T, R


def review_holes(
    holeset: HoleSet, input_path: str | Path, max_display_points: int = 500_000
) -> HoleSet:
    """Show detected holes over a top-down scatter; click holes to toggle/keep them.

    Holes are drawn as numbered circles in the up-rotated (top-down) view. Left
    click toggles the nearest hole on/off; close the window to accept. Returns a
    new :class:`HoleSet` with only the kept holes.
    """
    import matplotlib.pyplot as plt

    input_path = Path(input_path)
    work, R = _rotated_sample(input_path, holeset.up, max_display_points)
    cent_local = holeset.centroids()
    cent_work = cent_local @ R.T if len(cent_local) else cent_local
    keep = [True] * len(holeset.holes)

    fig, ax = plt.subplots(figsize=(10, 8))
    if len(work):
        ax.scatter(work[:, 0], work[:, 1], c=work[:, 2], s=0.5, marker=".", linewidths=0,
                   cmap="viridis")
    ax.set_aspect("equal")
    ax.set_title(f"{input_path.name} — click a circle to keep/drop it, then close")
    ax.set_xlabel("up-plane a"); ax.set_ylabel("up-plane b")

    artists = []
    for i, (cx, cy) in enumerate(cent_work[:, :2] if len(cent_work) else []):
        circ = plt.Circle((cx, cy), holeset.holes[i].radius, fill=False,
                          color="red", lw=2)
        ax.add_patch(circ)
        txt = ax.text(cx, cy, str(i), color="red", ha="center", va="center")
        artists.append((circ, txt))

    def on_click(event):
        if event.inaxes != ax or event.xdata is None or not len(cent_work):
            return
        d = np.hypot(cent_work[:, 0] - event.xdata, cent_work[:, 1] - event.ydata)
        i = int(np.argmin(d))
        keep[i] = not keep[i]
        artists[i][0].set_color("green" if not keep[i] else "red")
        artists[i][0].set_linestyle(":" if not keep[i] else "-")
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("button_press_event", on_click)
    plt.show()

    kept = [h for h, k in zip(holeset.holes, keep) if k]
    for i, h in enumerate(kept):
        h.id = i
    return HoleSet(kept, holeset.origin, holeset.up, holeset.mode,
                   holeset.resolution, holeset.source_path)


def pick_holes(
    input_path: str | Path,
    up: tuple[float, float, float] | None = None,
    resolution: float = DEFAULT_RESOLUTION,
    radius: float = 5.0,
    max_display_points: int = 500_000,
) -> HoleSet:
    """Manually place skylights by clicking a top-down view (the fallback path).

    Each left click drops a skylight at the cursor; its Z is taken from the local
    surface level. Useful when automatic detection misses or invents holes. Close
    the window to finish.
    """
    import matplotlib.pyplot as plt

    input_path = Path(input_path)
    with BinaryPcdReader(input_path) as reader:
        origin = reader.origin
    up_vec = (tuple(normalize(np.asarray(up, dtype=np.float64))) if up is not None
              else estimate_up(_load_sample(input_path, _UP_SAMPLE)))
    work, R = _rotated_sample(input_path, up_vec, max_display_points)
    tree = cKDTree(work[:, :2]) if len(work) else None

    picks: list[tuple[float, float, float]] = []
    fig, ax = plt.subplots(figsize=(10, 8))
    if len(work):
        ax.scatter(work[:, 0], work[:, 1], c=work[:, 2], s=0.5, marker=".", linewidths=0,
                   cmap="viridis")
    ax.set_aspect("equal")
    ax.set_title(f"{input_path.name} — click each skylight centre, then close")

    def on_click(event):
        if event.inaxes != ax or event.xdata is None or tree is None:
            return
        # Local surface level: median z of nearby points.
        nbr = tree.query_ball_point([event.xdata, event.ydata], radius)
        wz = float(np.median(work[nbr, 2])) if nbr else 0.0
        picks.append((event.xdata, event.ydata, wz))
        ax.plot(event.xdata, event.ydata, "r+", markersize=12, mew=2)
        ax.add_patch(plt.Circle((event.xdata, event.ydata), radius, fill=False, color="red"))
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("button_press_event", on_click)
    plt.show()

    holes: list[Hole] = []
    for i, (wx, wy, wz) in enumerate(picks):
        local = R.T @ np.array([wx, wy, wz])
        holes.append(Hole(
            id=i,
            centroid=(float(local[0]), float(local[1]), float(local[2])),
            radius=float(radius),
            area=float(np.pi * radius * radius),
            n_cells=0,
        ))
    return HoleSet(holes, origin, up_vec, "manual", resolution, str(input_path))
