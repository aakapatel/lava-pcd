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
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from lava_pcd.geometry import normalize, rotation_align
from lava_pcd.io.pcd_reader import BinaryPcdReader

DEFAULT_RESOLUTION = 1.0          # grid cell size in coordinate units (metres)
DEFAULT_MIN_AREA = 4.0            # ignore voids smaller than this (m^2)
DEFAULT_MIN_DENSITY = 1.0         # cells with fewer points than this count as empty
DEFAULT_SMOOTH = 0.0             # Gaussian sigma (cells) applied to counts before threshold
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
class Occupancy:
    """A 2-D occupancy histogram of a cloud, taken looking along its up-axis.

    ``counts[ix, iy]`` is the number of points whose up-rotated XY falls in cell
    ``(ix, iy)``; this is the image you inspect to see skylights (low/zero count
    regions). ``zmax`` is the per-cell ceiling height (used by ``ceiling`` mode).
    The cloud was rotated by ``R`` (up -> +Z) before binning, so ``R`` maps the
    cloud's local frame into this working frame and ``R.T`` maps back.
    """

    counts: np.ndarray   # (nx, ny) int
    zsum: np.ndarray     # (nx, ny) float, sum of z (working frame)
    zmax: np.ndarray     # (nx, ny) float, max z (-inf where empty)
    xmin: float
    ymin: float
    res: float
    R: np.ndarray = field(default_factory=lambda: np.eye(3))
    up: tuple[float, float, float] = (0.0, 0.0, 1.0)
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    mode: str = "aerial"

    @property
    def meanz(self) -> np.ndarray:
        out = np.full(self.counts.shape, np.nan)
        occ = self.counts > 0
        out[occ] = self.zsum[occ] / self.counts[occ]
        return out

    @property
    def extent(self) -> list[float]:
        """``[x0, x1, y0, y1]`` in working-frame coords (for ``imshow``)."""
        nx, ny = self.counts.shape
        return [self.xmin, self.xmin + nx * self.res,
                self.ymin, self.ymin + ny * self.res]


# Backwards-compatible alias (older name used internally).
_Grid = Occupancy


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
def _occupied_mask(
    grid: Occupancy, min_density: float, smooth: float
) -> np.ndarray:
    """Boolean "occupied ground" mask from the density image.

    A cell counts as occupied when it holds at least ``min_density`` points
    (optionally after Gaussian-smoothing the count image by ``smooth`` cells, so
    a few stray returns inside a hole don't keep it "occupied"). Raising
    ``min_density`` turns *less-occupied* skylights into voids, not just empty
    ones.
    """
    counts = grid.counts.astype(np.float64)
    if smooth > 0:
        counts = ndi.gaussian_filter(counts, smooth)
    return counts >= min_density


def _hole_cells_aerial(
    grid: Occupancy, min_density: float, smooth: float
) -> tuple[np.ndarray, np.ndarray]:
    """Enclosed low-density cells of the occupancy grid (the void skylights).

    Returns ``(hole_cells, occ)`` so callers can visualise the occupied mask too.
    """
    occ = _occupied_mask(grid, min_density, smooth)
    filled = ndi.binary_fill_holes(occ)
    return filled & ~occ, occ


def _hole_cells_ceiling(
    grid: Occupancy, jump: float, min_density: float, smooth: float, win: int = 7
) -> tuple[np.ndarray, np.ndarray]:
    """Enclosed cells where the ceiling is missing or deviates from its neighbours."""
    occ = _occupied_mask(grid, min_density, smooth)
    H = grid.zmax.copy()
    # Fill empty cells with the global ceiling median so the local filter is stable.
    H_filled = np.where(occ, H, np.median(H[occ]) if occ.any() else 0.0)
    local = ndi.median_filter(H_filled, size=win, mode="nearest")
    anomaly = occ & (np.abs(H - local) > jump)
    empty_enclosed = ndi.binary_fill_holes(occ) & ~occ
    return (anomaly | empty_enclosed), occ


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


def _resolve_up(
    input_path: Path, mode: str, up: tuple[float, float, float] | None
) -> tuple[float, float, float]:
    if up is not None:
        return tuple(normalize(np.asarray(up, dtype=np.float64)))
    if mode == "ceiling":
        return estimate_up(_load_sample(input_path, _UP_SAMPLE))
    return (0.0, 0.0, 1.0)


def build_occupancy(
    input_path: str | Path,
    mode: str = "aerial",
    up: tuple[float, float, float] | None = None,
    resolution: float = DEFAULT_RESOLUTION,
) -> Occupancy:
    """Project ``input_path`` to a 2-D occupancy histogram (looking along up).

    This is the image to inspect when tuning detection: skylights show up as
    low/zero-count regions. ``mode`` ``"aerial"`` bins top-down (up = ``+Z``);
    ``"ceiling"`` rotates by ``up`` (given or estimated) first. Returns an
    :class:`Occupancy` carrying the count image plus the frame it was taken in.
    """
    input_path = Path(input_path)
    if input_path.suffix.lower() != ".pcd":
        raise ValueError("hole detection works on a .pcd input")
    if mode not in ("aerial", "ceiling"):
        raise ValueError(f"mode must be 'aerial' or 'ceiling', got {mode!r}")
    with BinaryPcdReader(input_path) as reader:
        origin = reader.origin
    up_vec = _resolve_up(input_path, mode, up)
    R = rotation_align(np.asarray(up_vec), np.array([0.0, 0.0, 1.0]))
    grid = _accumulate_grid(input_path, R, resolution)
    grid.R = R
    grid.up = up_vec
    grid.origin = origin
    grid.mode = mode
    return grid


def detect_skylights(
    grid: Occupancy,
    min_area: float = DEFAULT_MIN_AREA,
    max_area: float | None = None,
    min_density: float = DEFAULT_MIN_DENSITY,
    smooth: float = DEFAULT_SMOOTH,
    ceiling_jump: float = DEFAULT_CEILING_JUMP,
) -> tuple[list[Hole], np.ndarray]:
    """Find skylights in a prebuilt :class:`Occupancy`. Returns ``(holes, mask)``.

    ``mask`` is the boolean grid of detected hole cells (handy to overlay on the
    occupancy image). See :func:`detect_holes` for the parameter meanings.
    """
    if grid.mode == "aerial":
        hole_cells, _ = _hole_cells_aerial(grid, min_density, smooth)
    else:
        hole_cells, _ = _hole_cells_ceiling(grid, ceiling_jump, min_density, smooth)
    holes = _components_to_holes(hole_cells, grid, grid.R, min_area, max_area)
    for i, h in enumerate(holes):  # renumber after area filtering
        h.id = i
    return holes, hole_cells


def detect_holes(
    input_path: str | Path,
    mode: str = "aerial",
    up: tuple[float, float, float] | None = None,
    resolution: float = DEFAULT_RESOLUTION,
    min_area: float = DEFAULT_MIN_AREA,
    max_area: float | None = None,
    min_density: float = DEFAULT_MIN_DENSITY,
    smooth: float = DEFAULT_SMOOTH,
    ceiling_jump: float = DEFAULT_CEILING_JUMP,
) -> HoleSet:
    """Detect skylight holes in ``input_path`` and return a :class:`HoleSet`.

    ``mode`` is ``"aerial"`` (enclosed low-density regions in a top-down
    occupancy grid, up = ``+Z``) or ``"ceiling"`` (the tube ceiling, viewed along
    ``up`` -- given or estimated). A cell counts as ground when it holds at least
    ``min_density`` points (optionally after Gaussian smoothing the count image by
    ``smooth`` cells), so *less-occupied* skylights are caught, not only empty
    ones. ``min_area`` / ``max_area`` filter candidates by footprint. Centroids
    are returned in the cloud's **local** frame, alongside ``up`` and ``origin``.
    """
    grid = build_occupancy(input_path, mode=mode, up=up, resolution=resolution)
    holes, _ = detect_skylights(
        grid, min_area=min_area, max_area=max_area, min_density=min_density,
        smooth=smooth, ceiling_jump=ceiling_jump,
    )
    return HoleSet(
        holes=holes,
        origin=grid.origin,
        up=grid.up,
        mode=mode,
        resolution=resolution,
        source_path=str(Path(input_path)),
    )


# --------------------------------------------------------------------------- #
# interactive review / manual placement
# --------------------------------------------------------------------------- #
def _rotated_sample(input_path: Path, up_vec: tuple[float, float, float], max_points: int):
    sample = _load_sample(input_path, max_points)
    R = rotation_align(np.asarray(up_vec), np.array([0.0, 0.0, 1.0]))
    return sample @ R.T, R


def show_occupancy(
    grid: Occupancy,
    holeset: HoleSet | None = None,
    hole_cells: np.ndarray | None = None,
    interactive: bool = False,
    log: bool = True,
    title: str = "",
) -> HoleSet | None:
    """Render the 2-D occupancy histogram, optionally with detected skylights.

    The count image is shown with ``imshow`` (log scale by default, so sparse
    returns inside holes stay visible); detected hole cells are tinted red and
    each :class:`Hole` is drawn as a numbered circle. With ``interactive=True``
    and a ``holeset``, left-clicking a circle toggles that hole on/off and the
    kept :class:`HoleSet` is returned when the window is closed. Otherwise this
    is purely a viewer and returns ``holeset`` unchanged (or ``None``).
    """
    import matplotlib.pyplot as plt

    nx, ny = grid.counts.shape
    counts = grid.counts.astype(np.float64)
    img = (np.log1p(counts) if log else counts).T  # transpose -> (row=y, col=x)

    fig, ax = plt.subplots(figsize=(11, 8))
    im = ax.imshow(img, origin="lower", extent=grid.extent, cmap="viridis",
                   aspect="equal", interpolation="nearest")
    fig.colorbar(im, ax=ax,
                 label=("log(1 + points)" if log else "points") + " per cell")
    ax.set_xlabel("up-plane a"); ax.set_ylabel("up-plane b")
    ax.set_title(title or "occupancy histogram "
                 f"({grid.mode}, res {grid.res:g})")

    if hole_cells is not None and hole_cells.any():
        rgba = np.zeros((ny, nx, 4))
        rgba[..., 0] = 1.0                       # red
        rgba[..., 3] = hole_cells.T.astype(float) * 0.45
        ax.imshow(rgba, origin="lower", extent=grid.extent, aspect="equal",
                  interpolation="nearest")

    keep: list[bool] = []
    artists = []
    cent_work = np.empty((0, 2))
    if holeset is not None and holeset.holes:
        cent_work = (holeset.centroids() @ grid.R.T)[:, :2]
        keep = [True] * len(holeset.holes)
        for i, (cx, cy) in enumerate(cent_work):
            circ = plt.Circle((cx, cy), max(holeset.holes[i].radius, grid.res),
                              fill=False, color="red", lw=2)
            ax.add_patch(circ)
            txt = ax.text(cx, cy, str(i), color="white", ha="center", va="center",
                          fontsize=8)
            artists.append((circ, txt))

    if interactive and len(cent_work):
        ax.set_title(title or "click a circle to keep/drop it, then close")

        def on_click(event):
            if event.inaxes != ax or event.xdata is None:
                return
            d = np.hypot(cent_work[:, 0] - event.xdata, cent_work[:, 1] - event.ydata)
            i = int(np.argmin(d))
            keep[i] = not keep[i]
            artists[i][0].set_color("lime" if not keep[i] else "red")
            artists[i][0].set_linestyle(":" if not keep[i] else "-")
            fig.canvas.draw_idle()

        fig.canvas.mpl_connect("button_press_event", on_click)

    plt.show()

    if holeset is None:
        return None
    if not interactive:
        return holeset
    kept = [h for h, k in zip(holeset.holes, keep) if k]
    for i, h in enumerate(kept):
        h.id = i
    return HoleSet(kept, holeset.origin, holeset.up, holeset.mode,
                   holeset.resolution, holeset.source_path)


def review_holes(
    grid: Occupancy, holeset: HoleSet, hole_cells: np.ndarray | None = None
) -> HoleSet:
    """Interactive review over the occupancy image (see :func:`show_occupancy`)."""
    return show_occupancy(grid, holeset, hole_cells, interactive=True)


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
