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
DEFAULT_RELATIVE_WINDOW = 21      # window (cells) for the local reference density
DEFAULT_EDGE_MARGIN = 0.0         # reject holes within this distance of the boundary
DEFAULT_MERGE_FACTOR = 1.0        # scale on ellipse radii for the overlap-merge test
DEFAULT_CEILING_JUMP = 1.0        # ceiling-height deviation flagged as "open" (m)
_MAX_GRID_CELLS = 50_000_000      # guard against an absurdly fine grid
_UP_SAMPLE = 200_000             # points sampled for up-axis estimation


@dataclass
class Hole:
    """A single detected skylight, in the cloud's local frame.

    Shape is described both as an equivalent-circle ``radius`` and as an
    equivalent **ellipse** (same second moments): ``semi_major``/``semi_minor``
    axes and ``orientation`` -- the major-axis angle in the up-plane, in radians
    in ``[0, pi)`` (0 = the up-plane's first axis). The ellipse captures the
    skylight's rough elongation and bearing.
    """

    id: int
    centroid: tuple[float, float, float]
    radius: float
    area: float
    n_cells: int
    semi_major: float = 0.0
    semi_minor: float = 0.0
    orientation: float = 0.0


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
                semi_major=float(h.get("semi_major", h["radius"])),
                semi_minor=float(h.get("semi_minor", h["radius"])),
                orientation=float(h.get("orientation", 0.0)),
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
    grid: Occupancy,
    min_density: float,
    smooth: float,
    relative: float | None = None,
    relative_window: int = DEFAULT_RELATIVE_WINDOW,
) -> np.ndarray:
    """Boolean "occupied ground" mask from the density image.

    Two thresholding modes:

    * **absolute** (default) -- a cell is occupied when it holds at least
      ``min_density`` points. Simple, but assumes roughly uniform ground density.
    * **relative** (``relative`` set) -- a cell is occupied when its count is at
      least ``relative`` times the *local* reference density (a median over a
      ``relative_window``-cell window). This adapts to ground density that varies
      across the map, so a hole in a sparse area is still caught and dense areas
      aren't over-flagged. ``min_density`` is ignored in this mode.

    In both modes the count image may first be Gaussian-smoothed by ``smooth``
    cells so a few stray returns inside a hole don't keep it "occupied".
    """
    counts = grid.counts.astype(np.float64)
    if smooth > 0:
        counts = ndi.gaussian_filter(counts, smooth)
    if relative is not None:
        ref = ndi.median_filter(counts, size=int(relative_window), mode="nearest")
        return counts >= relative * ref
    return counts >= min_density


def _hole_cells_aerial(
    grid: Occupancy, min_density: float, smooth: float,
    relative: float | None = None, relative_window: int = DEFAULT_RELATIVE_WINDOW,
) -> tuple[np.ndarray, np.ndarray]:
    """Enclosed low-density cells of the occupancy grid (the void skylights).

    Returns ``(hole_cells, occ)`` so callers can visualise the occupied mask too.
    """
    occ = _occupied_mask(grid, min_density, smooth, relative, relative_window)
    filled = ndi.binary_fill_holes(occ)
    return filled & ~occ, occ


def _hole_cells_ceiling(
    grid: Occupancy, jump: float, min_density: float, smooth: float,
    relative: float | None = None, relative_window: int = DEFAULT_RELATIVE_WINDOW,
    win: int = 7,
) -> tuple[np.ndarray, np.ndarray]:
    """Enclosed cells where the ceiling is missing or deviates from its neighbours."""
    occ = _occupied_mask(grid, min_density, smooth, relative, relative_window)
    H = grid.zmax.copy()
    # Fill empty cells with the global ceiling median so the local filter is stable.
    H_filled = np.where(occ, H, np.median(H[occ]) if occ.any() else 0.0)
    local = ndi.median_filter(H_filled, size=win, mode="nearest")
    anomaly = occ & (np.abs(H - local) > jump)
    empty_enclosed = ndi.binary_fill_holes(occ) & ~occ
    return (anomaly | empty_enclosed), occ


def _fit_ellipse(
    coords: np.ndarray, grid: _Grid
) -> tuple[float, float, float]:
    """Equivalent ellipse (same 2nd moments) of a set of cell indices.

    ``coords`` is ``(k, 2)`` of ``(ix, iy)``. Returns ``(semi_major, semi_minor,
    orientation)`` in coordinate units / radians, computed in the up-plane
    (working) frame. For a uniform filled ellipse with semi-axes ``a, b`` the
    coordinate variance along each axis is ``a^2/4``, so ``semi = 2*sqrt(eig)``.
    Returns zeros for a degenerate (too few cells) region.
    """
    if len(coords) < 3:
        return 0.0, 0.0, 0.0
    wx = grid.xmin + (coords[:, 0] + 0.5) * grid.res
    wy = grid.ymin + (coords[:, 1] + 0.5) * grid.res
    cov = np.cov(np.vstack([wx, wy]))
    w, v = np.linalg.eigh(cov)  # ascending eigenvalues
    semi_minor = 2.0 * np.sqrt(max(w[0], 0.0))
    semi_major = 2.0 * np.sqrt(max(w[1], 0.0))
    major_vec = v[:, 1]
    orientation = float(np.arctan2(major_vec[1], major_vec[0]) % np.pi)
    return float(semi_major), float(semi_minor), orientation


def _ellipse_radius(a: float, b: float, dphi: float) -> float:
    """Radius of an ellipse (semi-axes ``a``, ``b``) in direction ``dphi`` from its
    major axis."""
    denom = np.hypot(b * np.cos(dphi), a * np.sin(dphi))
    return float(a * b / denom) if denom > 1e-9 else float(max(a, b))


def _ellipses_overlap(p1: dict, p2: dict, factor: float) -> bool:
    """Approximate ellipse-ellipse overlap test.

    Each ``p`` is ``{cx, cy, a, b, th}``. The ellipses are deemed to overlap when
    the centre distance is within ``factor`` times the sum of each ellipse's radius
    measured along the line joining the centres (exact for circles; a good, cheap
    approximation for ellipses).
    """
    dx, dy = p2["cx"] - p1["cx"], p2["cy"] - p1["cy"]
    d = float(np.hypot(dx, dy))
    if d < 1e-9:
        return True
    phi = np.arctan2(dy, dx)
    r1 = _ellipse_radius(p1["a"], p1["b"], phi - p1["th"])
    r2 = _ellipse_radius(p2["a"], p2["b"], phi - p2["th"])
    return d <= factor * (r1 + r2)


def _merge_components(comps: list[dict], grid: _Grid, factor: float) -> list[dict]:
    """Union-find merge of components whose fitted ellipses overlap.

    Overlap is transitive, so a chain of overlapping ellipses collapses to one.
    Merged components pool their cells and the ellipse/centroid/level are recomputed
    from the union.
    """
    m = len(comps)
    parent = list(range(m))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(m):
        for j in range(i + 1, m):
            if find(i) != find(j) and _ellipses_overlap(comps[i], comps[j], factor):
                parent[find(i)] = find(j)

    groups: dict[int, list[int]] = {}
    for i in range(m):
        groups.setdefault(find(i), []).append(i)

    merged: list[dict] = []
    for members in groups.values():
        if len(members) == 1:
            merged.append(comps[members[0]])
            continue
        coords = np.vstack([comps[k]["coords"] for k in members])
        n_cells = int(sum(comps[k]["n_cells"] for k in members))
        wz = sum(comps[k]["wz"] * comps[k]["n_cells"] for k in members) / n_cells
        merged.append(_component(coords, n_cells, float(wz), grid))
    return merged


def _component(coords: np.ndarray, n_cells: int, wz: float, grid: _Grid) -> dict:
    """Bundle a hole component's geometry (working frame) for filtering/merging."""
    cix, ciy = coords[:, 0].mean(), coords[:, 1].mean()
    a, b, th = _fit_ellipse(coords, grid)
    radius = float(np.sqrt(n_cells * grid.res * grid.res / np.pi))
    if a == 0.0:  # too few cells to fit -> circle
        a = b = radius
    return {
        "coords": coords, "n_cells": n_cells, "wz": wz, "radius": radius,
        "cx": grid.xmin + (cix + 0.5) * grid.res,
        "cy": grid.ymin + (ciy + 0.5) * grid.res,
        "a": a, "b": b, "th": th,
    }


def _components_to_holes(
    hole_cells: np.ndarray, grid: _Grid, R: np.ndarray, min_area: float,
    max_area: float | None, interior: np.ndarray | None = None,
    merge_overlap: bool = True, merge_factor: float = 1.0,
) -> tuple[list[Hole], np.ndarray]:
    """Label connected hole cells and turn each into a :class:`Hole` (local frame).

    Returns ``(holes, kept_mask)`` where ``kept_mask`` marks only the cells of the
    components that pass the filters -- i.e. the mask that corresponds to the
    returned holes (handy to overlay as a faithful preview). When ``interior`` is
    given, a component is kept only if its centroid lies inside it (used to reject
    holes hugging the cloud boundary). When ``merge_overlap`` is set, kept
    components whose fitted ellipses overlap (scaled by ``merge_factor``) are fused
    into one hole.
    """
    labels, n = ndi.label(hole_cells)
    occ = grid.counts > 0
    meanz = grid.meanz
    cell_area = grid.res * grid.res

    # Pass 1: collect components that pass the area + interior filters.
    comps: list[dict] = []
    for lab in range(1, n + 1):
        cells = labels == lab
        n_cells = int(cells.sum())
        area = n_cells * cell_area
        if area < min_area or (max_area is not None and area > max_area):
            continue
        coords = np.argwhere(cells)  # (k, 2) of (ix, iy)
        cix, ciy = coords[:, 0].mean(), coords[:, 1].mean()
        # Drop holes whose centre sits within the boundary margin (edge nicks).
        if interior is not None:
            ci, cj = int(round(cix)), int(round(ciy))
            if not (0 <= ci < interior.shape[0] and 0 <= cj < interior.shape[1]
                    and interior[ci, cj]):
                continue
        # Surrounding surface level: median z of occupied cells bordering the hole.
        border = ndi.binary_dilation(cells) & occ
        zvals = meanz[border]
        zvals = zvals[np.isfinite(zvals)]
        wz = float(np.median(zvals)) if zvals.size else float(np.nanmedian(meanz))
        comps.append(_component(coords, n_cells, wz, grid))

    if merge_overlap and len(comps) > 1:
        comps = _merge_components(comps, grid, merge_factor)

    # Pass 2: emit holes.
    holes: list[Hole] = []
    kept_mask = np.zeros_like(hole_cells, dtype=bool)
    for c in comps:
        coords = c["coords"]
        kept_mask[coords[:, 0], coords[:, 1]] = True
        local = R.T @ np.array([c["cx"], c["cy"], c["wz"]])  # un-rotate to cloud frame
        holes.append(Hole(
            id=len(holes),
            centroid=(float(local[0]), float(local[1]), float(local[2])),
            radius=c["radius"],
            area=float(c["n_cells"] * cell_area),
            n_cells=c["n_cells"],
            semi_major=c["a"],
            semi_minor=c["b"],
            orientation=c["th"],
        ))
    return holes, kept_mask


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
    relative: float | None = None,
    relative_window: int = DEFAULT_RELATIVE_WINDOW,
    edge_margin: float = DEFAULT_EDGE_MARGIN,
    merge_overlap: bool = True,
    merge_factor: float = DEFAULT_MERGE_FACTOR,
    ceiling_jump: float = DEFAULT_CEILING_JUMP,
) -> tuple[list[Hole], np.ndarray]:
    """Find skylights in a prebuilt :class:`Occupancy`. Returns ``(holes, mask)``.

    ``mask`` is the boolean grid of the cells that became holes -- already
    area-filtered, so it matches the returned ``holes`` (handy to overlay on the
    occupancy image). See :func:`detect_holes` for the parameter meanings.
    """
    if grid.mode == "aerial":
        hole_cells, occ = _hole_cells_aerial(
            grid, min_density, smooth, relative, relative_window)
    else:
        hole_cells, occ = _hole_cells_ceiling(
            grid, ceiling_jump, min_density, smooth, relative, relative_window)
    interior = None
    if edge_margin > 0:
        # Keep only holes whose centre is >= edge_margin inside the cloud boundary:
        # erode the solid (filled) footprint inward by that many cells.
        iters = max(1, int(round(edge_margin / grid.res)))
        interior = ndi.binary_erosion(ndi.binary_fill_holes(occ), iterations=iters)
    holes, kept_mask = _components_to_holes(
        hole_cells, grid, grid.R, min_area, max_area, interior,
        merge_overlap=merge_overlap, merge_factor=merge_factor)
    for i, h in enumerate(holes):  # renumber after area filtering
        h.id = i
    return holes, kept_mask


def detect_holes(
    input_path: str | Path,
    mode: str = "aerial",
    up: tuple[float, float, float] | None = None,
    resolution: float = DEFAULT_RESOLUTION,
    min_area: float = DEFAULT_MIN_AREA,
    max_area: float | None = None,
    min_density: float = DEFAULT_MIN_DENSITY,
    smooth: float = DEFAULT_SMOOTH,
    relative: float | None = None,
    relative_window: int = DEFAULT_RELATIVE_WINDOW,
    edge_margin: float = DEFAULT_EDGE_MARGIN,
    merge_overlap: bool = True,
    merge_factor: float = DEFAULT_MERGE_FACTOR,
    ceiling_jump: float = DEFAULT_CEILING_JUMP,
) -> HoleSet:
    """Detect skylight holes in ``input_path`` and return a :class:`HoleSet`.

    ``mode`` is ``"aerial"`` (enclosed low-density regions in a top-down
    occupancy grid, up = ``+Z``) or ``"ceiling"`` (the tube ceiling, viewed along
    ``up`` -- given or estimated). A cell counts as ground when it holds at least
    ``min_density`` points; or, if ``relative`` is set, at least ``relative`` times
    the local median density (a ``relative_window``-cell window) -- which adapts to
    ground density that varies across the map. The count image may first be
    Gaussian-smoothed by ``smooth`` cells. ``min_area`` / ``max_area`` filter
    candidates by footprint, and ``edge_margin`` rejects holes within that distance
    of the cloud boundary (e.g. the ragged tube-ceiling rim). When ``merge_overlap``
    is set (default), holes whose fitted ellipses overlap (scaled by ``merge_factor``)
    are fused. Centroids are returned in the cloud's **local** frame, alongside
    ``up`` and ``origin``.
    """
    grid = build_occupancy(input_path, mode=mode, up=up, resolution=resolution)
    holes, _ = detect_skylights(
        grid, min_area=min_area, max_area=max_area, min_density=min_density,
        smooth=smooth, relative=relative, relative_window=relative_window,
        edge_margin=edge_margin, merge_overlap=merge_overlap,
        merge_factor=merge_factor, ceiling_jump=ceiling_jump,
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
    each :class:`Hole` is drawn as its fitted ellipse, numbered. With
    ``interactive=True`` and a ``holeset``, left-clicking an ellipse toggles that
    hole on/off and the kept :class:`HoleSet` is returned when the window is
    closed. Otherwise this is purely a viewer and returns ``holeset`` unchanged
    (or ``None``).
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

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
            h = holeset.holes[i]
            a = h.semi_major if h.semi_major > 0 else h.radius
            b = h.semi_minor if h.semi_minor > 0 else h.radius
            ell = Ellipse(
                (cx, cy), width=2 * max(a, grid.res), height=2 * max(b, grid.res),
                angle=np.degrees(h.orientation), fill=False, color="red", lw=2,
            )
            ax.add_patch(ell)
            txt = ax.text(cx, cy, str(i), color="white", ha="center", va="center",
                          fontsize=8)
            artists.append((ell, txt))

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
            semi_major=float(radius),
            semi_minor=float(radius),
            orientation=0.0,
        ))
    return HoleSet(holes, origin, up_vec, "manual", resolution, str(input_path))
