"""Apply a registration transform, refine it on the rims, and merge two clouds.

Given a :class:`~lava_pcd.register.Transform` (tube-local -> aerial-local), this
module:

* :func:`apply_transform` -- streams a cloud and writes it rotated+translated.
* :func:`icp_refine` -- a small point-to-point ICP run **only on the matched
  skylight rims** (global overlap is too small for global ICP), to tighten the
  landmark-based alignment.
* :func:`merge_clouds` -- transforms the tube into the aerial frame and writes a
  single combined cloud (optionally with a ``source`` channel so the two inputs
  can be told apart in a viewer).

Everything streams in chunks; the ICP works on small rim patches, so memory stays
bounded.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

from lava_pcd.geometry import (
    as_matrix,
    normalize,
    rotation_about_axis,
    transform_points,
)
from lava_pcd.io.pcd_reader import BinaryPcdReader
from lava_pcd.io.pcd_writer import BinaryPcdWriter
from lava_pcd.register import Transform

DEFAULT_RIM_RADIUS = 15.0   # horizontal radius around each skylight used for ICP
DEFAULT_RIM_HEIGHT = 6.0    # vertical band around each opening used for ICP (coord units)
DEFAULT_ICP_ITERS = 30
DEFAULT_CMAP = "viridis"    # elevation colormap for the tube in a merged cloud


@dataclass
class MergeResult:
    """Outcome of a :func:`merge_clouds` call."""

    point_count: int
    aerial_count: int
    tube_count: int
    refined: bool
    rms: float
    output_path: Path
    fields: tuple[str, ...]
    origin: tuple[float, float, float]


def apply_transform(
    input_path: str | Path,
    output_path: str | Path,
    matrix: np.ndarray,
    chunk_size: int = 5_000_000,
    show_progress: bool = True,
) -> Path:
    """Stream ``input_path``, apply the 4x4 ``matrix`` to XYZ, write ``output_path``.

    Non-coordinate fields and the local origin are preserved unchanged.
    """
    input_path, output_path = Path(input_path), Path(output_path)
    if input_path.suffix.lower() != ".pcd" or output_path.suffix.lower() != ".pcd":
        raise ValueError("transform works on .pcd input and output")
    matrix = np.asarray(matrix, dtype=np.float64)
    with BinaryPcdReader(input_path) as reader:
        if input_path.resolve() == output_path.resolve():
            raise ValueError("input and output must be different files")
        fields, origin, total = reader.fields, reader.origin, reader.point_count
        bar = tqdm(total=total, unit="pts", unit_scale=True, desc=input_path.name,
                   disable=not show_progress)
        with BinaryPcdWriter(output_path, max_points=total, fields=fields,
                             origin=origin) as writer, bar:
            for chunk in reader.chunks(chunk_size):
                out = chunk.copy()
                out[:, :3] = transform_points(chunk[:, :3].astype(np.float64), matrix)
                writer.write_chunk(np.ascontiguousarray(out, dtype=np.float32))
                bar.update(len(chunk))
    return output_path


def _plane_basis(up: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Right-handed in-plane basis ``(e1, e2)`` with ``e1 x e2 = up``."""
    up = normalize(up)
    seed = np.array([1.0, 0.0, 0.0]) if abs(up[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = normalize(seed - up * np.dot(seed, up))
    e2 = np.cross(up, e1)
    return e1, e2


def _fit_4dof(src: np.ndarray, dst: np.ndarray, up: np.ndarray) -> np.ndarray:
    """Best **yaw-about-up + translation** mapping ``src`` onto ``dst`` (4x4).

    Constrained rigid fit: the rotation is only about ``up`` (so the tube's tilt
    can't change), solved as a 2-D Procrustes in the up-plane plus a 1-D height
    offset. This is what stops ICP from laying the tube flat on the ground.
    """
    up = normalize(up)
    e1, e2 = _plane_basis(up)
    s2 = np.column_stack([src @ e1, src @ e2])
    d2 = np.column_stack([dst @ e1, dst @ e2])
    cs, cd = s2.mean(axis=0), d2.mean(axis=0)
    H = (s2 - cs).T @ (d2 - cd)
    U, _, Vt = np.linalg.svd(H)
    D = np.diag([1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    R2 = Vt.T @ D @ U.T
    ang = float(np.arctan2(R2[1, 0], R2[0, 0]))
    t2 = cd - R2 @ cs
    dz = float(np.mean(dst @ up - src @ up))
    R = rotation_about_axis(up, ang)
    t = e1 * t2[0] + e2 * t2[1] + up * dz
    return as_matrix(R, t)


def _gather_rim(
    input_path: Path, anchors: np.ndarray, up: np.ndarray, radius: float,
    height: float, matrix: np.ndarray | None = None,
) -> np.ndarray:
    """Stream a cloud and keep XYZ inside a **cylinder** around the nearest anchor:
    horizontal distance ``<= radius`` and vertical (along ``up``) ``<= height``.

    The height band keeps only points near each opening's level, excluding the deep
    tube body and far ground that point-to-point ICP would otherwise collapse
    together. If ``matrix`` is given, points are transformed before the test but the
    **original** coordinates are returned (so ICP can re-apply an evolving matrix).
    """
    up = normalize(up)
    tree = cKDTree(anchors)
    kept: list[np.ndarray] = []
    with BinaryPcdReader(input_path) as reader:
        for chunk in reader.chunks():
            xyz = chunk[:, :3].astype(np.float64)
            probe = transform_points(xyz, matrix) if matrix is not None else xyz
            _, j = tree.query(probe, k=1, workers=-1)
            dvec = probe - anchors[j]
            dz = dvec @ up
            horiz = np.linalg.norm(dvec - np.outer(dz, up), axis=1)
            mask = (horiz <= radius) & (np.abs(dz) <= height)
            if mask.any():
                kept.append(xyz[mask])
    return np.concatenate(kept) if kept else np.empty((0, 3), dtype=np.float64)


def icp_refine(
    aerial_path: str | Path,
    tube_path: str | Path,
    transform: Transform,
    radius: float = DEFAULT_RIM_RADIUS,
    rim_height: float = DEFAULT_RIM_HEIGHT,
    max_iters: int = DEFAULT_ICP_ITERS,
    max_shift: float | None = None,
    tol: float = 1e-4,
) -> Transform:
    """Refine ``transform`` with a **constrained** ICP on the matched skylight rims.

    Only points within a cylinder (``radius`` horizontally, ``rim_height``
    vertically) of each matched rim centre (``transform.anchors``, aerial frame) are
    used, and the fit is **4-DOF** (yaw about the aerial up-axis + translation), so
    the tube's tilt is fixed and it cannot be flattened onto the ground. Far
    correspondences are trimmed each iteration. As a safety net the result is
    **rejected** -- the landmark alignment kept, with a warning -- if it moves the
    matched skylights by more than ``max_shift`` (default: ``radius``) or makes the
    rim RMS worse.
    """
    anchors = np.asarray(transform.anchors, dtype=np.float64)
    if len(anchors) == 0:
        return transform
    up = normalize(np.asarray(transform.aerial_up, dtype=np.float64))
    M0 = transform.array
    if max_shift is None:
        max_shift = radius

    aerial_rim = _gather_rim(Path(aerial_path), anchors, up, radius, rim_height)
    tube_rim = _gather_rim(Path(tube_path), anchors, up, radius, rim_height, matrix=M0)
    if len(aerial_rim) < 4 or len(tube_rim) < 4:
        return transform  # not enough rim overlap to refine; keep the landmark fit

    tree = cKDTree(aerial_rim)

    def rim_rms(M: np.ndarray) -> float:
        d, _ = tree.query(transform_points(tube_rim, M), k=1, workers=-1)
        return float(np.sqrt(np.mean(np.minimum(d, radius) ** 2)))

    old_rms = rim_rms(M0)
    M = M0
    for _ in range(max_iters):
        dist, idx = tree.query(transform_points(tube_rim, M), k=1, workers=-1)
        # Trim: keep the closer correspondences, ignore far (non-overlapping) pairs.
        thr = min(radius, float(np.percentile(dist, 70)))
        keep = dist <= max(thr, 1e-6)
        if keep.sum() < 4:
            break
        M_new = _fit_4dof(tube_rim[keep], aerial_rim[idx[keep]], up)
        if np.linalg.norm(M_new - M) < tol:
            M = M_new
            break
        M = M_new

    # Safety: how far did the refinement move the (already-aligned) skylights?
    delta = M @ np.linalg.inv(M0)
    moved = float(np.max(np.linalg.norm(transform_points(anchors, delta) - anchors, axis=1)))
    new_rms = rim_rms(M)

    warnings = list(transform.warnings)
    if moved > max_shift or new_rms > old_rms + 1e-9:
        warnings.append(
            f"ICP refine rejected (skylights would move {moved:.2f}, rim RMS "
            f"{old_rms:.2f} -> {new_rms:.2f}); kept the landmark alignment."
        )
        return replace(transform, warnings=warnings)

    return replace(
        transform,
        matrix=[list(map(float, row)) for row in M],
        rms=new_rms,
        mode=transform.mode + "+icp",
        warnings=warnings,
    )


def merge_clouds(
    aerial_path: str | Path,
    tube_path: str | Path,
    output_path: str | Path,
    transform: Transform,
    refine: bool = False,
    color: bool = True,
    elevation_cmap: str = DEFAULT_CMAP,
    source_field: bool = False,
    rim_radius: float = DEFAULT_RIM_RADIUS,
    rim_height: float = DEFAULT_RIM_HEIGHT,
    chunk_size: int = 5_000_000,
    show_progress: bool = True,
) -> MergeResult:
    """Transform the tube into the aerial frame and write a single merged ``.pcd``.

    The output is in **aerial-local** coordinates (the aerial origin), since the
    transform maps tube-local -> aerial-local. With ``color`` (default) the merged
    cloud carries a packed ``rgb`` field: the **aerial** points keep their own RGB
    (a grey fallback if the aerial cloud has none), while the **tube** points are
    coloured by **elevation** (output-frame Z) with the ``elevation_cmap`` colormap
    -- so the photographic aerial surface and the depth-shaded tube read distinctly
    in a viewer. Set ``color=False`` for a plain ``x y z`` cloud. ``source_field``
    adds a ``source`` channel (0 = aerial, 1 = tube) either way.
    """
    aerial_path, tube_path, output_path = (
        Path(aerial_path), Path(tube_path), Path(output_path)
    )
    for p in (aerial_path, tube_path):
        if p.suffix.lower() != ".pcd":
            raise ValueError("merge works on .pcd inputs")
    if output_path.suffix.lower() != ".pcd":
        raise ValueError("merge output must be a .pcd")
    if output_path.resolve() in (aerial_path.resolve(), tube_path.resolve()):
        raise ValueError("output must differ from both inputs")

    if refine:
        transform = icp_refine(aerial_path, tube_path, transform,
                               radius=rim_radius, rim_height=rim_height)
    M = transform.array

    fields = ("x", "y", "z")
    if color:
        fields += ("rgb",)
    if source_field:
        fields += ("source",)

    with BinaryPcdReader(aerial_path) as ra, BinaryPcdReader(tube_path) as rb:
        origin = ra.origin
        a_fields = ra.fields
        a_total, b_total = ra.point_count, rb.point_count
    max_points = a_total + b_total
    a_rgb_idx = a_fields.index("rgb") if "rgb" in a_fields else None

    # Elevation colour needs the tube's output-frame Z range up front (one pass).
    cmap = _load_cmap(elevation_cmap) if color else None
    zlo, zhi = 0.0, 1.0
    if color:
        zmin, zmax = np.inf, -np.inf
        with BinaryPcdReader(tube_path) as rb:
            for chunk in rb.chunks(chunk_size):
                zt = transform_points(chunk[:, :3].astype(np.float64), M)[:, 2]
                if len(zt):
                    zmin = min(zmin, float(zt.min()))
                    zmax = max(zmax, float(zt.max()))
        if np.isfinite(zmin):
            zlo, zhi = zmin, zmax

    bar = tqdm(total=max_points, unit="pts", unit_scale=True, desc="merge",
               disable=not show_progress)
    with BinaryPcdWriter(output_path, max_points=max_points, fields=fields,
                         origin=origin) as writer, bar:
        # Aerial cloud: identity (already in the target frame), keep its RGB.
        with BinaryPcdReader(aerial_path) as ra:
            for chunk in ra.chunks(chunk_size):
                xyz = chunk[:, :3].astype(np.float64)
                rgb = None
                if color:
                    rgb = (chunk[:, a_rgb_idx] if a_rgb_idx is not None
                           else _const_rgb(len(chunk)))
                writer.write_chunk(_assemble(xyz, rgb, 0 if source_field else None))
                bar.update(len(chunk))
        # Tube cloud: transform into the aerial frame, colour by elevation.
        with BinaryPcdReader(tube_path) as rb:
            for chunk in rb.chunks(chunk_size):
                xyz = transform_points(chunk[:, :3].astype(np.float64), M)
                rgb = _elevation_rgb(xyz[:, 2], zlo, zhi, cmap) if color else None
                writer.write_chunk(_assemble(xyz, rgb, 1 if source_field else None))
                bar.update(len(chunk))
        written = writer._written

    return MergeResult(
        point_count=written,
        aerial_count=a_total,
        tube_count=b_total,
        refined=refine,
        rms=transform.rms,
        output_path=output_path,
        fields=fields,
        origin=origin,
    )


def _assemble(xyz: np.ndarray, rgb: np.ndarray | None, source: int | None) -> np.ndarray:
    """Build the PCD-layout chunk: xyz (+ packed rgb) (+ a constant source column)."""
    cols = [np.ascontiguousarray(xyz, dtype=np.float32)]
    if rgb is not None:
        cols.append(np.ascontiguousarray(rgb, dtype=np.float32).reshape(-1, 1))
    if source is not None:
        cols.append(np.full((len(xyz), 1), float(source), dtype=np.float32))
    return np.ascontiguousarray(np.hstack(cols), dtype=np.float32)


def _pack_rgb(r: np.ndarray, g: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pack 0..255 r/g/b arrays into the PCL packed-``rgb`` float32 (0x00RRGGBB)."""
    r = np.clip(r, 0, 255).astype(np.uint32)
    g = np.clip(g, 0, 255).astype(np.uint32)
    b = np.clip(b, 0, 255).astype(np.uint32)
    packed = (r << 16) | (g << 8) | b
    return np.ascontiguousarray(packed).view(np.float32)


def _const_rgb(n: int, rgb: tuple[int, int, int] = (200, 200, 200)) -> np.ndarray:
    """A constant packed-rgb column of length ``n`` (grey by default)."""
    val = (rgb[0] << 16) | (rgb[1] << 8) | rgb[2]
    return np.full(n, val, dtype=np.uint32).view(np.float32)


def _elevation_rgb(z: np.ndarray, zlo: float, zhi: float, cmap) -> np.ndarray:
    """Packed-rgb column colouring ``z`` from ``zlo``..``zhi`` through ``cmap``."""
    rng = zhi - zlo
    t = np.clip((z - zlo) / rng, 0.0, 1.0) if rng > 1e-9 else np.zeros_like(z)
    cols = np.asarray(cmap(t))[:, :3] * 255.0
    return _pack_rgb(cols[:, 0], cols[:, 1], cols[:, 2])


def _load_cmap(name: str):
    """Look up a matplotlib colormap by name (lazy import)."""
    import matplotlib
    try:
        return matplotlib.colormaps[name]
    except KeyError as err:
        raise ValueError(f"unknown colormap {name!r}") from err
