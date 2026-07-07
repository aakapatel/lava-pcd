"""Apply a registration transform, refine it on the rims, and merge two clouds.

Given a :class:`~lava_pcd.register.Transform` (tube-local -> aerial-local), this
module:

* :func:`apply_transform` -- streams a cloud and writes it rotated+translated.
* :func:`icp_refine` -- a small **GICP** (via ``small_gicp``) run **only on the
  matched skylight rims** (global overlap is too small for global ICP), to tighten
  the landmark-based alignment; optionally constrained to 4-DOF (yaw + translation)
  so the tilted tube can't be flattened onto the ground.
* :func:`merge_clouds` -- transforms the tube into the aerial frame and writes a
  single combined cloud (optionally with a ``source`` channel so the two inputs
  can be told apart in a viewer).

Everything streams in chunks; the ICP works on small rim patches, so memory stays
bounded.
"""

from __future__ import annotations

import os
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

DEFAULT_RIM_RADIUS = 15.0   # fallback horizontal radius when a hole has no ellipse
DEFAULT_RIM_INFLATE = 1.5   # factor each hole's ellipse semi-axes are grown by for the rim
DEFAULT_RIM_HEIGHT = 6.0    # vertical band around each opening used for ICP (coord units)
DEFAULT_CMAP = "viridis"    # elevation colormap for the tube in a merged cloud

# --- small_gicp rim-refinement tunables (edit here to tune `merge --refine`) ---
DEFAULT_GICP_TYPE = "GICP"       # registration_type: 'ICP' | 'PLANE_ICP' | 'GICP' | 'VGICP'
DEFAULT_GICP_MAX_CORR = 1.0      # max_correspondence_distance (m): cap on rim-point matches
DEFAULT_GICP_DOWNSAMPLE = 0.5   # downsampling_resolution (m): voxel size the rims are reduced to
DEFAULT_GICP_VOXEL = 1.0         # voxel_resolution (m): correspondence voxels, VGICP only
DEFAULT_GICP_ITERS = 30          # max_iterations for the GICP optimisation
DEFAULT_GICP_THREADS = min(8, os.cpu_count() or 1)  # num_threads for GICP


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


def _import_small_gicp():
    """Import ``small_gicp`` lazily so the rest of the package works without it."""
    try:
        import small_gicp
    except ImportError as e:  # pragma: no cover - exercised only when missing
        raise ImportError(
            "icp_refine needs the 'small_gicp' package for GICP rim refinement; "
            "install it with `pip install small_gicp`."
        ) from e
    return small_gicp


def _project_to_4dof(matrix: np.ndarray, up: np.ndarray, src_centroid: np.ndarray) -> np.ndarray:
    """Project a 6-DOF rigid transform onto **yaw-about-up + translation** (4-DOF).

    Keeps only the rotation component about ``up`` -- dropping any pitch/roll, so the
    tube's tilt is fixed and it cannot be flattened onto the ground -- then adjusts
    the translation (``t4 = t + (R - R4) @ src_centroid``) so the constrained map is
    the least-squares-closest 4-DOF transform to ``matrix`` over points near the
    rim centroid.
    """
    up = normalize(up)
    R, t = matrix[:3, :3], matrix[:3, 3]
    e1, e2 = _plane_basis(up)
    Re1 = R @ e1
    yaw = float(np.arctan2(Re1 @ e2, Re1 @ e1))
    R4 = rotation_about_axis(up, yaw)
    t4 = t + (R - R4) @ src_centroid
    return as_matrix(R4, t4)


def _gather_rim(
    input_path: Path, anchors: np.ndarray, up: np.ndarray, radius: float,
    height: float, matrix: np.ndarray | None = None,
    axes: np.ndarray | None = None, inflate: float = 1.0,
) -> np.ndarray:
    """Stream a cloud and keep XYZ near the nearest anchor's opening.

    The vertical gate is always a band ``<= height`` (along ``up``) -- it keeps only
    points near each opening's level, excluding the deep tube body and far ground
    that GICP would otherwise try to collapse together. The horizontal gate is
    **per hole**: if ``axes`` is given (``(M, 5)`` rows ``[semi_major, semi_minor,
    mx, my, mz]`` from :attr:`Transform.anchor_axes`), a point must fall inside that
    anchor's ellipse grown by ``inflate``; otherwise it falls back to a circle of
    ``radius`` (so transforms without ellipse data still work).

    If ``matrix`` is given, points are transformed before the test but the
    **original** coordinates are returned (so the caller can re-apply an evolving
    matrix).
    """
    up = normalize(up)
    tree = cKDTree(anchors)
    if axes is not None:
        semis = np.maximum(axes[:, :2] * inflate, 1e-6)        # (M, 2) inflated semi-axes
        major = axes[:, 2:5]
        major = major / np.linalg.norm(major, axis=1, keepdims=True)
        minor = np.cross(up, major)
        minor = minor / np.linalg.norm(minor, axis=1, keepdims=True)
    kept: list[np.ndarray] = []
    with BinaryPcdReader(input_path) as reader:
        for chunk in reader.chunks():
            xyz = chunk[:, :3].astype(np.float64)
            probe = transform_points(xyz, matrix) if matrix is not None else xyz
            _, j = tree.query(probe, k=1, workers=-1)
            dvec = probe - anchors[j]
            dz = dvec @ up
            if axes is None:
                horiz = np.linalg.norm(dvec - np.outer(dz, up), axis=1)
                in_plane = horiz <= radius
            else:
                plane = dvec - dz[:, None] * up                # in-plane component
                u = np.einsum("ij,ij->i", plane, major[j])
                v = np.einsum("ij,ij->i", plane, minor[j])
                in_plane = (u / semis[j, 0]) ** 2 + (v / semis[j, 1]) ** 2 <= 1.0
            mask = in_plane & (np.abs(dz) <= height)
            if mask.any():
                kept.append(xyz[mask])
    return np.concatenate(kept) if kept else np.empty((0, 3), dtype=np.float64)


def _anchor_axes(transform: Transform, n_anchors: int) -> np.ndarray | None:
    """``(M, 5)`` ellipse array from ``transform``, or ``None`` if absent/malformed."""
    if not transform.anchor_axes:
        return None
    axes = np.asarray(transform.anchor_axes, dtype=np.float64)
    return axes if axes.shape == (n_anchors, 5) else None


def icp_refine(
    aerial_path: str | Path,
    tube_path: str | Path,
    transform: Transform,
    radius: float = DEFAULT_RIM_RADIUS,
    rim_height: float = DEFAULT_RIM_HEIGHT,
    rim_inflate: float = DEFAULT_RIM_INFLATE,
    dof: int = 4,
    reg_type: str = DEFAULT_GICP_TYPE,
    max_corr_dist: float = DEFAULT_GICP_MAX_CORR,
    downsample: float = DEFAULT_GICP_DOWNSAMPLE,
    voxel_resolution: float = DEFAULT_GICP_VOXEL,
    max_iters: int = DEFAULT_GICP_ITERS,
    num_threads: int = DEFAULT_GICP_THREADS,
    max_shift: float | None = None,
) -> Transform:
    """Refine ``transform`` with **GICP** (via ``small_gicp``) on the matched rims.

    The rim of each matched skylight is gathered per hole: points within that hole's
    ellipse (:attr:`Transform.anchor_axes`, aerial frame) grown by ``rim_inflate``
    and within ``rim_height`` of its opening level. Holes without ellipse data fall
    back to a ``radius`` circle. Only these rim patches are registered -- the global
    overlap is too small for full-cloud registration; the aerial rim is the GICP
    target, the tube rim the source, and the landmark transform the initial guess.

    ``dof`` selects the result's degrees of freedom: ``4`` (default) projects GICP's
    fit onto **yaw-about-up + translation**, so the tube's tilt is fixed and it
    cannot be flattened onto the ground; ``6`` keeps GICP's full rigid transform.
    The ``small_gicp`` knobs (``reg_type``, ``max_corr_dist``, ``downsample``,
    ``voxel_resolution``, ``max_iters``, ``num_threads``) default to the
    ``DEFAULT_GICP_*`` constants at the top of this module -- tune them there.

    As a safety net the result is **rejected** -- the landmark alignment kept, with a
    warning -- if it moves the matched skylights by more than ``max_shift``
    (default: ``radius``) or makes the rim RMS worse.
    """
    if dof not in (4, 6):
        raise ValueError("dof must be 4 (yaw + translation) or 6 (full rigid)")
    anchors = np.asarray(transform.anchors, dtype=np.float64)
    if len(anchors) == 0:
        return transform
    up = normalize(np.asarray(transform.aerial_up, dtype=np.float64))
    M0 = transform.array
    if max_shift is None:
        max_shift = radius
    axes = _anchor_axes(transform, len(anchors))

    aerial_rim = _gather_rim(Path(aerial_path), anchors, up, radius, rim_height,
                             axes=axes, inflate=rim_inflate)
    tube_rim = _gather_rim(Path(tube_path), anchors, up, radius, rim_height,
                           matrix=M0, axes=axes, inflate=rim_inflate)
    if len(aerial_rim) < 4 or len(tube_rim) < 4:
        return transform  # not enough rim overlap to refine; keep the landmark fit

    tree = cKDTree(aerial_rim)

    def rim_rms(M: np.ndarray) -> float:
        d, _ = tree.query(transform_points(tube_rim, M), k=1, workers=-1)
        return float(np.sqrt(np.mean(np.minimum(d, radius) ** 2)))

    old_rms = rim_rms(M0)

    # GICP (plane-to-plane) on the rim patches. The result maps source -> target,
    # i.e. tube-local -> aerial-local, just like the landmark transform M0.
    small_gicp = _import_small_gicp()
    result = small_gicp.align(
        aerial_rim, tube_rim,
        init_T_target_source=M0,
        registration_type=reg_type,
        voxel_resolution=voxel_resolution,
        downsampling_resolution=downsample,
        max_correspondence_distance=max_corr_dist,
        num_threads=num_threads,
        max_iterations=max_iters,
    )
    M = np.asarray(result.T_target_source, dtype=np.float64)
    if dof == 4:
        M = _project_to_4dof(M, up, tube_rim.mean(axis=0))

    # Safety: how far did the refinement move the (already-aligned) skylights?
    delta = M @ np.linalg.inv(M0)
    moved = float(np.max(np.linalg.norm(transform_points(anchors, delta) - anchors, axis=1)))
    new_rms = rim_rms(M)

    warnings = list(transform.warnings)
    if not result.converged:
        warnings.append("GICP refine did not fully converge.")
    if moved > max_shift or new_rms > old_rms + 1e-9:
        warnings.append(
            f"GICP refine rejected (skylights would move {moved:.2f}, rim RMS "
            f"{old_rms:.2f} -> {new_rms:.2f}); kept the landmark alignment."
        )
        return replace(transform, warnings=warnings)

    return replace(
        transform,
        matrix=[list(map(float, row)) for row in M],
        rms=new_rms,
        mode=f"{transform.mode}+gicp{dof}",
        warnings=warnings,
    )


def visualize_rims(
    aerial_path: str | Path,
    tube_path: str | Path,
    transform: Transform,
    refined: Transform | None = None,
    radius: float = DEFAULT_RIM_RADIUS,
    rim_height: float = DEFAULT_RIM_HEIGHT,
    rim_inflate: float = DEFAULT_RIM_INFLATE,
    title: str | None = None,
) -> None:
    """Scatter the rim points the GICP refinement actually operates on.

    Gathers the same per-hole rim patches as :func:`icp_refine` (each hole's ellipse
    grown by ``rim_inflate``, or a ``radius`` circle without ellipse data, within
    ``rim_height`` of the opening) and plots them in the aerial frame: **aerial** rim
    points blue, **tube** rim points orange, the matched anchors as black crosses.
    Each column shows a top-down view (the aerial up-plane) over an elevation view
    (an in-plane axis vs the up-axis), so both the horizontal overlap and the
    vertical rim-to-rim gap -- the offset GICP closes -- are visible. With ``refined``
    given a second column redraws the *same* tube rim under the refined transform, so
    before/after sit side by side.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    anchors = np.asarray(transform.anchors, dtype=np.float64)
    if len(anchors) == 0:
        raise ValueError("transform has no matched rim anchors to visualise")
    up = normalize(np.asarray(transform.aerial_up, dtype=np.float64))
    e1, e2 = _plane_basis(up)
    M0 = transform.array
    axes = _anchor_axes(transform, len(anchors))

    aerial_rim = _gather_rim(Path(aerial_path), anchors, up, radius, rim_height,
                             axes=axes, inflate=rim_inflate)
    # Gather the tube rim once (gated by M0, returned in tube-local coords) exactly
    # as icp_refine does, so the same points can be redrawn under either matrix.
    tube_local = _gather_rim(Path(tube_path), anchors, up, radius, rim_height,
                             matrix=M0, axes=axes, inflate=rim_inflate)

    columns = [("before refine", M0)]
    if refined is not None:
        columns.append((f"after refine  (rim RMS {refined.rms:.3f})", refined.array))

    fig, axes = plt.subplots(2, len(columns), figsize=(8 * len(columns), 10),
                             squeeze=False)
    centre = anchors.mean(axis=0)  # recentre so axes read as metres from the rims

    def proj(pts: np.ndarray, vy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        d = pts - centre
        return d @ e1, d @ vy

    for col, (subtitle, M) in enumerate(columns):
        tube_rim = transform_points(tube_local, M)
        for row, (vy, ylabel) in enumerate(
            [(e2, "up-plane b (m)"), (up, "height along up (m)")]
        ):
            ax = axes[row][col]
            ax.scatter(*proj(aerial_rim, vy), s=4, c="tab:blue", alpha=0.4)
            ax.scatter(*proj(tube_rim, vy), s=4, c="tab:orange", alpha=0.4)
            ax.scatter(*proj(anchors, vy), marker="x", s=80, c="k", zorder=3)
            ax.set_aspect("equal")
            ax.set_xlabel("up-plane a (m)")
            ax.set_ylabel(ylabel)
            if row == 0:
                ax.set_title(subtitle)

    fig.suptitle(title or (
        f"merge rims: {len(anchors)} skylight(s)   radius {radius:g}   "
        f"height {rim_height:g}"
    ))
    fig.legend(handles=[
        Line2D([0], [0], marker="o", color="tab:blue", lw=0, label="aerial rim"),
        Line2D([0], [0], marker="o", color="tab:orange", lw=0, label="tube rim"),
        Line2D([0], [0], marker="x", color="k", lw=0, label="matched anchor"),
    ], loc="lower center", ncol=3)
    plt.tight_layout(rect=(0, 0.03, 1, 0.97))
    plt.show()


def merge_clouds(
    aerial_path: str | Path,
    tube_path: str | Path,
    output_path: str | Path,
    transform: Transform,
    refine: bool = False,
    color: bool = True,
    elevation_cmap: str = DEFAULT_CMAP,
    z_offset: float = 0.0,
    source_field: bool = False,
    rim_radius: float = DEFAULT_RIM_RADIUS,
    rim_height: float = DEFAULT_RIM_HEIGHT,
    rim_inflate: float = DEFAULT_RIM_INFLATE,
    dof: int = 4,
    show_rims: bool = False,
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

    With ``refine`` the alignment is tightened by GICP on the matched rims (see
    :func:`icp_refine`); ``dof`` picks 4 (yaw + translation, default) or 6 (full
    rigid). Each rim is gathered per hole from its ellipse grown by ``rim_inflate``
    (falling back to a ``rim_radius`` circle without ellipse data). ``z_offset``
    slides the transformed tube along the aerial up-axis
    (metres, applied after any ``refine``) to set how deep its roof sits relative to
    the aerial surface -- useful when the matched skylights disagree on depth. With
    ``show_rims`` the rim patches GICP operates on are plotted (see
    :func:`visualize_rims`) before the cloud is written; with ``refine`` this shows
    the rims before and after refinement side by side.
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
    if dof not in (4, 6):
        raise ValueError("dof must be 4 (yaw + translation) or 6 (full rigid)")

    landmark = transform
    if refine:
        transform = icp_refine(aerial_path, tube_path, transform,
                               radius=rim_radius, rim_height=rim_height,
                               rim_inflate=rim_inflate, dof=dof)
    if show_rims:
        visualize_rims(aerial_path, tube_path, landmark,
                       refined=transform if refine else None,
                       radius=rim_radius, rim_height=rim_height, rim_inflate=rim_inflate)
    M = transform.array
    if z_offset:
        up = normalize(np.asarray(transform.aerial_up, dtype=np.float64))
        M = M.copy()
        M[:3, 3] += z_offset * up

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
