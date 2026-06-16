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

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

from lava_pcd.geometry import as_matrix, from_matrix, kabsch, transform_points
from lava_pcd.io.pcd_reader import BinaryPcdReader
from lava_pcd.io.pcd_writer import BinaryPcdWriter
from lava_pcd.register import Transform

DEFAULT_RIM_RADIUS = 15.0   # gather points within this of each skylight for ICP
DEFAULT_ICP_ITERS = 30


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


def _gather_near(
    input_path: Path, anchors: np.ndarray, radius: float, matrix: np.ndarray | None = None
) -> np.ndarray:
    """Stream a cloud and keep XYZ within ``radius`` of any anchor.

    If ``matrix`` is given, points are transformed before the distance test, but
    the **original** (untransformed) coordinates are returned (so ICP can keep
    re-applying an evolving transform to a fixed source set).
    """
    tree = cKDTree(anchors)
    kept: list[np.ndarray] = []
    with BinaryPcdReader(input_path) as reader:
        for chunk in reader.chunks():
            xyz = chunk[:, :3].astype(np.float64)
            probe = transform_points(xyz, matrix) if matrix is not None else xyz
            near = tree.query_ball_point(probe, radius, workers=-1)
            mask = np.array([len(n) > 0 for n in near], dtype=bool)
            if mask.any():
                kept.append(xyz[mask])
    return np.concatenate(kept) if kept else np.empty((0, 3), dtype=np.float64)


def icp_refine(
    aerial_path: str | Path,
    tube_path: str | Path,
    transform: Transform,
    radius: float = DEFAULT_RIM_RADIUS,
    max_iters: int = DEFAULT_ICP_ITERS,
    tol: float = 1e-4,
) -> Transform:
    """Refine ``transform`` with point-to-point ICP on the matched skylight rims.

    Only points within ``radius`` of the matched rim centres (``transform.anchors``,
    in the aerial frame) are used from each cloud. Returns a new :class:`Transform`
    with the refined matrix and a recomputed rim RMS; the original landmark
    correspondences/warnings are carried over.
    """
    anchors = np.asarray(transform.anchors, dtype=np.float64)
    if len(anchors) == 0:
        return transform
    M = transform.array

    aerial_rim = _gather_near(Path(aerial_path), anchors, radius)
    tube_rim = _gather_near(Path(tube_path), anchors, radius, matrix=M)
    if len(aerial_rim) < 3 or len(tube_rim) < 3:
        return transform  # not enough overlap to refine; keep the landmark fit

    tree = cKDTree(aerial_rim)
    for _ in range(max_iters):
        pred = transform_points(tube_rim, M)
        dist, idx = tree.query(pred, k=1, workers=-1)
        keep = dist <= radius
        if keep.sum() < 3:
            break
        R, t = kabsch(tube_rim[keep], aerial_rim[idx[keep]])
        M_new = as_matrix(R, t)
        if np.linalg.norm(M_new - M) < tol:
            M = M_new
            break
        M = M_new

    pred = transform_points(tube_rim, M)
    dist, _ = tree.query(pred, k=1, workers=-1)
    rms = float(np.sqrt(np.mean(np.minimum(dist, radius) ** 2)))
    R, t = from_matrix(M)
    return Transform(
        matrix=[list(map(float, row)) for row in M],
        inliers=transform.inliers,
        rms=rms,
        mode=transform.mode + "+icp",
        aerial_origin=transform.aerial_origin,
        tube_origin=transform.tube_origin,
        anchors=transform.anchors,
        warnings=transform.warnings,
    )


def merge_clouds(
    aerial_path: str | Path,
    tube_path: str | Path,
    output_path: str | Path,
    transform: Transform,
    refine: bool = False,
    source_field: bool = False,
    rim_radius: float = DEFAULT_RIM_RADIUS,
    chunk_size: int = 5_000_000,
    show_progress: bool = True,
) -> MergeResult:
    """Transform the tube into the aerial frame and write a single merged ``.pcd``.

    The output is in **aerial-local** coordinates (the aerial origin), since the
    transform maps tube-local -> aerial-local. The two clouds usually carry
    different fields (aerial RGB vs tube intensity), so the merged cloud keeps
    only ``x y z`` plus, when ``source_field`` is set, a ``source`` channel
    (0 = aerial, 1 = tube) to colour by origin in a viewer.
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
        transform = icp_refine(aerial_path, tube_path, transform, radius=rim_radius)
    M = transform.array

    fields = ("x", "y", "z", "source") if source_field else ("x", "y", "z")

    with BinaryPcdReader(aerial_path) as ra, BinaryPcdReader(tube_path) as rb:
        origin = ra.origin
        a_total, b_total = ra.point_count, rb.point_count
    max_points = a_total + b_total

    bar = tqdm(total=max_points, unit="pts", unit_scale=True, desc="merge",
               disable=not show_progress)
    with BinaryPcdWriter(output_path, max_points=max_points, fields=fields,
                         origin=origin) as writer, bar:
        # Aerial cloud: identity (already in the target frame), source 0.
        with BinaryPcdReader(aerial_path) as ra:
            for chunk in ra.chunks(chunk_size):
                writer.write_chunk(_emit(chunk[:, :3].astype(np.float64), 0, source_field))
                bar.update(len(chunk))
        # Tube cloud: transform into the aerial frame, source 1.
        with BinaryPcdReader(tube_path) as rb:
            for chunk in rb.chunks(chunk_size):
                xyz = transform_points(chunk[:, :3].astype(np.float64), M)
                writer.write_chunk(_emit(xyz, 1, source_field))
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


def _emit(xyz: np.ndarray, source: int, source_field: bool) -> np.ndarray:
    """Build the PCD-layout chunk: xyz (+ a constant source column)."""
    if source_field:
        col = np.full((len(xyz), 1), float(source), dtype=np.float64)
        xyz = np.hstack([xyz, col])
    return np.ascontiguousarray(xyz, dtype=np.float32)
