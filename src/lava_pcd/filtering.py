"""Outlier filtering for binary ``.pcd`` clouds.

Two configurable methods, mirroring the familiar PCL / Open3D semantics:

* **radius** -- remove points with too few neighbours within a fixed radius
  (``remove_radius_outlier``): a point is kept when at least ``min_neighbors``
  points (including itself) lie within ``radius``.
* **statistical** -- remove points whose mean distance to their ``k`` nearest
  neighbours is an outlier (``remove_statistical_outlier``): a point is kept
  when that mean distance is ``<= global_mean + std_ratio * global_std``.

Both build a KD-tree over the XYZ of the whole cloud, so the cloud is loaded
into memory (neighbour queries cross chunk boundaries -- they cannot be
streamed). Filter large clouds after cropping/downsampling. Every field and the
local origin (the ``# LAVA_PCD_ORIGIN`` comment) is preserved on the survivors.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

from lava_pcd.io.pcd_reader import BinaryPcdReader
from lava_pcd.io.pcd_writer import BinaryPcdWriter

DEFAULT_RADIUS = 0.5
DEFAULT_MIN_NEIGHBORS = 5
DEFAULT_K = 20
DEFAULT_STD_RATIO = 2.0
_QUERY_BLOCK = 1_000_000  # points per neighbour-query batch (bounds memory)


@dataclass
class FilterResult:
    """Outcome of a :func:`filter_pcd` call."""

    point_count: int
    source_count: int
    removed: int
    method: str
    origin: tuple[float, float, float]
    output_path: Path
    fields: tuple[str, ...] = ()
    params: dict = field(default_factory=dict)

    def __int__(self) -> int:
        return self.point_count


def radius_outlier_mask(
    xyz: np.ndarray,
    radius: float = DEFAULT_RADIUS,
    min_neighbors: int = DEFAULT_MIN_NEIGHBORS,
    show_progress: bool = True,
) -> np.ndarray:
    """Boolean keep-mask: ``True`` where >= ``min_neighbors`` points (incl. self)
    lie within ``radius``."""
    if radius <= 0:
        raise ValueError(f"radius must be > 0, got {radius}")
    if min_neighbors < 1:
        raise ValueError(f"min_neighbors must be >= 1, got {min_neighbors}")
    n = len(xyz)
    tree = cKDTree(xyz)
    counts = np.empty(n, dtype=np.int64)
    bar = tqdm(total=n, unit="pts", unit_scale=True, desc="radius filter",
               disable=not show_progress)
    with bar:
        for s in range(0, n, _QUERY_BLOCK):
            e = min(s + _QUERY_BLOCK, n)
            counts[s:e] = tree.query_ball_point(
                xyz[s:e], radius, return_length=True, workers=-1
            )
            bar.update(e - s)
    return counts >= min_neighbors


def statistical_outlier_mask(
    xyz: np.ndarray,
    k: int = DEFAULT_K,
    std_ratio: float = DEFAULT_STD_RATIO,
    show_progress: bool = True,
) -> np.ndarray:
    """Boolean keep-mask: ``True`` where the mean distance to the ``k`` nearest
    neighbours is ``<= mean + std_ratio * std`` over the whole cloud."""
    if std_ratio <= 0:
        raise ValueError(f"std_ratio must be > 0, got {std_ratio}")
    n = len(xyz)
    if n <= 1:
        return np.ones(n, dtype=bool)
    k = min(k, n - 1)  # can't ask for more neighbours than exist
    tree = cKDTree(xyz)
    mean_dist = np.empty(n, dtype=np.float64)
    bar = tqdm(total=n, unit="pts", unit_scale=True, desc="statistical filter",
               disable=not show_progress)
    with bar:
        for s in range(0, n, _QUERY_BLOCK):
            e = min(s + _QUERY_BLOCK, n)
            # k+1 because the first neighbour is the point itself (distance 0).
            dist, _ = tree.query(xyz[s:e], k=k + 1, workers=-1)
            mean_dist[s:e] = dist[:, 1:].mean(axis=1)
            bar.update(e - s)
    threshold = mean_dist.mean() + std_ratio * mean_dist.std()
    return mean_dist <= threshold


def _load_cloud(
    input_path: Path,
) -> tuple[np.ndarray, tuple[str, ...], tuple[float, float, float]]:
    with BinaryPcdReader(input_path) as reader:
        fields = reader.fields
        origin = reader.origin
        parts = [chunk for chunk in reader.chunks()]
    if parts:
        data = np.concatenate(parts)
    else:
        data = np.empty((0, len(fields)), dtype=np.float32)
    return data, fields, origin


def filter_pcd(
    input_path: str | Path,
    output_path: str | Path,
    method: str = "statistical",
    radius: float = DEFAULT_RADIUS,
    min_neighbors: int = DEFAULT_MIN_NEIGHBORS,
    k: int = DEFAULT_K,
    std_ratio: float = DEFAULT_STD_RATIO,
    show_progress: bool = True,
) -> FilterResult:
    """Remove outliers from ``input_path`` and write the survivors to ``output_path``.

    ``method`` is ``"radius"`` (uses ``radius`` + ``min_neighbors``) or
    ``"statistical"`` (uses ``k`` + ``std_ratio``). All fields and the local
    origin are preserved. Returns a :class:`FilterResult`.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    if input_path.suffix.lower() != ".pcd" or output_path.suffix.lower() != ".pcd":
        raise ValueError("filter works on .pcd input and output")
    if input_path.resolve() == output_path.resolve():
        raise ValueError("input and output must be different files")
    if method not in ("radius", "statistical"):
        raise ValueError(f"method must be 'radius' or 'statistical', got {method!r}")

    data, fields, origin = _load_cloud(input_path)
    source_count = len(data)
    xyz = np.ascontiguousarray(data[:, :3], dtype=np.float64)

    if method == "radius":
        mask = radius_outlier_mask(xyz, radius, min_neighbors, show_progress)
        params = {"radius": radius, "min_neighbors": min_neighbors}
    else:
        mask = statistical_outlier_mask(xyz, k, std_ratio, show_progress)
        params = {"k": k, "std_ratio": std_ratio}

    kept = data[mask]
    with BinaryPcdWriter(
        output_path, max_points=len(kept), fields=fields, origin=origin
    ) as writer:
        if len(kept):
            writer.write_chunk(np.ascontiguousarray(kept, dtype=np.float32))

    return FilterResult(
        point_count=len(kept),
        source_count=source_count,
        removed=source_count - len(kept),
        method=method,
        origin=origin,
        output_path=output_path,
        fields=fields,
        params=params,
    )
