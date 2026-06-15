"""Chunked reading of ``.laz`` / ``.las`` point clouds via :mod:`laspy`.

Reading in chunks keeps memory bounded for very large clouds: only
``chunk_size`` points are materialised at a time instead of the whole file.

A ``origin`` (local-coordinate shift) is subtracted from XYZ **in float64**
before the values are downcast to float32. This is essential for survey-grade
data (e.g. UTM coordinates in the millions of metres): float32 has only ~7
significant digits, so storing raw UTM values would lose centimetre/metre
precision. Subtracting a nearby origin keeps the stored values small and exact.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import laspy
import numpy as np

DEFAULT_CHUNK_SIZE = 5_000_000

# Internal per-point column layouts for each named schema. RGB channels are
# carried as separate 0..255 float columns (so averaging during downsampling is
# valid); they get packed into the single PCD ``rgb`` float field at write time.
_SCHEMAS: dict[str, tuple[str, ...]] = {
    "xyz": ("x", "y", "z"),
    "intensity": ("x", "y", "z", "intensity"),
    "rgb": ("x", "y", "z", "r", "g", "b"),
    "all": ("x", "y", "z", "intensity", "r", "g", "b"),
}


def in_bounds_mask(
    gx: np.ndarray,
    gy: np.ndarray,
    gz: np.ndarray,
    mins: tuple[float, float, float],
    maxs: tuple[float, float, float],
    eps: float = 1e-3,
) -> np.ndarray:
    """Boolean mask of points within the (mins, maxs) box (inclusive, +/-eps).

    LAZ files can contain invalid/sentinel points (e.g. raw coords at INT32
    limits) that fall outside the header's declared bounding box; these would
    otherwise wreck the cloud's scale in a viewer.
    """
    return (
        (gx >= mins[0] - eps) & (gx <= maxs[0] + eps)
        & (gy >= mins[1] - eps) & (gy <= maxs[1] + eps)
        & (gz >= mins[2] - eps) & (gz <= maxs[2] + eps)
    )


def _resolve_backend(parallel: bool) -> "laspy.LazBackend | None":
    """Pick a LAZ backend. Default to single-threaded lazrs.

    The parallel lazrs backend is known to panic on some files
    (``range end index ... out of range for slice of length 0``), so it is
    opt-in only.
    """
    backends = laspy.LazBackend
    if parallel and backends.LazrsParallel.is_available():
        return backends.LazrsParallel
    if backends.Lazrs.is_available():
        return backends.Lazrs
    return None  # let laspy choose whatever is available


class LazChunkReader:
    """Stream the XYZ + intensity of a LAZ/LAS file as ``(k, 4)`` float32 chunks.

    Usage::

        with LazChunkReader("scan.laz") as reader:
            print(reader.point_count, reader.header_offset)
            for chunk in reader.chunks(origin=reader.header_offset):
                ...  # chunk is np.ndarray, shape (k, 4): x y z intensity
    """

    def __init__(
        self,
        path: str | Path,
        chunk_size: int = DEFAULT_CHUNK_SIZE,
        parallel: bool = False,
        filter_bounds: bool = True,
    ) -> None:
        self.path = Path(path)
        if chunk_size <= 0:
            raise ValueError(f"chunk_size must be positive, got {chunk_size}")
        self.chunk_size = chunk_size
        self.parallel = parallel
        self.filter_bounds = filter_bounds
        self.dropped = 0  # points discarded as out-of-bounds during reading
        self._reader: laspy.LasReader | None = None

    def __enter__(self) -> "LazChunkReader":
        self._reader = laspy.open(self.path, laz_backend=_resolve_backend(self.parallel))
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def close(self) -> None:
        if self._reader is not None:
            self._reader.close()
            self._reader = None

    @property
    def point_count(self) -> int:
        """Total number of points, available without reading point data."""
        return int(self._require_reader().header.point_count)

    @property
    def header_offset(self) -> tuple[float, float, float]:
        """The LAS header XYZ offset -- a natural local origin near the data."""
        off = self._require_reader().header.offsets
        return (float(off[0]), float(off[1]), float(off[2]))

    @property
    def header_bounds(self) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        """The LAS header (mins, maxs) bounding box of valid coordinates."""
        h = self._require_reader().header
        return (
            (float(h.mins[0]), float(h.mins[1]), float(h.mins[2])),
            (float(h.maxs[0]), float(h.maxs[1]), float(h.maxs[2])),
        )

    @property
    def has_rgb(self) -> bool:
        """Whether the point format carries red/green/blue channels."""
        dims = self._require_reader().header.point_format.dimension_names
        return {"red", "green", "blue"}.issubset(set(dims))

    def resolve_fields(self, fields: str) -> tuple[str, ...]:
        """Map a schema name (incl. ``"auto"``) to internal column names.

        ``"auto"`` picks ``rgb`` when the file has colour, else ``intensity``
        (LiDAR intensity is meaningless/empty on many colourised products).
        """
        if fields == "auto":
            fields = "rgb" if self.has_rgb else "intensity"
        if fields not in _SCHEMAS:
            raise ValueError(
                f"unknown fields schema {fields!r}; "
                f"choose from {sorted(_SCHEMAS) + ['auto']}"
            )
        cols = _SCHEMAS[fields]
        if ("r" in cols) and not self.has_rgb:
            raise ValueError(
                f"schema {fields!r} needs RGB but this file has none "
                f"(point format {self._require_reader().header.point_format.id})"
            )
        return cols

    def chunks(
        self,
        columns: tuple[str, ...],
        origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> Iterator[np.ndarray]:
        """Yield ``(k, len(columns))`` float32 arrays for the given columns.

        ``origin`` is subtracted from the scaled coordinates **in float64**
        before downcasting to float32, preserving precision for large
        (e.g. UTM) coordinates. RGB channels are scaled from 16-bit to 0..255.
        When ``filter_bounds`` is set, points outside the header bounding box
        are dropped (counted in ``self.dropped``), so chunk sizes may vary.
        """
        reader = self._require_reader()
        ox, oy, oz = origin
        shifts = {"x": ox, "y": oy, "z": oz}
        mins, maxs = self.header_bounds
        self.dropped = 0
        for points in reader.chunk_iterator(self.chunk_size):
            gx = np.asarray(points.x, dtype=np.float64)
            gy = np.asarray(points.y, dtype=np.float64)
            gz = np.asarray(points.z, dtype=np.float64)
            if self.filter_bounds:
                mask = in_bounds_mask(gx, gy, gz, mins, maxs)
                self.dropped += int(len(gx) - mask.sum())
            else:
                mask = slice(None)
            globals_xyz = {"x": gx, "y": gy, "z": gz}
            k = int(mask.sum()) if self.filter_bounds else len(gx)
            out = np.empty((k, len(columns)), dtype=np.float32)
            for i, name in enumerate(columns):
                if name in ("x", "y", "z"):
                    out[:, i] = globals_xyz[name][mask] - shifts[name]
                elif name == "intensity":
                    out[:, i] = np.asarray(points.intensity)[mask]
                else:  # r, g, b: 16-bit -> 0..255
                    channel = {"r": "red", "g": "green", "b": "blue"}[name]
                    out[:, i] = np.asarray(getattr(points, channel))[mask] / 257.0
            yield out

    def _require_reader(self) -> laspy.LasReader:
        if self._reader is None:
            raise RuntimeError(
                "LazChunkReader must be used as a context manager "
                "(`with LazChunkReader(...) as reader:`)."
            )
        return self._reader
