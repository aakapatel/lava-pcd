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
    ) -> None:
        self.path = Path(path)
        if chunk_size <= 0:
            raise ValueError(f"chunk_size must be positive, got {chunk_size}")
        self.chunk_size = chunk_size
        self.parallel = parallel
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

    def chunks(
        self, origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    ) -> Iterator[np.ndarray]:
        """Yield ``(k, 4)`` float32 arrays of columns ``x, y, z, intensity``.

        ``origin`` is subtracted from the scaled coordinates **in float64**
        before downcasting to float32, preserving precision for large
        (e.g. UTM) coordinates. Intensity (LAZ ``uint16``) is cast to float32.
        """
        reader = self._require_reader()
        ox, oy, oz = origin
        for points in reader.chunk_iterator(self.chunk_size):
            out = np.empty((len(points), 4), dtype=np.float32)
            out[:, 0] = np.asarray(points.x, dtype=np.float64) - ox
            out[:, 1] = np.asarray(points.y, dtype=np.float64) - oy
            out[:, 2] = np.asarray(points.z, dtype=np.float64) - oz
            out[:, 3] = points.intensity
            yield out

    def _require_reader(self) -> laspy.LasReader:
        if self._reader is None:
            raise RuntimeError(
                "LazChunkReader must be used as a context manager "
                "(`with LazChunkReader(...) as reader:`)."
            )
        return self._reader
