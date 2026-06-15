"""Chunked reading of ``.laz`` / ``.las`` point clouds via :mod:`laspy`.

Reading in chunks keeps memory bounded for very large clouds: only
``chunk_size`` points are materialised at a time instead of the whole file.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import laspy
import numpy as np

DEFAULT_CHUNK_SIZE = 5_000_000


class LazChunkReader:
    """Stream the XYZ + intensity of a LAZ/LAS file as ``(k, 4)`` float32 chunks.

    Usage::

        with LazChunkReader("scan.laz") as reader:
            print(reader.point_count)
            for chunk in reader.chunks():
                ...  # chunk is np.ndarray, shape (k, 4), columns x y z intensity
    """

    def __init__(self, path: str | Path, chunk_size: int = DEFAULT_CHUNK_SIZE) -> None:
        self.path = Path(path)
        if chunk_size <= 0:
            raise ValueError(f"chunk_size must be positive, got {chunk_size}")
        self.chunk_size = chunk_size
        self._reader: laspy.LasReader | None = None

    def __enter__(self) -> "LazChunkReader":
        self._reader = laspy.open(self.path)
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

    def chunks(self) -> Iterator[np.ndarray]:
        """Yield ``(k, 4)`` float32 arrays of columns ``x, y, z, intensity``.

        Coordinates are the *scaled* values (laspy applies scale + offset) and
        are downcast to float32. Intensity (LAZ ``uint16``) is cast to float32.
        """
        reader = self._require_reader()
        for points in reader.chunk_iterator(self.chunk_size):
            out = np.empty((len(points), 4), dtype=np.float32)
            out[:, 0] = points.x
            out[:, 1] = points.y
            out[:, 2] = points.z
            out[:, 3] = points.intensity
            yield out

    def _require_reader(self) -> laspy.LasReader:
        if self._reader is None:
            raise RuntimeError(
                "LazChunkReader must be used as a context manager "
                "(`with LazChunkReader(...) as reader:`)."
            )
        return self._reader
