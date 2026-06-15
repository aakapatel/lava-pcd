"""Streaming writer for binary PCD (Point Cloud Data) files.

Writes a PCD v0.7 file with fields ``x y z intensity`` (all float32). Open3D's
own PCD writer only supports points/colors/normals, so it cannot emit an
``intensity`` field that ``pcl_viewer`` recognises -- hence this direct writer.
The total point count is written into the header up front, so it must be known
before the data is streamed (the LAZ header provides it).
"""

from __future__ import annotations

from pathlib import Path
from types import TracebackType

import numpy as np

FIELDS = ("x", "y", "z", "intensity")


def _header(num_points: int, origin: tuple[float, float, float]) -> str:
    fields = " ".join(FIELDS)
    sizes = " ".join("4" for _ in FIELDS)
    types = " ".join("F" for _ in FIELDS)
    counts = " ".join("1" for _ in FIELDS)
    # Record the local-origin shift (global = local + origin) as a comment so the
    # cloud can be georeferenced back. Kept out of VIEWPOINT so viewers render the
    # small local coordinates without float jitter.
    origin_comment = (
        f"# LAVA_PCD_ORIGIN {origin[0]!r} {origin[1]!r} {origin[2]!r}\n"
    )
    return (
        "# .PCD v0.7 - Point Cloud Data file format\n"
        + origin_comment
        + "VERSION 0.7\n"
        f"FIELDS {fields}\n"
        f"SIZE {sizes}\n"
        f"TYPE {types}\n"
        f"COUNT {counts}\n"
        f"WIDTH {num_points}\n"
        "HEIGHT 1\n"
        "VIEWPOINT 0 0 0 1 0 0 0\n"
        f"POINTS {num_points}\n"
        "DATA binary\n"
    )


class BinaryPcdWriter:
    """Stream ``(k, 4)`` float32 chunks into a binary PCD file.

    Usage::

        with BinaryPcdWriter("out.pcd", num_points=N) as writer:
            for chunk in chunks:           # each (k, 4) float32: x y z intensity
                writer.write_chunk(chunk)

    ``num_points`` must equal the total number of points written, since it is
    baked into the header.
    """

    def __init__(
        self,
        path: str | Path,
        num_points: int,
        origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> None:
        if num_points < 0:
            raise ValueError(f"num_points must be non-negative, got {num_points}")
        self.path = Path(path)
        self.num_points = num_points
        self.origin = origin
        self._written = 0
        self._fh = None

    def __enter__(self) -> "BinaryPcdWriter":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(self.path, "wb")
        self._fh.write(_header(self.num_points, self.origin).encode("ascii"))
        return self

    def write_chunk(self, points: np.ndarray) -> None:
        """Append a ``(k, 4)`` array of ``x y z intensity`` (cast to float32)."""
        if self._fh is None:
            raise RuntimeError(
                "BinaryPcdWriter must be used as a context manager "
                "(`with BinaryPcdWriter(...) as writer:`)."
            )
        if points.ndim != 2 or points.shape[1] != len(FIELDS):
            raise ValueError(
                f"expected a (k, {len(FIELDS)}) array, got shape {points.shape}"
            )
        # Ensure C-contiguous float32 so .tobytes() matches the declared layout.
        data = np.ascontiguousarray(points, dtype=np.float32)
        self._fh.write(data.tobytes())
        self._written += len(data)

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None
        # Only validate the count on a clean exit; don't mask an in-flight error.
        if exc_type is None and self._written != self.num_points:
            raise ValueError(
                f"wrote {self._written} points but header declared "
                f"{self.num_points}; PCD file is inconsistent."
            )
