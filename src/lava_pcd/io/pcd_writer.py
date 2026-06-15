"""Streaming writer for binary PCD (Point Cloud Data) files.

Writes a PCD v0.7 file from internal per-point columns (``x y z`` plus any of
``intensity`` and the colour channels ``r g b``). Open3D's own PCD writer only
supports points/colors/normals and cannot emit an ``intensity`` field that
``pcl_viewer`` recognises -- hence this direct writer.

Colour is emitted as PCL's packed ``rgb`` field: a single float32 whose bits
hold ``0x00RRGGBB``. ``pcl_viewer`` colours by this field automatically. The
total point count is written into the header up front, so it must be known
before the data is streamed (the LAZ header provides it, or downsampling does).
"""

from __future__ import annotations

from pathlib import Path
from types import TracebackType

import numpy as np


def pcd_fields_for(columns: tuple[str, ...]) -> tuple[str, ...]:
    """Map internal columns to PCD field names (``r g b`` -> packed ``rgb``)."""
    out: list[str] = []
    for name in columns:
        if name == "r":
            out.append("rgb")
        elif name in ("g", "b"):
            continue  # folded into the single rgb field
        else:
            out.append(name)
    return tuple(out)


def pack_columns(data: np.ndarray, columns: tuple[str, ...]) -> np.ndarray:
    """Convert internal ``(k, len(columns))`` data to the PCD output layout.

    ``x y z`` and ``intensity`` pass through; ``r g b`` (each 0..255) are packed
    into one float32 ``rgb`` value. Returns a C-contiguous float32 array.
    """
    idx = {name: i for i, name in enumerate(columns)}
    out_cols: list[np.ndarray] = []
    for name in columns:
        if name in ("g", "b"):
            continue
        if name == "r":
            r = np.clip(data[:, idx["r"]], 0, 255).astype(np.uint32)
            g = np.clip(data[:, idx["g"]], 0, 255).astype(np.uint32)
            b = np.clip(data[:, idx["b"]], 0, 255).astype(np.uint32)
            packed = (r << 16) | (g << 8) | b
            out_cols.append(packed.view(np.float32))
        else:
            out_cols.append(data[:, idx[name]].astype(np.float32))
    return np.ascontiguousarray(np.column_stack(out_cols), dtype=np.float32)


def _header(
    count_field: str, fields: tuple[str, ...], origin: tuple[float, float, float]
) -> tuple[str, int, int]:
    """Build the PCD header. Returns (header, width_offset, points_offset).

    ``count_field`` is the point count rendered to a fixed width; the two byte
    offsets locate it inside WIDTH/POINTS so the value can be rewritten on close
    without changing the header length (the binary data offset stays valid).
    """
    field_str = " ".join(fields)
    sizes = " ".join("4" for _ in fields)
    types = " ".join("F" for _ in fields)
    counts = " ".join("1" for _ in fields)
    # Record the local-origin shift (global = local + origin) as a comment so the
    # cloud can be georeferenced back. Kept out of VIEWPOINT so viewers render the
    # small local coordinates without float jitter.
    origin_comment = f"# LAVA_PCD_ORIGIN {origin[0]!r} {origin[1]!r} {origin[2]!r}\n"
    header = (
        "# .PCD v0.7 - Point Cloud Data file format\n"
        + origin_comment
        + "VERSION 0.7\n"
        f"FIELDS {field_str}\n"
        f"SIZE {sizes}\n"
        f"TYPE {types}\n"
        f"COUNT {counts}\n"
        f"WIDTH {count_field}\n"
        "HEIGHT 1\n"
        "VIEWPOINT 0 0 0 1 0 0 0\n"
        f"POINTS {count_field}\n"
        "DATA binary\n"
    )
    width_off = header.index("WIDTH ") + len("WIDTH ")
    points_off = header.index("POINTS ") + len("POINTS ")
    return header, width_off, points_off


class BinaryPcdWriter:
    """Stream float32 chunks (already in PCD column layout) into a binary PCD.

    ``fields`` are the PCD field names (e.g. ``("x", "y", "z", "rgb")``).
    Each chunk passed to :meth:`write_chunk` must have ``len(fields)`` columns.

    ``max_points`` is an upper bound on the number of points (e.g. the LAZ
    header count). The actual count -- which may be smaller if points were
    filtered out -- is written into the header on close. The header reserves a
    fixed-width count field so this rewrite never shifts the binary data.
    """

    def __init__(
        self,
        path: str | Path,
        max_points: int,
        fields: tuple[str, ...] = ("x", "y", "z", "intensity"),
        origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> None:
        if max_points < 0:
            raise ValueError(f"max_points must be non-negative, got {max_points}")
        self.path = Path(path)
        self.max_points = max_points
        self.fields = fields
        self.origin = origin
        self._width = len(str(max_points))
        self._written = 0
        self._width_off = 0
        self._points_off = 0
        self._fh = None

    def __enter__(self) -> "BinaryPcdWriter":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(self.path, "wb")
        # Reserve the count field at full width; finalised in __exit__.
        placeholder = str(self.max_points).rjust(self._width)
        header, self._width_off, self._points_off = _header(
            placeholder, self.fields, self.origin
        )
        self._fh.write(header.encode("ascii"))
        return self

    def write_chunk(self, points: np.ndarray) -> None:
        """Append a ``(k, len(fields))`` float32 array (already PCD-laid-out)."""
        if self._fh is None:
            raise RuntimeError(
                "BinaryPcdWriter must be used as a context manager "
                "(`with BinaryPcdWriter(...) as writer:`)."
            )
        if points.ndim != 2 or points.shape[1] != len(self.fields):
            raise ValueError(
                f"expected a (k, {len(self.fields)}) array, got shape {points.shape}"
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
        if self._fh is None:
            return
        # On a clean exit, write the true count back into WIDTH/POINTS.
        if exc_type is None:
            if self._written > self.max_points:
                self._fh.close()
                self._fh = None
                raise ValueError(
                    f"wrote {self._written} points but reserved only "
                    f"{self.max_points}; header count field too narrow."
                )
            value = str(self._written).rjust(self._width).encode("ascii")
            self._fh.seek(self._width_off)
            self._fh.write(value)
            self._fh.seek(self._points_off)
            self._fh.write(value)
        self._fh.close()
        self._fh = None
