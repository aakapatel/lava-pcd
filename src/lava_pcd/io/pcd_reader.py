"""Chunked reading of the binary PCD files written by :mod:`lava_pcd`.

Only the dialect this package emits is supported: ``DATA binary`` with every
field stored as a single float32 (``SIZE 4``, ``TYPE F``, ``COUNT 1``). That
covers ``x y z`` plus any of ``intensity`` and the packed ``rgb`` field. Other
PCD encodings (ASCII, binary_compressed, non-float fields) raise a clear error.

Reading is chunked so memory stays bounded for very large clouds, mirroring the
LAZ reader. The local-origin shift recorded in the ``# LAVA_PCD_ORIGIN`` header
comment is recovered so it can be carried through to derived outputs.
"""

from __future__ import annotations

from pathlib import Path
from types import TracebackType
from typing import Iterator

import numpy as np

_ORIGIN_TAG = "# LAVA_PCD_ORIGIN"


def pcd_fields_to_columns(fields: tuple[str, ...]) -> tuple[str, ...]:
    """Map PCD field names to internal columns (packed ``rgb`` -> ``r g b``)."""
    out: list[str] = []
    for name in fields:
        if name == "rgb":
            out.extend(("r", "g", "b"))
        else:
            out.append(name)
    return tuple(out)


def unpack_columns(data: np.ndarray, fields: tuple[str, ...]) -> np.ndarray:
    """Convert a PCD-layout ``(k, len(fields))`` array to internal columns.

    The inverse of :func:`lava_pcd.io.pcd_writer.pack_columns`: ``x y z`` and
    ``intensity`` pass through, while a packed ``rgb`` float is unpacked into
    three 0..255 ``r g b`` columns (so averaging during downsampling is valid).
    """
    out_cols: list[np.ndarray] = []
    for i, name in enumerate(fields):
        col = np.ascontiguousarray(data[:, i], dtype=np.float32)
        if name == "rgb":
            packed = col.view(np.uint32)
            out_cols.append(((packed >> 16) & 0xFF).astype(np.float32))
            out_cols.append(((packed >> 8) & 0xFF).astype(np.float32))
            out_cols.append((packed & 0xFF).astype(np.float32))
        else:
            out_cols.append(col)
    return np.ascontiguousarray(np.column_stack(out_cols), dtype=np.float32)


def _parse_origin(comment: str) -> tuple[float, float, float]:
    parts = comment[len(_ORIGIN_TAG):].split()
    if len(parts) == 3:
        try:
            return (float(parts[0]), float(parts[1]), float(parts[2]))
        except ValueError:
            pass
    return (0.0, 0.0, 0.0)


class BinaryPcdReader:
    """Stream a binary float32 PCD file as ``(k, len(fields))`` chunks.

    Usage::

        with BinaryPcdReader("cloud.pcd") as reader:
            print(reader.fields, reader.point_count, reader.origin)
            for chunk in reader.chunks():
                ...  # chunk is np.ndarray, shape (k, len(reader.fields))
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.fields: tuple[str, ...] = ()
        self.point_count = 0
        self.origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
        self._point_step = 0
        self._data_offset = 0
        self._fh = None

    def __enter__(self) -> "BinaryPcdReader":
        if not self.path.exists():
            raise FileNotFoundError(f"input file not found: {self.path}")
        self._fh = open(self.path, "rb")
        self._read_header()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None

    def _read_header(self) -> None:
        assert self._fh is not None
        meta: dict[str, str] = {}
        while True:
            raw = self._fh.readline()
            if not raw:
                raise ValueError(f"{self.path.name}: reached EOF before DATA line")
            line = raw.decode("ascii", "replace").rstrip("\n")
            if line.startswith(_ORIGIN_TAG):
                self.origin = _parse_origin(line)
                continue
            if line.startswith("#") or not line.strip():
                continue
            key, _, value = line.partition(" ")
            meta[key.upper()] = value.strip()
            if key.upper() == "DATA":
                break

        if meta.get("DATA", "").lower() != "binary":
            raise ValueError(
                f"{self.path.name}: only 'DATA binary' PCD is supported, "
                f"got DATA {meta.get('DATA')!r}"
            )
        fields = tuple(meta.get("FIELDS", "").split())
        sizes = meta.get("SIZE", "").split()
        types = meta.get("TYPE", "").split()
        counts = meta.get("COUNT", "").split() or ["1"] * len(fields)
        if not fields:
            raise ValueError(f"{self.path.name}: missing FIELDS in header")
        if not (all(s == "4" for s in sizes) and all(t == "F" for t in types)
                and all(c == "1" for c in counts)):
            raise ValueError(
                f"{self.path.name}: only all-float32 (SIZE 4 / TYPE F / COUNT 1) "
                "PCD files are supported"
            )
        self.fields = fields
        self.point_count = int(meta.get("POINTS", meta.get("WIDTH", "0")))
        self._point_step = len(fields) * 4
        self._data_offset = self._fh.tell()

    def chunks(self, chunk_size: int = 5_000_000) -> Iterator[np.ndarray]:
        """Yield ``(k, len(fields))`` float32 arrays, ``chunk_size`` points each."""
        if self._fh is None:
            raise RuntimeError(
                "BinaryPcdReader must be used as a context manager "
                "(`with BinaryPcdReader(...) as reader:`)."
            )
        if chunk_size <= 0:
            raise ValueError(f"chunk_size must be positive, got {chunk_size}")
        ncols = len(self.fields)
        self._fh.seek(self._data_offset)
        remaining = self.point_count
        while remaining > 0:
            k = min(chunk_size, remaining)
            buf = self._fh.read(k * self._point_step)
            got = len(buf) // self._point_step
            if got == 0:
                break
            yield np.frombuffer(buf, dtype=np.float32, count=got * ncols).reshape(got, ncols)
            remaining -= got
