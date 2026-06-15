"""High-level conversion routines."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from tqdm import tqdm

from lava_pcd.io.laz_reader import DEFAULT_CHUNK_SIZE, LazChunkReader
from lava_pcd.io.pcd_writer import BinaryPcdWriter

_LAZ_SUFFIXES = {".laz", ".las"}


@dataclass
class ConvertResult:
    """Outcome of a :func:`laz_to_pcd` call."""

    point_count: int
    origin: tuple[float, float, float]
    output_path: Path
    sidecar_path: Path | None

    def __int__(self) -> int:  # backwards-compatible: len-like usage
        return self.point_count


def _resolve_origin(
    origin: str | tuple[float, float, float],
    header_offset: tuple[float, float, float],
) -> tuple[float, float, float]:
    if origin == "header":
        return header_offset
    if origin == "none":
        return (0.0, 0.0, 0.0)
    if isinstance(origin, (tuple, list)) and len(origin) == 3:
        return (float(origin[0]), float(origin[1]), float(origin[2]))
    raise ValueError(
        f"origin must be 'header', 'none', or an (x, y, z) tuple, got {origin!r}"
    )


def laz_to_pcd(
    input_path: str | Path,
    output_path: str | Path,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    origin: str | tuple[float, float, float] = "header",
    parallel: bool = False,
    write_sidecar: bool = True,
    show_progress: bool = True,
) -> ConvertResult:
    """Convert a ``.laz`` / ``.las`` file to a binary ``.pcd`` (x y z intensity).

    Reads the input in chunks of ``chunk_size`` points to keep memory bounded,
    streaming each chunk into the output PCD.

    ``origin`` controls the local-coordinate shift subtracted from XYZ (in
    float64, before the float32 cast) so large survey coordinates keep their
    precision:

    * ``"header"`` (default) -- use the LAS header offset, a natural origin
      near the data.
    * ``"none"`` -- no shift (only safe for clouds already near the origin).
    * ``(x, y, z)`` -- an explicit origin.

    The chosen origin is recorded in the PCD header comment and, when
    ``write_sidecar`` is set, in a ``<output>.origin.json`` file, so points can
    be georeferenced back via ``global = local + origin``.

    Set ``parallel=True`` to try the multi-threaded LAZ backend (faster, but it
    panics on some files; off by default).

    Returns a :class:`ConvertResult`.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.exists():
        raise FileNotFoundError(f"input file not found: {input_path}")
    if input_path.suffix.lower() not in _LAZ_SUFFIXES:
        raise ValueError(
            f"expected a .laz or .las input, got '{input_path.suffix}' ({input_path})"
        )

    with LazChunkReader(input_path, chunk_size=chunk_size, parallel=parallel) as reader:
        total = reader.point_count
        resolved_origin = _resolve_origin(origin, reader.header_offset)
        with BinaryPcdWriter(output_path, num_points=total, origin=resolved_origin) as writer:
            progress = tqdm(
                total=total,
                unit="pts",
                unit_scale=True,
                desc=input_path.name,
                disable=not show_progress,
            )
            with progress:
                for chunk in reader.chunks(origin=resolved_origin):
                    writer.write_chunk(chunk)
                    progress.update(len(chunk))

    sidecar_path: Path | None = None
    if write_sidecar:
        sidecar_path = output_path.with_suffix(output_path.suffix + ".origin.json")
        sidecar_path.write_text(
            json.dumps(
                {
                    "source": str(input_path),
                    "point_count": total,
                    "origin_xyz": list(resolved_origin),
                    "note": "global_xyz = local_xyz + origin_xyz",
                },
                indent=2,
            )
        )

    return ConvertResult(
        point_count=total,
        origin=resolved_origin,
        output_path=output_path,
        sidecar_path=sidecar_path,
    )
