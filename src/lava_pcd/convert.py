"""High-level conversion routines."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from tqdm import tqdm

from lava_pcd.io.laz_reader import DEFAULT_CHUNK_SIZE, LazChunkReader
from lava_pcd.io.pcd_reader import BinaryPcdReader, pcd_fields_to_columns, unpack_columns
from lava_pcd.io.pcd_writer import BinaryPcdWriter, pack_columns, pcd_fields_for
from lava_pcd.voxel import VoxelDownsampler

_LAZ_SUFFIXES = {".laz", ".las"}


@dataclass
class ConvertResult:
    """Outcome of a :func:`laz_to_pcd` call."""

    point_count: int
    origin: tuple[float, float, float]
    output_path: Path
    sidecar_path: Path | None
    source_count: int = 0
    voxel_size: float = 0.0
    fields: tuple[str, ...] = ()
    dropped: int = 0
    reproject: str | None = None

    def __int__(self) -> int:  # backwards-compatible: len-like usage
        return self.point_count


@dataclass
class DownsampleResult:
    """Outcome of a :func:`downsample_pcd` call."""

    point_count: int
    source_count: int
    voxel_size: float
    origin: tuple[float, float, float]
    output_path: Path
    fields: tuple[str, ...] = ()

    def __int__(self) -> int:
        return self.point_count


def _resolve_origin(
    origin: str | tuple[float, float, float],
    default_origin: tuple[float, float, float],
) -> tuple[float, float, float]:
    if origin == "header":
        return default_origin
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
    voxel_size: float = 0.0,
    fields: str = "auto",
    reproject: str | None = None,
    parallel: bool = False,
    filter_bounds: bool = True,
    write_sidecar: bool = True,
    show_progress: bool = True,
) -> ConvertResult:
    """Convert a ``.laz`` / ``.las`` file to a binary ``.pcd``.

    Reads the input in chunks of ``chunk_size`` points to keep memory bounded,
    streaming each chunk into the output PCD.

    ``fields`` selects which attributes to export:

    * ``"auto"`` (default) -- RGB if the file has colour, else intensity.
    * ``"xyz"`` / ``"intensity"`` / ``"rgb"`` / ``"all"`` -- explicit. ``rgb``
      is written as PCL's packed ``rgb`` float field (coloured by ``pcl_viewer``).

    ``voxel_size`` > 0 voxel-downsamples the cloud (cubic voxels of that edge
    length, in coordinate units) before writing; each output point is the
    centroid of its voxel. ``voxel_size`` == 0 disables downsampling.

    ``origin`` controls the local-coordinate shift subtracted from XYZ (in
    float64, before the float32 cast) so large survey coordinates keep their
    precision:

    * ``"header"`` (default) -- a natural origin near the data (the LAS header
      offset, or the reprojected header-bbox corner when ``reproject`` is set).
    * ``"none"`` -- no shift (only safe for clouds already near the origin).
    * ``(x, y, z)`` -- an explicit origin (in the output CRS).

    ``reproject`` is an optional target CRS (e.g. ``"EPSG:32627"``). When given,
    X/Y are transformed from the file's CRS to it (Z unchanged), per chunk,
    before the origin shift.

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

    if voxel_size < 0:
        raise ValueError(f"voxel_size must be >= 0, got {voxel_size}")

    with LazChunkReader(
        input_path, chunk_size=chunk_size, parallel=parallel,
        filter_bounds=filter_bounds, reproject=reproject,
    ) as reader:
        source_count = reader.point_count
        resolved_origin = _resolve_origin(origin, reader.suggested_origin)
        columns = reader.resolve_fields(fields)
        pcd_fields = pcd_fields_for(columns)
        progress = tqdm(
            total=source_count,
            unit="pts",
            unit_scale=True,
            desc=input_path.name,
            disable=not show_progress,
        )

        def consume(chunk: "object", prev_dropped: int) -> int:
            # Advance the bar by input points consumed (kept + newly dropped).
            advance = len(chunk) + (reader.dropped - prev_dropped)
            progress.update(advance)
            return reader.dropped

        if voxel_size > 0:
            # Aggregate voxels across all chunks; the final count is only known
            # once everything has been read, so we write after accumulating.
            downsampler = VoxelDownsampler(voxel_size)
            prev_dropped = 0
            with progress:
                for chunk in reader.chunks(columns, origin=resolved_origin):
                    downsampler.add(chunk)
                    prev_dropped = consume(chunk, prev_dropped)
            points = downsampler.result()
            written = len(points)
            with BinaryPcdWriter(
                output_path, max_points=written, fields=pcd_fields, origin=resolved_origin
            ) as writer:
                writer.write_chunk(pack_columns(points, columns))
        else:
            prev_dropped = 0
            with BinaryPcdWriter(
                output_path, max_points=source_count, fields=pcd_fields,
                origin=resolved_origin,
            ) as writer:
                with progress:
                    for chunk in reader.chunks(columns, origin=resolved_origin):
                        writer.write_chunk(pack_columns(chunk, columns))
                        prev_dropped = consume(chunk, prev_dropped)
                written = writer._written

        dropped = reader.dropped

    sidecar_path: Path | None = None
    if write_sidecar:
        sidecar_path = output_path.with_suffix(output_path.suffix + ".origin.json")
        sidecar_path.write_text(
            json.dumps(
                {
                    "source": str(input_path),
                    "source_point_count": source_count,
                    "point_count": written,
                    "dropped_out_of_bounds": dropped,
                    "voxel_size": voxel_size,
                    "fields": list(pcd_fields),
                    "reproject": reproject,
                    "origin_xyz": list(resolved_origin),
                    "note": "global_xyz = local_xyz + origin_xyz (origin in output CRS)",
                },
                indent=2,
            )
        )

    return ConvertResult(
        point_count=written,
        origin=resolved_origin,
        output_path=output_path,
        sidecar_path=sidecar_path,
        source_count=source_count,
        voxel_size=voxel_size,
        fields=pcd_fields,
        dropped=dropped,
        reproject=reproject,
    )


def downsample_pcd(
    input_path: str | Path,
    output_path: str | Path,
    voxel_size: float,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    show_progress: bool = True,
) -> DownsampleResult:
    """Voxel-downsample an existing binary ``.pcd`` into a new ``.pcd``.

    Reads the input in chunks (memory stays bounded by the *downsampled* size),
    aggregates points into cubic voxels of edge length ``voxel_size`` -- each
    output point is the centroid of its voxel, with any ``intensity`` and colour
    averaged the same way -- and writes the result. The input's fields and its
    local-origin shift (the ``# LAVA_PCD_ORIGIN`` header comment) are preserved.

    Returns a :class:`DownsampleResult`.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    if input_path.suffix.lower() != ".pcd":
        raise ValueError(
            f"expected a .pcd input, got '{input_path.suffix}' ({input_path})"
        )
    if voxel_size <= 0:
        raise ValueError(f"voxel_size must be > 0, got {voxel_size}")

    with BinaryPcdReader(input_path) as reader:
        if input_path.resolve() == output_path.resolve():
            raise ValueError("input and output must be different files")
        pcd_fields = reader.fields
        columns = pcd_fields_to_columns(pcd_fields)
        source_count = reader.point_count
        origin = reader.origin

        downsampler = VoxelDownsampler(voxel_size)
        progress = tqdm(
            total=source_count, unit="pts", unit_scale=True,
            desc=input_path.name, disable=not show_progress,
        )
        with progress:
            for chunk in reader.chunks(chunk_size):
                downsampler.add(unpack_columns(chunk, pcd_fields))
                progress.update(len(chunk))

    points = downsampler.result()
    written = len(points)
    with BinaryPcdWriter(
        output_path, max_points=written, fields=pcd_fields, origin=origin
    ) as writer:
        writer.write_chunk(pack_columns(points, columns))

    return DownsampleResult(
        point_count=written,
        source_count=source_count,
        voxel_size=voxel_size,
        origin=origin,
        output_path=output_path,
        fields=pcd_fields,
    )
