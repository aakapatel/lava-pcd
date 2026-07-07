"""High-level conversion routines."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import laspy
import numpy as np
from tqdm import tqdm

from lava_pcd.io.laz_reader import DEFAULT_CHUNK_SIZE, LazChunkReader
from lava_pcd.io.pcd_reader import BinaryPcdReader, pcd_fields_to_columns, unpack_columns
from lava_pcd.io.pcd_writer import BinaryPcdWriter, pack_columns, pcd_fields_for
from lava_pcd.voxel import VoxelDownsampler

_LAZ_SUFFIXES = {".laz", ".las"}
# Output formats for pcd_to_las: .las (uncompressed) and .laz (laszip/lazrs).
_LAS_OUT_SUFFIXES = {".las", ".laz"}
# PCD fields that map onto native LAS dimensions; anything else becomes an
# ``ExtraBytes`` dimension (e.g. the ``source`` channel that ``merge`` adds).
_STANDARD_PCD_FIELDS = {"x", "y", "z", "intensity", "rgb"}


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


@dataclass
class ExportResult:
    """Outcome of a :func:`pcd_to_las` call."""

    point_count: int
    output_path: Path
    origin: tuple[float, float, float]
    fields: tuple[str, ...] = ()
    scale: float = 0.0
    crs: str | None = None
    extra_dims: tuple[str, ...] = ()

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


def pcd_to_las(
    input_path: str | Path,
    output_path: str | Path,
    crs: str | None = None,
    scale: float = 0.001,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    show_progress: bool = True,
) -> ExportResult:
    """Convert a binary ``.pcd`` back into a ``.las`` / ``.laz`` file.

    The reverse of :func:`laz_to_pcd`. The output suffix selects the encoding:
    ``.las`` is uncompressed, ``.laz`` is laszip/lazrs-compressed. Reads the
    input in chunks of ``chunk_size`` points so memory stays bounded.

    Coordinates are written **georeferenced**: the cloud's local-origin shift
    (the ``# LAVA_PCD_ORIGIN`` header comment, i.e. ``global = local + origin``)
    is added back, and that origin is stored as the LAS header offset so the
    quantised integer coordinates stay small and exact.

    ``scale`` is the LAS coordinate quantisation step in coordinate units
    (default ``0.001`` = 1 mm); coordinates are stored as
    ``round((global - offset) / scale)`` in int32.

    Field mapping:

    * ``intensity`` -> LAS ``intensity`` (rounded/clamped to uint16).
    * packed ``rgb`` -> LAS ``red``/``green``/``blue`` (8-bit channels scaled to
      16-bit by ``x257``; point format 2). Without colour, point format 0.
    * any other field (e.g. ``merge``'s ``source`` channel) -> a float32
      ``ExtraBytes`` dimension of the same name.

    ``crs`` optionally embeds a coordinate reference system in the header
    (e.g. ``"EPSG:32627"``); ``.pcd`` files carry no CRS of their own, so pass
    the one the cloud is in (the output CRS used during the original convert).

    Returns an :class:`ExportResult`.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    if input_path.suffix.lower() != ".pcd":
        raise ValueError(
            f"expected a .pcd input, got '{input_path.suffix}' ({input_path})"
        )
    if output_path.suffix.lower() not in _LAS_OUT_SUFFIXES:
        raise ValueError(
            f"expected a .las or .laz output, got '{output_path.suffix}' ({output_path})"
        )
    if scale <= 0:
        raise ValueError(f"scale must be > 0, got {scale}")

    with BinaryPcdReader(input_path) as reader:
        fields = reader.fields
        for axis in ("x", "y", "z"):
            if axis not in fields:
                raise ValueError(
                    f"{input_path.name}: PCD is missing the required '{axis}' field"
                )
        has_intensity = "intensity" in fields
        has_rgb = "rgb" in fields
        extra = tuple(f for f in fields if f not in _STANDARD_PCD_FIELDS)
        origin = reader.origin
        ox, oy, oz = origin
        point_count = reader.point_count

        header = laspy.LasHeader(point_format=2 if has_rgb else 0, version="1.4")
        header.scales = [scale, scale, scale]
        header.offsets = [ox, oy, oz]
        for name in extra:
            header.add_extra_dim(laspy.ExtraBytesParams(name=name, type=np.float32))
        if crs is not None:
            import pyproj  # lazy: only needed when embedding a CRS

            header.add_crs(pyproj.CRS.from_user_input(crs))

        # unpack_columns expands a packed rgb field to r/g/b columns; map names.
        internal = pcd_fields_to_columns(fields)
        cidx = {name: i for i, name in enumerate(internal)}

        output_path.parent.mkdir(parents=True, exist_ok=True)
        written = 0
        progress = tqdm(
            total=point_count, unit="pts", unit_scale=True,
            desc=input_path.name, disable=not show_progress,
        )
        with laspy.open(output_path, mode="w", header=header) as writer, progress:
            for chunk in reader.chunks(chunk_size):
                cols = unpack_columns(chunk, fields)
                k = len(cols)
                record = laspy.ScaleAwarePointRecord.zeros(k, header=header)
                # Restore global coords in float64; laspy quantises with offset/scale.
                record.x = cols[:, cidx["x"]].astype(np.float64) + ox
                record.y = cols[:, cidx["y"]].astype(np.float64) + oy
                record.z = cols[:, cidx["z"]].astype(np.float64) + oz
                if has_intensity:
                    record.intensity = (
                        np.clip(np.round(cols[:, cidx["intensity"]]), 0, 65535)
                        .astype(np.uint16)
                    )
                if has_rgb:
                    for channel, dim in (("r", "red"), ("g", "green"), ("b", "blue")):
                        setattr(
                            record, dim,
                            np.clip(np.round(cols[:, cidx[channel]] * 257.0), 0, 65535)
                            .astype(np.uint16),
                        )
                for name in extra:
                    setattr(record, name, cols[:, cidx[name]].astype(np.float32))
                writer.write_points(record)
                written += k
                progress.update(k)

    return ExportResult(
        point_count=written,
        output_path=output_path,
        origin=origin,
        fields=fields,
        scale=scale,
        crs=crs,
        extra_dims=extra,
    )
