"""High-level conversion routines."""

from __future__ import annotations

from pathlib import Path

from tqdm import tqdm

from lava_pcd.io.laz_reader import DEFAULT_CHUNK_SIZE, LazChunkReader
from lava_pcd.io.pcd_writer import BinaryPcdWriter

_LAZ_SUFFIXES = {".laz", ".las"}


def laz_to_pcd(
    input_path: str | Path,
    output_path: str | Path,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    show_progress: bool = True,
) -> int:
    """Convert a ``.laz`` / ``.las`` file to a binary ``.pcd`` (x y z intensity).

    Reads the input in chunks of ``chunk_size`` points to keep memory bounded,
    streaming each chunk into the output PCD. Returns the number of points
    written.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.exists():
        raise FileNotFoundError(f"input file not found: {input_path}")
    if input_path.suffix.lower() not in _LAZ_SUFFIXES:
        raise ValueError(
            f"expected a .laz or .las input, got '{input_path.suffix}' ({input_path})"
        )

    with LazChunkReader(input_path, chunk_size=chunk_size) as reader:
        total = reader.point_count
        with BinaryPcdWriter(output_path, num_points=total) as writer:
            progress = tqdm(
                total=total,
                unit="pts",
                unit_scale=True,
                desc=input_path.name,
                disable=not show_progress,
            )
            with progress:
                for chunk in reader.chunks():
                    writer.write_chunk(chunk)
                    progress.update(len(chunk))

    return total
