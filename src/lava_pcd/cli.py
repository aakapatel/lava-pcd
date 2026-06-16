"""Command-line interface for lava-pcd."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from lava_pcd import __version__
from lava_pcd.convert import downsample_pcd, laz_to_pcd
from lava_pcd.crop import DEFAULT_MAX_DISPLAY, crop_pcd, select_rectangle
from lava_pcd.filtering import (
    DEFAULT_K,
    DEFAULT_MIN_NEIGHBORS,
    DEFAULT_RADIUS,
    DEFAULT_STD_RATIO,
    filter_pcd,
)
from lava_pcd.io.laz_reader import DEFAULT_CHUNK_SIZE

app = typer.Typer(
    help="Process large point cloud files (convert, downsample, voxelize, crop, merge).",
    add_completion=False,
)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"lava-pcd {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Optional[bool] = typer.Option(
        None, "--version", callback=_version_callback, is_eager=True,
        help="Show the version and exit.",
    ),
) -> None:
    """lava-pcd command-line interface."""


@app.command()
def convert(
    input: Path = typer.Argument(..., help="Input .laz/.las point cloud."),
    output: Path = typer.Argument(..., help="Output .pcd file."),
    chunk_size: int = typer.Option(
        DEFAULT_CHUNK_SIZE, "--chunk-size", "-c", min=1,
        help="Points read per chunk (lower = less memory).",
    ),
    voxel: float = typer.Option(
        0.0, "--voxel", "-v", min=0.0, prompt="Voxel size (0 = no downsampling)",
        help="Voxel-downsample resolution in coordinate units (0 disables).",
    ),
    fields: str = typer.Option(
        "auto", "--fields", "-f",
        help="Attributes to export: auto, xyz, intensity, rgb, all.",
    ),
    reproject: str = typer.Option(
        None, "--reproject", "-r", metavar="CRS",
        help="Reproject X/Y to this CRS, e.g. EPSG:32627 (UTM 27N).",
    ),
    origin: str = typer.Option(
        "header", "--origin", "-o",
        help="Local-origin shift: 'header' (LAS offset), 'none', or 'x,y,z'.",
    ),
    parallel: bool = typer.Option(
        False, "--parallel/--no-parallel",
        help="Use the multi-threaded LAZ backend (faster; panics on some files).",
    ),
    keep_invalid: bool = typer.Option(
        False, "--keep-invalid",
        help="Keep points outside the LAS header bounding box (off by default).",
    ),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress the progress bar."),
) -> None:
    """Convert a .laz/.las file to a binary .pcd (auto-selects RGB or intensity)."""
    origin_arg: str | tuple[float, float, float] = origin
    if origin not in ("header", "none"):
        try:
            parts = [float(p) for p in origin.replace(" ", "").split(",")]
        except ValueError:
            parts = []
        if len(parts) != 3:
            typer.secho(
                f"error: --origin must be 'header', 'none', or 'x,y,z', got {origin!r}",
                fg=typer.colors.RED, err=True,
            )
            raise typer.Exit(code=1)
        origin_arg = (parts[0], parts[1], parts[2])

    try:
        result = laz_to_pcd(
            input, output, chunk_size=chunk_size, origin=origin_arg,
            voxel_size=voxel, fields=fields, reproject=reproject,
            parallel=parallel, filter_bounds=not keep_invalid,
            show_progress=not quiet,
        )
    except (FileNotFoundError, ValueError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    ox, oy, oz = result.origin
    typer.secho(
        f"wrote {result.point_count:,} points -> {result.output_path}",
        fg=typer.colors.GREEN,
    )
    typer.echo(f"fields: {' '.join(result.fields)}")
    if result.reproject:
        typer.echo(f"reprojected to: {result.reproject}")
    if result.dropped:
        typer.echo(f"dropped {result.dropped:,} out-of-bounds points")
    if result.voxel_size > 0:
        ratio = result.source_count / result.point_count if result.point_count else 0.0
        typer.echo(
            f"voxel downsample @ {result.voxel_size}: "
            f"{result.source_count:,} -> {result.point_count:,} points "
            f"({ratio:.1f}x reduction)"
        )
    typer.echo(f"local origin (global = local + origin): {ox} {oy} {oz}")
    if result.sidecar_path is not None:
        typer.echo(f"origin metadata: {result.sidecar_path}")


@app.command()
def downsample(
    input: Path = typer.Argument(..., help="Input .pcd file."),
    output: Path = typer.Argument(..., help="Output (downsampled) .pcd file."),
    voxel: float = typer.Option(
        ..., "--voxel", "-v", min=0.0, prompt="Voxel size (coordinate units)",
        help="Voxel-downsample resolution in coordinate units (must be > 0).",
    ),
    chunk_size: int = typer.Option(
        DEFAULT_CHUNK_SIZE, "--chunk-size", "-c", min=1,
        help="Points read per chunk (lower = less memory).",
    ),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress the progress bar."),
) -> None:
    """Voxel-downsample an existing binary .pcd into a new .pcd."""
    try:
        result = downsample_pcd(
            input, output, voxel_size=voxel, chunk_size=chunk_size,
            show_progress=not quiet,
        )
    except (FileNotFoundError, ValueError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    ratio = result.source_count / result.point_count if result.point_count else 0.0
    typer.secho(
        f"wrote {result.point_count:,} points -> {result.output_path}",
        fg=typer.colors.GREEN,
    )
    typer.echo(f"fields: {' '.join(result.fields)}")
    typer.echo(
        f"voxel downsample @ {result.voxel_size}: "
        f"{result.source_count:,} -> {result.point_count:,} points "
        f"({ratio:.1f}x reduction)"
    )
    ox, oy, oz = result.origin
    typer.echo(f"local origin (global = local + origin): {ox} {oy} {oz}")


@app.command()
def crop(
    input: Path = typer.Argument(..., help="Input .pcd file."),
    output: Path = typer.Argument(..., help="Output (cropped) .pcd file."),
    bounds: str = typer.Option(
        None, "--bounds", "-b", metavar="MIN_A,MAX_A,MIN_B,MAX_B",
        help="Explicit rectangle over --axes (skips the interactive selector).",
    ),
    axes: str = typer.Option(
        "xy", "--axes", "-a",
        help="Axis pair the rectangle spans: xy (top-down), xz, or yz.",
    ),
    use_global: bool = typer.Option(
        False, "--global",
        help="Interpret --bounds in global coords (origin is subtracted).",
    ),
    max_display: int = typer.Option(
        DEFAULT_MAX_DISPLAY, "--max-display", min=1,
        help="Max points shown in the interactive selector (subsampled).",
    ),
    chunk_size: int = typer.Option(
        DEFAULT_CHUNK_SIZE, "--chunk-size", "-c", min=1,
        help="Points read per chunk (lower = less memory).",
    ),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress the progress bar."),
) -> None:
    """Crop a .pcd to a rectangle (interactive top-down selection by default)."""
    try:
        from lava_pcd.io.pcd_reader import BinaryPcdReader

        if bounds is not None:
            try:
                parts = [float(p) for p in bounds.replace(" ", "").split(",")]
            except ValueError:
                parts = []
            if len(parts) != 4:
                raise ValueError(
                    f"--bounds must be 'min_a,max_a,min_b,max_b', got {bounds!r}"
                )
            rect = (parts[0], parts[1], parts[2], parts[3])
            if use_global:
                with BinaryPcdReader(input) as r:
                    ox, oy, oz = r.origin
                o = {"x": ox, "y": oy, "z": oz}
                sa, sb = o[axes[0].lower()], o[axes[1].lower()]
                rect = (rect[0] - sa, rect[1] - sa, rect[2] - sb, rect[3] - sb)
        else:
            typer.echo("Opening interactive selector — drag a rectangle, then close the window.")
            rect = select_rectangle(input, axes=axes, max_display_points=max_display)

        result = crop_pcd(
            input, output, bounds=rect, axes=axes,
            chunk_size=chunk_size, show_progress=not quiet,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    lo_a, hi_a, lo_b, hi_b = result.bounds
    pct = 100.0 * result.point_count / result.source_count if result.source_count else 0.0
    typer.secho(
        f"wrote {result.point_count:,} points -> {result.output_path}",
        fg=typer.colors.GREEN,
    )
    typer.echo(f"fields: {' '.join(result.fields)}")
    typer.echo(
        f"crop {result.axes}: {result.axes[0]} [{lo_a:.3f}, {hi_a:.3f}]  "
        f"{result.axes[1]} [{lo_b:.3f}, {hi_b:.3f}]  (local coords)"
    )
    typer.echo(
        f"kept {result.point_count:,} of {result.source_count:,} points ({pct:.1f}%)"
    )
    ox, oy, oz = result.origin
    typer.echo(f"local origin (global = local + origin): {ox} {oy} {oz}")


@app.command(name="filter")
def filter_outliers(
    input: Path = typer.Argument(..., help="Input .pcd file."),
    output: Path = typer.Argument(..., help="Output (filtered) .pcd file."),
    method: str = typer.Option(
        "statistical", "--method", "-m",
        help="Outlier method: 'radius' or 'statistical'.",
    ),
    radius: float = typer.Option(
        DEFAULT_RADIUS, "--radius", "-r", min=0.0,
        help="[radius] neighbourhood radius in coordinate units.",
    ),
    min_neighbors: int = typer.Option(
        DEFAULT_MIN_NEIGHBORS, "--min-neighbors", "-n", min=1,
        help="[radius] min points (incl. self) within --radius to keep a point.",
    ),
    neighbors: int = typer.Option(
        DEFAULT_K, "--neighbors", "-k", min=1,
        help="[statistical] number of nearest neighbours for the mean distance.",
    ),
    std_ratio: float = typer.Option(
        DEFAULT_STD_RATIO, "--std-ratio", "-s", min=0.0,
        help="[statistical] keep points within mean + std_ratio*std of the mean distance.",
    ),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress the progress bar."),
) -> None:
    """Remove outliers from a .pcd (radius or statistical)."""
    try:
        result = filter_pcd(
            input, output, method=method, radius=radius,
            min_neighbors=min_neighbors, k=neighbors, std_ratio=std_ratio,
            show_progress=not quiet,
        )
    except (FileNotFoundError, ValueError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    pct = 100.0 * result.removed / result.source_count if result.source_count else 0.0
    typer.secho(
        f"wrote {result.point_count:,} points -> {result.output_path}",
        fg=typer.colors.GREEN,
    )
    typer.echo(f"fields: {' '.join(result.fields)}")
    params = "  ".join(f"{k}={v}" for k, v in result.params.items())
    typer.echo(f"method: {result.method} ({params})")
    typer.echo(
        f"removed {result.removed:,} of {result.source_count:,} points "
        f"({pct:.1f}%) as outliers"
    )
    ox, oy, oz = result.origin
    typer.echo(f"local origin (global = local + origin): {ox} {oy} {oz}")


if __name__ == "__main__":
    app()
