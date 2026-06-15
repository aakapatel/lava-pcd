"""Command-line interface for lava-pcd."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from lava_pcd import __version__
from lava_pcd.convert import laz_to_pcd
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
    origin: str = typer.Option(
        "header", "--origin", "-o",
        help="Local-origin shift: 'header' (LAS offset), 'none', or 'x,y,z'.",
    ),
    parallel: bool = typer.Option(
        False, "--parallel/--no-parallel",
        help="Use the multi-threaded LAZ backend (faster; panics on some files).",
    ),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress the progress bar."),
) -> None:
    """Convert a .laz/.las file to a binary .pcd (x y z intensity)."""
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
            parallel=parallel, show_progress=not quiet,
        )
    except (FileNotFoundError, ValueError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    ox, oy, oz = result.origin
    typer.secho(
        f"wrote {result.point_count:,} points -> {result.output_path}",
        fg=typer.colors.GREEN,
    )
    typer.echo(f"local origin (global = local + origin): {ox} {oy} {oz}")
    if result.sidecar_path is not None:
        typer.echo(f"origin metadata: {result.sidecar_path}")


if __name__ == "__main__":
    app()
