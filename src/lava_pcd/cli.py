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
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress the progress bar."),
) -> None:
    """Convert a .laz/.las file to a binary .pcd (x y z intensity)."""
    try:
        count = laz_to_pcd(
            input, output, chunk_size=chunk_size, show_progress=not quiet
        )
    except (FileNotFoundError, ValueError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    typer.secho(
        f"wrote {count:,} points -> {output}", fg=typer.colors.GREEN
    )


if __name__ == "__main__":
    app()
