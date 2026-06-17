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
from lava_pcd.holes import (
    DEFAULT_CEILING_JUMP,
    DEFAULT_EDGE_MARGIN,
    DEFAULT_MERGE_FACTOR,
    DEFAULT_MIN_AREA,
    DEFAULT_MIN_DENSITY,
    DEFAULT_RELATIVE_WINDOW,
    DEFAULT_RESOLUTION,
    DEFAULT_SMOOTH,
    HoleSet,
    build_occupancy,
    detect_skylights,
    pick_holes,
    review_holes,
    show_occupancy,
)
from lava_pcd.io.laz_reader import DEFAULT_CHUNK_SIZE
from lava_pcd.merge import (
    DEFAULT_CMAP,
    DEFAULT_RIM_HEIGHT,
    DEFAULT_RIM_INFLATE,
    DEFAULT_RIM_RADIUS,
    apply_transform,
    merge_clouds,
)
from lava_pcd.register import (
    DEFAULT_TOLERANCE,
    Transform,
    match_constellations,
    vertical_residuals,
    visualize_match,
)

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


def _parse_up(up: str | None) -> tuple[float, float, float] | None:
    if up is None:
        return None
    try:
        parts = [float(p) for p in up.replace(" ", "").split(",")]
    except ValueError:
        parts = []
    if len(parts) != 3:
        raise typer.BadParameter(f"--up must be 'x,y,z', got {up!r}")
    return (parts[0], parts[1], parts[2])


@app.command()
def holes(
    input: Path = typer.Argument(..., help="Input .pcd file."),
    output: Path = typer.Argument(..., help="Output skylight .json file."),
    mode: str = typer.Option(
        "aerial", "--mode", "-m",
        help="'aerial' (top-down voids) or 'ceiling' (tube ceiling along --up).",
    ),
    up: str = typer.Option(
        None, "--up", metavar="X,Y,Z",
        help="Up-axis for 'ceiling'/'manual' (estimated if omitted).",
    ),
    res: float = typer.Option(
        DEFAULT_RESOLUTION, "--res", "-r", min=1e-6,
        help="Grid cell size in coordinate units.",
    ),
    min_density: float = typer.Option(
        DEFAULT_MIN_DENSITY, "--min-density", min=0.0,
        help="Cells with fewer points than this count as empty (raise to catch "
             "less-occupied holes). Ignored when --relative is set.",
    ),
    relative: float = typer.Option(
        None, "--relative", min=0.0,
        help="Adaptive threshold: a cell is empty below this fraction of its local "
             "median density (e.g. 0.3). Use when ground density varies across the map.",
    ),
    relative_window: int = typer.Option(
        DEFAULT_RELATIVE_WINDOW, "--relative-window", min=3,
        help="[--relative] window in cells for the local reference density.",
    ),
    smooth: float = typer.Option(
        DEFAULT_SMOOTH, "--smooth", min=0.0,
        help="Gaussian sigma (cells) to smooth the count image before thresholding.",
    ),
    min_area: float = typer.Option(
        DEFAULT_MIN_AREA, "--min-area", min=0.0,
        help="Ignore voids smaller than this area (coord units squared).",
    ),
    max_area: float = typer.Option(
        None, "--max-area", help="Ignore voids larger than this area (optional)."
    ),
    edge_margin: float = typer.Option(
        DEFAULT_EDGE_MARGIN, "--edge-margin", min=0.0,
        help="Reject holes within this distance (coord units) of the cloud boundary "
             "-- drops ragged-rim false positives on the tube ceiling.",
    ),
    merge_overlap: bool = typer.Option(
        True, "--merge-overlap/--no-merge-overlap",
        help="Fuse holes whose fitted ellipses overlap into one.",
    ),
    merge_factor: float = typer.Option(
        DEFAULT_MERGE_FACTOR, "--merge-factor", min=0.0,
        help="Scale on the ellipse radii for the overlap test (>1 merges near holes).",
    ),
    ceiling_jump: float = typer.Option(
        DEFAULT_CEILING_JUMP, "--ceiling-jump", min=0.0,
        help="[ceiling] ceiling-height deviation flagged as an opening.",
    ),
    show: bool = typer.Option(
        False, "--show/--no-show",
        help="Show the occupancy image with detections and let you toggle them.",
    ),
    manual: bool = typer.Option(
        False, "--manual", help="Place skylights by hand instead of detecting."
    ),
) -> None:
    """Detect (or place) skylight holes in a .pcd and save them as JSON."""
    try:
        up_vec = _parse_up(up)
        if manual:
            holeset = pick_holes(input, up=up_vec, resolution=res)
        else:
            grid = build_occupancy(input, mode=mode, up=up_vec, resolution=res)
            holes, hole_cells = detect_skylights(
                grid, min_area=min_area, max_area=max_area,
                min_density=min_density, smooth=smooth, relative=relative,
                relative_window=relative_window, edge_margin=edge_margin,
                merge_overlap=merge_overlap, merge_factor=merge_factor,
                ceiling_jump=ceiling_jump,
            )
            holeset = HoleSet(holes, grid.origin, grid.up, mode, res, str(input))
            if show:
                holeset = review_holes(grid, holeset, hole_cells)
        holeset.to_json(output)
    except (FileNotFoundError, ValueError, RuntimeError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    ux, uy, uz = holeset.up
    typer.secho(
        f"found {len(holeset.holes)} skylight(s) -> {output}", fg=typer.colors.GREEN
    )
    typer.echo(f"mode: {holeset.mode}   up: {ux:.3f} {uy:.3f} {uz:.3f}")
    for h in holeset.holes:
        cx, cy, cz = h.centroid
        deg = h.orientation * 180.0 / 3.141592653589793
        typer.echo(
            f"  #{h.id}: centre ({cx:.2f}, {cy:.2f}, {cz:.2f})  "
            f"ellipse a~{h.semi_major:.2f} b~{h.semi_minor:.2f} θ~{deg:.0f}°  "
            f"area {h.area:.1f}"
        )


@app.command()
def occupancy(
    input: Path = typer.Argument(..., help="Input .pcd file."),
    mode: str = typer.Option(
        "aerial", "--mode", "-m", help="'aerial' (top-down) or 'ceiling' (along --up)."
    ),
    up: str = typer.Option(
        None, "--up", metavar="X,Y,Z", help="Up-axis for 'ceiling' (estimated if omitted)."
    ),
    res: float = typer.Option(
        DEFAULT_RESOLUTION, "--res", "-r", min=1e-6, help="Grid cell size in coord units."
    ),
    min_density: float = typer.Option(
        None, "--min-density", min=0.0,
        help="If set, also preview the detected holes at this points-per-cell threshold.",
    ),
    relative: float = typer.Option(
        None, "--relative", min=0.0,
        help="Preview with the adaptive threshold (fraction of local median density).",
    ),
    relative_window: int = typer.Option(
        DEFAULT_RELATIVE_WINDOW, "--relative-window", min=3,
        help="[--relative] window in cells for the local reference density.",
    ),
    smooth: float = typer.Option(
        DEFAULT_SMOOTH, "--smooth", min=0.0,
        help="Gaussian sigma (cells) to smooth the count image before thresholding.",
    ),
    min_area: float = typer.Option(
        DEFAULT_MIN_AREA, "--min-area", min=0.0,
        help="[preview] ignore voids smaller than this area (coord units squared).",
    ),
    max_area: float = typer.Option(
        None, "--max-area", help="[preview] ignore voids larger than this area (optional)."
    ),
    edge_margin: float = typer.Option(
        DEFAULT_EDGE_MARGIN, "--edge-margin", min=0.0,
        help="[preview] reject holes within this distance of the cloud boundary.",
    ),
    merge_overlap: bool = typer.Option(
        True, "--merge-overlap/--no-merge-overlap",
        help="[preview] fuse holes whose fitted ellipses overlap.",
    ),
    merge_factor: float = typer.Option(
        DEFAULT_MERGE_FACTOR, "--merge-factor", min=0.0,
        help="[preview] scale on the ellipse radii for the overlap test.",
    ),
    linear: bool = typer.Option(
        False, "--linear", help="Use a linear colour scale (default is log)."
    ),
) -> None:
    """Show the 2-D occupancy histogram of a cloud (to pick --res / --min-density)."""
    try:
        up_vec = _parse_up(up)
        grid = build_occupancy(input, mode=mode, up=up_vec, resolution=res)
        hole_cells = None
        if min_density is not None or relative is not None:
            # Preview matches `holes`: the overlay is the area-filtered hole mask.
            _, hole_cells = detect_skylights(
                grid, min_area=min_area, max_area=max_area,
                min_density=min_density if min_density is not None else DEFAULT_MIN_DENSITY,
                smooth=smooth, relative=relative, relative_window=relative_window,
                edge_margin=edge_margin, merge_overlap=merge_overlap,
                merge_factor=merge_factor,
            )
    except (FileNotFoundError, ValueError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    nx, ny = grid.counts.shape
    typer.echo(
        f"{input.name}: {nx}x{ny} cells @ res {res}  "
        f"(occupied {(grid.counts > 0).sum():,} / {nx * ny:,})"
    )
    show_occupancy(grid, hole_cells=hole_cells, log=not linear,
                   title=f"{input.name} occupancy ({mode}, res {res:g})")


@app.command()
def register(
    aerial_holes: Path = typer.Argument(..., help="Aerial skylight .json."),
    tube_holes: Path = typer.Argument(..., help="Lava-tube skylight .json."),
    output: Path = typer.Option(
        ..., "--output", "-o", help="Output transform .json (tube -> aerial)."
    ),
    mode: str = typer.Option(
        "auto", "--mode", "-m",
        help="'auto' (4dof then 6dof; default), '4dof' (up-assisted) or '6dof' (Kabsch).",
    ),
    tolerance: float = typer.Option(
        DEFAULT_TOLERANCE, "--tolerance", "-t", min=0.0,
        help="Max landmark mismatch (coord units) to count as an inlier.",
    ),
    min_inliers: int = typer.Option(
        3, "--min-inliers", "-n", min=2,
        help="Min mutually-consistent skylights required (the outlier-rejection lever).",
    ),
    show: bool = typer.Option(
        False, "--show", help="Plot the aligned constellations (matches + outliers)."
    ),
) -> None:
    """Match two skylight constellations into a rigid transform (tube -> aerial)."""
    try:
        aerial = HoleSet.from_json(aerial_holes)
        tube = HoleSet.from_json(tube_holes)
        transform = match_constellations(
            aerial, tube, mode=mode, tolerance=tolerance, min_inliers=min_inliers
        )
        transform.to_json(output)
    except (FileNotFoundError, ValueError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    typer.secho(
        f"matched {transform.n_inliers} skylight(s), RMS {transform.rms:.3f} "
        f"-> {output}", fg=typer.colors.GREEN,
    )
    typer.echo(
        f"mode: {transform.mode}   uniqueness margin: {transform.margin}   "
        f"shape score: {transform.shape_score:.3f}"
    )
    typer.echo(f"matches (tube -> aerial hole id): {transform.inliers}")
    vres = vertical_residuals(aerial, tube, transform)
    typer.echo(
        "per-skylight vertical residual (m): "
        + " ".join(f"{v:+.2f}" for v in vres)
        + "   (large spread => skylights disagree on depth; tune `merge --z-offset`)"
    )
    typer.echo("transform (tube-local -> aerial-local):")
    for row in transform.array:
        typer.echo("  " + "  ".join(f"{v: .4f}" for v in row))
    for w in transform.warnings:
        typer.secho(f"warning: {w}", fg=typer.colors.YELLOW)
    if show:
        visualize_match(aerial, tube, transform)


@app.command()
def merge(
    aerial: Path = typer.Argument(..., help="Aerial .pcd (target frame)."),
    tube: Path = typer.Argument(..., help="Lava-tube .pcd (to be transformed)."),
    output: Path = typer.Argument(..., help="Output merged .pcd."),
    transform: Path = typer.Option(
        ..., "--transform", "-t", help="Transform .json from `register`."
    ),
    refine: bool = typer.Option(
        False, "--refine", help="GICP-refine on the matched rims before merging."
    ),
    dof: int = typer.Option(
        4, "--dof",
        help="[--refine] GICP degrees of freedom: 4 (yaw + translation, keeps the "
             "tube from tilting/flattening) or 6 (full rigid).",
    ),
    color: bool = typer.Option(
        True, "--color/--no-color",
        help="Keep aerial RGB and elevation-colour the tube (rgb field).",
    ),
    cmap: str = typer.Option(
        DEFAULT_CMAP, "--cmap",
        help="Matplotlib colormap for the tube's elevation colour.",
    ),
    z_offset: float = typer.Option(
        0.0, "--z-offset", "-z",
        help="Slide the tube vertically (along aerial up) to set its roof depth.",
    ),
    source_field: bool = typer.Option(
        False, "--source-field", "-s",
        help="Also add a 'source' channel (0=aerial, 1=tube).",
    ),
    rim_inflate: float = typer.Option(
        DEFAULT_RIM_INFLATE, "--rim-inflate", min=1.0,
        help="[--refine] grow each skylight's ellipse by this factor to gather its "
             "rim points (so the patch fits each hole's size/shape automatically).",
    ),
    rim_radius: float = typer.Option(
        DEFAULT_RIM_RADIUS, "--rim-radius", min=0.0,
        help="[--refine] fallback horizontal radius for holes that have no ellipse "
             "data; also the GICP correspondence-clamp scale.",
    ),
    rim_height: float = typer.Option(
        DEFAULT_RIM_HEIGHT, "--rim-height", min=0.0,
        help="[--refine] vertical band around each opening for ICP (excludes the "
             "deep tube body / far ground that would collapse the fit).",
    ),
    show_rims: bool = typer.Option(
        False, "--show-rims",
        help="Plot the rim points the ICP operates on (aerial vs tube, before/after).",
    ),
    chunk_size: int = typer.Option(
        DEFAULT_CHUNK_SIZE, "--chunk-size", "-c", min=1,
        help="Points read per chunk (lower = less memory).",
    ),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress the progress bar."),
) -> None:
    """Apply a transform to the tube cloud and write a merged .pcd."""
    try:
        tf = Transform.from_json(transform)
        result = merge_clouds(
            aerial, tube, output, tf, refine=refine, color=color,
            elevation_cmap=cmap, z_offset=z_offset, source_field=source_field,
            rim_radius=rim_radius, rim_height=rim_height, rim_inflate=rim_inflate,
            dof=dof, show_rims=show_rims, chunk_size=chunk_size, show_progress=not quiet,
        )
    except (FileNotFoundError, ValueError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    typer.secho(
        f"wrote {result.point_count:,} points -> {result.output_path}",
        fg=typer.colors.GREEN,
    )
    typer.echo(f"fields: {' '.join(result.fields)}")
    typer.echo(
        f"aerial {result.aerial_count:,} + tube {result.tube_count:,} points"
        + (f"   (ICP-refined, rim RMS {result.rms:.3f})" if result.refined else "")
    )
    ox, oy, oz = result.origin
    typer.echo(f"local origin (global = local + origin): {ox} {oy} {oz}")


@app.command(name="transform")
def transform_cmd(
    input: Path = typer.Argument(..., help="Input .pcd file."),
    output: Path = typer.Argument(..., help="Output (transformed) .pcd file."),
    transform_json: Path = typer.Argument(..., help="Transform .json from `register`."),
    chunk_size: int = typer.Option(
        DEFAULT_CHUNK_SIZE, "--chunk-size", "-c", min=1,
        help="Points read per chunk (lower = less memory).",
    ),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress the progress bar."),
) -> None:
    """Apply a registration transform to a single .pcd (preserves fields)."""
    try:
        tf = Transform.from_json(transform_json)
        apply_transform(input, output, tf.array, chunk_size=chunk_size,
                        show_progress=not quiet)
    except (FileNotFoundError, ValueError) as err:
        typer.secho(f"error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    typer.secho(f"wrote -> {output}", fg=typer.colors.GREEN)


if __name__ == "__main__":
    app()
