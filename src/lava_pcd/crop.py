"""Crop a binary ``.pcd`` to an axis-aligned rectangle.

Two ways to choose the rectangle:

* **Interactive** -- open a top-down (or side) scatter of the cloud and drag a
  rectangle with the mouse (matplotlib's ``RectangleSelector``). The cloud is
  randomly subsampled for display only; the crop is applied to every point.
* **Explicit bounds** -- pass ``(min_a, max_a, min_b, max_b)`` directly, for
  headless / scripted use.

The crop is a 2-D rectangle over a pair of axes (``xy`` top-down by default);
the third axis is kept in full, so an ``xy`` crop is a vertical "cookie-cutter"
column. Reading and writing are chunked, so memory stays bounded. The input's
fields and its local origin (the ``# LAVA_PCD_ORIGIN`` comment) are preserved.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from tqdm import tqdm

from lava_pcd.io.pcd_reader import BinaryPcdReader
from lava_pcd.io.pcd_writer import BinaryPcdWriter

DEFAULT_MAX_DISPLAY = 500_000

_AXES = {"x": 0, "y": 1, "z": 2}


@dataclass
class CropResult:
    """Outcome of a :func:`crop_pcd` call."""

    point_count: int
    source_count: int
    bounds: tuple[float, float, float, float]
    axes: str
    origin: tuple[float, float, float]
    output_path: Path
    fields: tuple[str, ...] = ()

    def __int__(self) -> int:
        return self.point_count


def _axis_indices(axes: str) -> tuple[int, int]:
    axes = axes.lower()
    if len(axes) != 2 or any(a not in _AXES for a in axes) or axes[0] == axes[1]:
        raise ValueError(f"axes must be two distinct axes from x/y/z, got {axes!r}")
    return _AXES[axes[0]], _AXES[axes[1]]


def _normalise_bounds(
    bounds: tuple[float, float, float, float]
) -> tuple[float, float, float, float]:
    a0, a1, b0, b1 = bounds
    lo_a, hi_a = (a0, a1) if a0 <= a1 else (a1, a0)
    lo_b, hi_b = (b0, b1) if b0 <= b1 else (b1, b0)
    return (lo_a, hi_a, lo_b, hi_b)


def crop_pcd(
    input_path: str | Path,
    output_path: str | Path,
    bounds: tuple[float, float, float, float],
    axes: str = "xy",
    chunk_size: int = 5_000_000,
    show_progress: bool = True,
) -> CropResult:
    """Crop ``input_path`` to ``bounds`` over ``axes`` and write ``output_path``.

    ``bounds`` is ``(min_a, max_a, min_b, max_b)`` in the *local* PCD frame
    (global = local + origin), where ``a`` / ``b`` are the two ``axes`` (e.g.
    ``"xy"`` -> a=x, b=y). Points whose two selected coordinates fall inside the
    rectangle (inclusive) are kept; the third axis is unconstrained.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    if input_path.suffix.lower() != ".pcd" or output_path.suffix.lower() != ".pcd":
        raise ValueError("crop works on .pcd input and output")

    ia, ib = _axis_indices(axes)
    lo_a, hi_a, lo_b, hi_b = _normalise_bounds(bounds)

    with BinaryPcdReader(input_path) as reader:
        if input_path.resolve() == output_path.resolve():
            raise ValueError("input and output must be different files")
        fields = reader.fields
        source_count = reader.point_count
        origin = reader.origin

        progress = tqdm(
            total=source_count, unit="pts", unit_scale=True,
            desc=input_path.name, disable=not show_progress,
        )
        with BinaryPcdWriter(
            output_path, max_points=source_count, fields=fields, origin=origin
        ) as writer, progress:
            for chunk in reader.chunks(chunk_size):
                a, b = chunk[:, ia], chunk[:, ib]
                mask = (a >= lo_a) & (a <= hi_a) & (b >= lo_b) & (b <= hi_b)
                kept = chunk[mask]
                if len(kept):
                    writer.write_chunk(np.ascontiguousarray(kept, dtype=np.float32))
                progress.update(len(chunk))
            written = writer._written

    return CropResult(
        point_count=written,
        source_count=source_count,
        bounds=(lo_a, hi_a, lo_b, hi_b),
        axes=axes.lower(),
        origin=origin,
        output_path=output_path,
        fields=fields,
    )


def _load_display_sample(
    input_path: Path, max_points: int, seed: int = 0
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Stream the cloud and keep a random subsample of at most ``max_points``."""
    rng = np.random.default_rng(seed)
    with BinaryPcdReader(input_path) as reader:
        total = max(reader.point_count, 1)
        keep_frac = min(1.0, max_points / total)
        fields = reader.fields
        parts: list[np.ndarray] = []
        for chunk in reader.chunks():
            if keep_frac < 1.0:
                chunk = chunk[rng.random(len(chunk)) < keep_frac]
            if len(chunk):
                parts.append(chunk)
    sample = np.concatenate(parts) if parts else np.empty((0, len(fields)), np.float32)
    return sample, fields


def _display_colours(sample: np.ndarray, fields: tuple[str, ...], ib: int) -> object:
    """Per-point colours for the scatter: unpacked RGB if present, else height."""
    if "rgb" in fields:
        packed = np.ascontiguousarray(sample[:, fields.index("rgb")], np.float32).view(np.uint32)
        r = ((packed >> 16) & 0xFF) / 255.0
        g = ((packed >> 8) & 0xFF) / 255.0
        b = (packed & 0xFF) / 255.0
        return np.column_stack([r, g, b])
    # No colour: shade by the out-of-plane axis (e.g. z for a top-down xy view).
    other = ({0, 1, 2} - {0, 1}) if ib in (0, 1) else {2}
    zcol = sample[:, next(iter(other))] if sample.shape[1] > 2 else sample[:, ib]
    return zcol


def select_rectangle(
    input_path: str | Path,
    axes: str = "xy",
    max_display_points: int = DEFAULT_MAX_DISPLAY,
) -> tuple[float, float, float, float]:
    """Open an interactive view and return the dragged rectangle's bounds.

    Returns ``(min_a, max_a, min_b, max_b)`` in the local PCD frame. Raises
    ``RuntimeError`` if the window is closed without a selection.
    """
    import matplotlib.pyplot as plt
    from matplotlib.widgets import RectangleSelector

    input_path = Path(input_path)
    ia, ib = _axis_indices(axes)
    sample, fields = _load_display_sample(input_path, max_display_points)
    if len(sample) == 0:
        raise RuntimeError(f"{input_path.name} has no points to display")

    a, b = sample[:, ia], sample[:, ib]
    colours = _display_colours(sample, fields, ib)

    selection: dict[str, tuple[float, float, float, float]] = {}

    def on_select(eclick, erelease) -> None:
        if None in (eclick.xdata, erelease.xdata, eclick.ydata, erelease.ydata):
            return
        selection["bounds"] = _normalise_bounds(
            (eclick.xdata, erelease.xdata, eclick.ydata, erelease.ydata)
        )

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(a, b, c=colours, s=0.5, marker=".", linewidths=0)
    ax.set_aspect("equal")
    ax.set_xlabel(f"{axes[0]} (local)")
    ax.set_ylabel(f"{axes[1]} (local)")
    ax.set_title(
        f"{input_path.name} — drag a rectangle, adjust handles, "
        "then close the window to crop"
    )
    selector = RectangleSelector(
        ax, on_select, useblit=True, interactive=True, button=[1],
        props=dict(facecolor="red", edgecolor="red", alpha=0.2, fill=True),
    )
    fig._lava_selector = selector  # keep a reference so it isn't GC'd
    plt.show()

    if "bounds" not in selection:
        raise RuntimeError("no rectangle was selected")
    return selection["bounds"]
