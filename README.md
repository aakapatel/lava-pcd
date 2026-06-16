# lava-pcd

Tools for processing large point cloud files — convert, downsample, voxelize, crop, merge.

The first capability is converting `.laz`/`.las` LiDAR clouds into **binary `.pcd`** files
(XYZ plus RGB or intensity) that work both in Python and with `pcl_viewer`, with optional
voxel downsampling and CRS reprojection.

## Install

```bash
pip install -e .
```

This pulls `laspy[lazrs,laszip]` (both LAZ backends), `pyproj` (reprojection), `numpy`,
`typer`, and `tqdm`.

## Convert a file

```bash
lava-pcd convert input.laz output.pcd
# options:
#   -f / --fields      attributes to export: auto, xyz, intensity, rgb, all (default auto)
#   -v / --voxel       voxel-downsample resolution in coord units (0 = off).
#                      Prompted interactively if not given.
#   -r / --reproject   reproject X/Y to a CRS, e.g. EPSG:32627 (UTM 27N)
#   -c / --chunk-size  points read per chunk (lower = less memory; default 5,000,000)
#   -o / --origin      local-origin shift: 'header' (default), 'none', or 'x,y,z'
#       --parallel     use the multi-threaded LAZ backend (faster; panics on some files)
#       --keep-invalid keep points outside the LAS header bbox (off by default)
#   -q / --quiet       suppress the progress bar
```

`--fields auto` exports **RGB** when the file is colourised (e.g. photogrammetry /
aerial products, where LiDAR intensity is often empty), otherwise **intensity**. RGB is
written as PCL's packed `rgb` float field, which `pcl_viewer` and Open3D colour by
automatically. Use `intensity`, `rgb`, `xyz`, or `all` to force a choice.

`--voxel`/`-v` downsamples to one centroid per cubic voxel of the given edge length
(in coordinate units, e.g. metres). `0` disables it. Downsampling is streamed, so it
stays memory-bounded even for hundred-million-point clouds.

`--reproject`/`-r` transforms X/Y from the file's CRS into the given one (Z unchanged),
using `pyproj`. Useful when the source is geographic (lon/lat degrees) and you need
metres — e.g. `-r EPSG:32627` for UTM zone 27N.

From Python:

```python
from lava_pcd import laz_to_pcd

res = laz_to_pcd(
    "input.laz", "output.pcd", voxel_size=0.1, reproject="EPSG:32627"
)
print(f"{res.source_count} -> {res.point_count} points, origin {res.origin}")
```

## Downsample an existing .pcd

Voxel-downsample a `.pcd` that was already written (e.g. to make a lighter copy
for viewing) without re-reading the source LAZ:

```bash
lava-pcd downsample input.pcd output.pcd -v 0.5
# options:
#   -v / --voxel       voxel resolution in coord units (prompted if not given)
#   -c / --chunk-size  points read per chunk (lower = less memory; default 5,000,000)
#   -q / --quiet       suppress the progress bar
```

Each output point is the centroid of its cubic voxel; `intensity` and RGB are
averaged the same way. Reading is chunked, so memory stays bounded by the
*downsampled* size. The input's fields and its local origin (the
`# LAVA_PCD_ORIGIN` header comment) are preserved. Only the binary float32 PCD
this package writes is accepted as input.

```python
from lava_pcd import downsample_pcd

res = downsample_pcd("input.pcd", "output.pcd", voxel_size=0.5)
print(f"{res.source_count} -> {res.point_count} points")
```

## Crop a .pcd to a rectangle

Cut out a rectangular region of a `.pcd`. By default it opens an interactive
top-down view — drag a rectangle with the mouse, adjust the handles, then close
the window and the crop is applied to the full cloud:

```bash
lava-pcd crop input.pcd output.pcd
# options:
#   -a / --axes        axis pair the rectangle spans: xy (top-down), xz, yz (default xy)
#   -b / --bounds      MIN_A,MAX_A,MIN_B,MAX_B — skip the GUI, crop headless
#       --global       interpret --bounds in global coords (origin is subtracted)
#       --max-display  max points drawn in the selector (subsampled; default 500,000)
#   -c / --chunk-size  points read per chunk (lower = less memory)
#   -q / --quiet       suppress the progress bar
```

The rectangle spans two axes (top-down `xy` by default); the third axis is kept
in full, so an `xy` crop is a vertical "cookie-cutter" column. The selector
subsamples the cloud for display only — every point is tested against the
rectangle when writing. For scripted/headless use, give `--bounds` directly:

```bash
# local coords (as shown in the selector / PCD header):
lava-pcd crop input.pcd out.pcd -b 120,260,80,210
# or global coords (e.g. UTM), with the origin subtracted automatically:
lava-pcd crop input.pcd out.pcd -b 500120,500260,7000080,7000210 --global
```

From Python:

```python
from lava_pcd import crop_pcd, select_rectangle

rect = select_rectangle("input.pcd", axes="xy")     # interactive; returns bounds
res = crop_pcd("input.pcd", "out.pcd", bounds=rect)  # or pass explicit bounds
print(f"kept {res.point_count} of {res.source_count} points")
```

## Filter outliers

Remove stray/noise points from a `.pcd` with one of two methods:

```bash
# statistical: drop points whose mean distance to their k nearest neighbours
# is an outlier (mean + std_ratio*std over the cloud). Defaults k=20, std=2.0.
lava-pcd filter input.pcd output.pcd -m statistical -k 20 -s 2.0

# radius: drop points with fewer than --min-neighbors points (incl. self)
# within --radius (coordinate units).
lava-pcd filter input.pcd output.pcd -m radius -r 0.5 -n 5
```

```
# options:
#   -m / --method         radius | statistical (default statistical)
#   -r / --radius         [radius] neighbourhood radius in coord units
#   -n / --min-neighbors  [radius] min points (incl. self) within radius to keep
#   -k / --neighbors      [statistical] number of nearest neighbours
#   -s / --std-ratio      [statistical] keep within mean + std_ratio*std
#   -q / --quiet          suppress the progress bar
```

These match the usual PCL / Open3D semantics. Both build a KD-tree over the XYZ
of the **whole** cloud, so the cloud is loaded into memory (neighbour queries
can't be streamed) — filter large clouds *after* cropping/downsampling. All
fields and the local origin are preserved on the survivors.

```python
from lava_pcd import filter_pcd

res = filter_pcd("input.pcd", "out.pcd", method="statistical", k=20, std_ratio=2.0)
print(f"removed {res.removed} of {res.source_count} points")
```

## Viewing

```bash
pcl_viewer output.pcd      # colours by the rgb/intensity field automatically
```

## Notes

- **Reading is chunked**, so memory stays bounded for very large clouds.
- **Local origin:** coordinates are stored as **float32**, so a local origin (the LAS
  header offset, or the reprojected bbox corner) is subtracted first to keep survey-grade
  (e.g. UTM) values precise. The origin is recorded in the PCD header comment and a
  `<output>.origin.json` sidecar — recover global coords via `global = local + origin`.
- **Invalid points:** points outside the LAS header bounding box (e.g. INT32 sentinels, or
  garbage from a partially-corrupt LAZ) are dropped by default; use `--keep-invalid` to keep
  them. A small number of stray points at extreme coordinates will otherwise wreck the
  cloud's scale in a viewer.
- **Corrupt LAZ:** if a viewer shows the cloud collapsed to a speck at huge axis values, the
  source LAZ may be corrupt. Both `lazrs` and the reference `laszip` backend failing at the
  same point confirms file corruption rather than a decoder bug.

## Development

```bash
pip install -e ".[dev]"
pytest -q
```
