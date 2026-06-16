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

## Merge an aerial map with a lava-tube map (via skylights)

The headline feature: merge two clouds that barely overlap — an **aerial** surface
map and an underground **lava-tube** map. Their only shared geometry is the set of
**skylights** (collapse holes), which appear as voids in the aerial ground and as
holes in the tube ceiling. The skylight openings coincide in the world, so they form
a sparse landmark *constellation* common to both clouds; matching the two
constellations recovers the rigid transform even with almost no other overlap.

It is a **staged, inspectable** pipeline — detection is the fragile part, so you
confirm the landmarks before registering:

```bash
# 1. detect skylights in each cloud (saved as JSON landmark sets)
lava-pcd holes aerial.pcd aerial_holes.json --mode aerial --show
lava-pcd holes tube.pcd   tube_holes.json   --mode ceiling --up 0,0.2,0.98 --show

# 2. match the two constellations -> rigid transform (tube -> aerial)
lava-pcd register aerial_holes.json tube_holes.json -o transform.json

# 3. apply + write the merged cloud (optionally ICP-refined on the rims)
lava-pcd merge aerial.pcd tube.pcd merged.pcd -t transform.json --refine --source-field
```

### `holes` — detect (or place) skylights
```
#   -m / --mode          aerial (top-down voids) | ceiling (tube roof, along --up)
#        --up X,Y,Z       up-axis for ceiling/manual mode (estimated if omitted)
#   -r / --res            grid cell size in coordinate units (default 1.0)
#        --min-area       ignore voids smaller than this (coord units squared)
#        --max-area       ignore voids larger than this (optional)
#        --ceiling-jump   [ceiling] roof-height deviation flagged as an opening
#        --show           review/toggle detections interactively before saving
#        --manual         place skylights by hand (click each centre) — the fallback
```
`aerial` finds **enclosed empty regions** of a top-down occupancy grid (the laser
passes through a hole and gives no return). `ceiling` first rotates the tube so the
`--up` axis points up, then flags enclosed cells where the ceiling is missing or
jumps away from its neighbours. The tube is in an arbitrary SLAM frame, so give a
known `--up` when you have one; otherwise it is estimated (approximate). `--show`
lets you drop false positives; `--manual` lets you click skylights directly when
automatic detection struggles.

### `register` — match constellations
```
#   -o / --output     output transform .json (required)
#   -m / --mode       4dof (up-assisted, robust; default) | 6dof (Kabsch RANSAC)
#   -t / --tolerance  max landmark mismatch (coord units) to count as an inlier
```
`4dof` uses the up-axis from each hole set to reduce the match to yaw + translation —
solvable from **2 holes** and not degenerate for a straight (collinear) tube. `6dof`
is the fallback when no up-axis is trustworthy; it needs **≥3 non-collinear** matches.
The command prints the transform, the residual RMS, and warnings (too few / collinear
landmarks).

### `merge` — apply and combine
```
#   -t / --transform   transform .json from `register` (required)
#        --refine       ICP-refine on the matched rims before merging
#        --rim-radius   [--refine] radius around each skylight used for ICP
#   -s / --source-field add a 'source' channel (0=aerial, 1=tube) to colour by origin
#   -c / --chunk-size / -q / --quiet  as elsewhere
```
The tube is transformed into the **aerial** frame and the two clouds are written as a
single `.pcd`. They usually carry different fields (aerial RGB vs tube intensity), so
the merged cloud keeps `x y z` plus, with `--source-field`, a `source` channel to
colour by origin. `--refine` runs a small point-to-point ICP **only on the matched
skylight rims** (global overlap is too small for global ICP). There is also a thin
`lava-pcd transform IN OUT transform.json` to apply a transform to one cloud.

From Python:

```python
from lava_pcd import detect_holes, match_constellations, merge_clouds

aerial = detect_holes("aerial.pcd", mode="aerial")
tube = detect_holes("tube.pcd", mode="ceiling", up=(0, 0.2, 0.98))
tf = match_constellations(aerial, tube, mode="4dof")   # tube-local -> aerial-local
print(f"matched {len(tf.inliers)} skylights, RMS {tf.rms:.3f}")
merge_clouds("aerial.pcd", "tube.pcd", "merged.pcd", tf, refine=True, source_field=True)
```

**Caveats** (registration is hard when overlap is tiny): you need **≥2** matched
skylights for `4dof` (**≥3 non-collinear** for `6dof`); nearly collinear skylights
leave rotation about that line weakly constrained (mitigated by the up-axis and
`--refine`); and the tube up-axis is only estimated in an arbitrary frame, so prefer
`--up` or the `--manual` picker when results look off.

## Viewing

```bash
pcl_viewer output.pcd      # colours by the rgb/intensity field automatically
pcl_viewer merged.pcd      # press 2 to colour merged clouds by the 'source' field
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
