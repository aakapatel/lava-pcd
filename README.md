# lava-pcd

Analysis code for the manuscript "Autonomous interior survey of a lava tube
for planetary habitat assessment" (Patel, Stathoulopoulos, Llewellin,
Oskarsson, Nikolakopoulos; under review, 2026), plus the point-cloud tools
it grew out of. Project page: https://aakapatel.github.io/lavatube-survey/.
Data release (maps, surface model, derived products):
https://github.com/aakapatel/lavatube-survey/releases/tag/data-v1.

## Reproducing the paper

1. `pip install -e .` in a Python 3.12 environment (the `.venv` used for the
   paper is not committed).
2. Download the data release into `maps/` (the `*_ed_ortho.pcd` interior
   maps, the `aerial_isn16_ortho_*.pcd` surface products and the two DEM
   crops) and unzip `analysis_outputs_baseline_and_ed.zip` into the repo
   root (it restores `analysis_out/` and `analysis_out_ed/`, which hold the
   registration transforms and mission statistics the pipeline reads).
3. `bash analysis/run_v9_pipeline.sh` recomputes every number into
   `analysis_out_v9/` (`paper_numbers.json` is the collection point), and
   `bash analysis/make_v9_figures.sh` plus `analysis/make_nature_figures.py`,
   `analysis/make_nature_ed_figures.py` and
   `analysis/make_nature_ed4_deployments.py` regenerate the figures into
   `analysis_out_v10/figures_nature/`.
4. `analysis/behaviour_tree.py` draws the behaviour tree of Extended Data
   Fig. 2; `analysis/zenodo_deposit.py` creates the Zenodo record.

The committed `analysis_out_v9/` is the paper run, so the numbers in the
manuscript can be checked without rerunning anything.

## Point-cloud tools

Tools for processing large point cloud files — convert, downsample, voxelize, crop, merge.

The first capability is converting `.laz`/`.las` LiDAR clouds into **binary `.pcd`** files
(XYZ plus RGB or intensity) that work both in Python and with `pcl_viewer`, with optional
voxel downsampling and CRS reprojection.

## Install

```bash
pip install -e .
```

This pulls `laspy[lazrs,laszip]` (both LAZ backends), `pyproj` (reprojection), `numpy`,
`typer`, `tqdm`, `scipy` (KD-trees / grid filters), `matplotlib` (the interactive
selectors and plots), and `small_gicp` (the `merge --refine` GICP).

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

## Convert a .pcd back to .las/.laz

The reverse of `convert`: write a binary `.pcd` back out as a LAS/LAZ file. The
output suffix selects the encoding — `.las` is uncompressed, `.laz` is
laszip-compressed:

```bash
lava-pcd to-las input.pcd output.las
# options:
#   --crs              embed a CRS in the header, e.g. EPSG:32627 (.pcd carries none)
#   -s / --scale       LAS coordinate quantisation step in coord units (default 0.001 = 1 mm)
#   -c / --chunk-size  points read per chunk (lower = less memory; default 5,000,000)
#   -q / --quiet       suppress the progress bar
```

Coordinates are written **georeferenced**: the cloud's local-origin shift (the
`# LAVA_PCD_ORIGIN` header comment, `global = local + origin`) is added back and
stored as the LAS header offset, so the quantised integer coordinates stay small
and exact. `intensity` maps to LAS `intensity`, packed `rgb` to 16-bit
`red`/`green`/`blue` (point format 2), and any other field (e.g. `merge`'s
`source` channel) becomes a float32 `ExtraBytes` dimension. A `.pcd` stores no
CRS of its own, so pass `--crs` to label the output.

```python
from lava_pcd import pcd_to_las

res = pcd_to_las("input.pcd", "output.las", crs="EPSG:32627")
print(f"wrote {res.point_count} points, offset {res.origin}")
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

# 3. apply + write the merged cloud (optionally GICP-refined on the rims)
lava-pcd merge aerial.pcd tube.pcd merged.pcd -t transform.json --refine --source-field
```

### Detection: a 2-D occupancy image

Both modes work by projecting the cloud to a **2-D occupancy histogram** (points per
cell, looking along the up-axis) and finding **enclosed low-density regions** — that's
where a skylight is. Inspect that image first to pick good parameters:

```bash
lava-pcd occupancy aerial.pcd --res 1.0                 # view the count image (log scale)
lava-pcd occupancy aerial.pcd --res 1.0 --min-density 3 # overlay the holes you'd detect
```

With `--min-density`, the red overlay is exactly what `holes` would detect — it honours
the same `--min-density`, `--smooth`, `--min-area` and `--max-area`, so you can dial all
of them in here before running `holes`.

A skylight rarely reads as *perfectly empty* — a few stray returns land inside it. So a
cell counts as "ground" only when it holds at least `--min-density` points (raise it to
turn *less-occupied* holes into voids), optionally after `--smooth`-ing the count image
to wash out isolated stray points. Tune `--res`, `--min-density`, `--smooth` against the
`occupancy` view, then run `holes` with the same values.

`--min-density` is an *absolute* points-per-cell threshold, so it assumes the ground
density is roughly uniform. If density varies across the map (flight-line overlap, range,
incidence angle), use **`--relative FRAC`** instead: a cell is empty when its count drops
below `FRAC` of its *local* median density (over a `--relative-window` cells window, default
21). This catches a hole in a sparse area without over-flagging dense areas, and it
overrides `--min-density`. Preview it the same way: `lava-pcd occupancy ... --relative 0.3`.

### `holes` — detect (or place) skylights
```
#   -m / --mode          aerial (top-down voids) | ceiling (tube roof, along --up)
#        --up X,Y,Z       up-axis for ceiling/manual mode (estimated if omitted)
#   -r / --res            grid cell size in coordinate units (default 1.0)
#        --min-density    points/cell below which a cell is "empty" (default 1)
#        --relative       FRAC: empty below this fraction of the local median density
#        --relative-window  window (cells) for the local reference density (default 21)
#        --smooth         Gaussian sigma (cells) to smooth counts before thresholding
#        --min-area       ignore voids smaller than this (coord units squared)
#        --max-area       ignore voids larger than this (optional)
#        --edge-margin    reject holes within this distance of the cloud boundary
#        --merge-overlap / --no-merge-overlap  fuse holes whose ellipses overlap (on by default)
#        --merge-factor   scale on the ellipse radii for the overlap test (>1 merges near holes)
#        --ceiling-jump   [ceiling] roof-height deviation flagged as an opening
#        --show           show the occupancy image + detections; click to toggle holes
#        --manual         place skylights by hand (click each centre) — the fallback
```
`aerial` finds **enclosed low-density regions** of a top-down occupancy grid (the laser
passes through a hole and gives few/no returns). `ceiling` first rotates the tube so the
`--up` axis points up, then flags enclosed cells where the ceiling is missing or jumps
away from its neighbours. The tube is in an arbitrary SLAM frame, so give a known `--up`
when you have one; otherwise it is estimated (approximate). `--show` overlays the
detections on the occupancy image so you can drop false positives; `--manual` lets you
click skylights directly when automatic detection struggles.

Each saved skylight carries an **equivalent ellipse** fitted from its cell second
moments — `semi_major`, `semi_minor` and `orientation` (major-axis bearing in the
up-plane) — so you get its rough shape and direction, not just a radius. `--show` draws
those ellipses; the JSON stores them for later use in matching.

When a skylight gets split into pieces (a thin bridge of returns across it, a rim nick
touching a real hole), the detector returns several overlapping ellipses. By default
(**`--merge-overlap`**) any holes whose fitted ellipses overlap are fused into one — their
cells are pooled and the ellipse re-fitted to the union. `--merge-factor` scales the
overlap test (`>1` also fuses holes that are merely close); `--no-merge-overlap` turns it
off.

The tube ceiling is a thin ribbon with a **ragged rim**, which spawns lots of tiny false
holes along its sides. Use **`--edge-margin D`** to keep only holes whose centre is at
least `D` (coord units) inside the ceiling boundary — it erodes the ribbon footprint
inward by `D` and rejects anything in that border band. Raise it until the rim nicks
disappear, but keep it below the ribbon's half-width or you'll erode real skylights too;
combine with `--min-area` to drop the smallest specks.

### `register` — match constellations (outlier-robust)
```
#   -o / --output      output transform .json (required)
#   -m / --mode        auto (default) | 4dof (up-assisted) | 6dof (Kabsch)
#   -t / --tolerance   max landmark mismatch (coord units) to count as an inlier
#   -n / --min-inliers min mutually-consistent skylights required (default 3)
#        --show        plot the aligned constellations (matches + rejected outliers)
```
This is designed to handle **lots of outlier holes** — small noisy detections or holes
that exist in only one map. It finds correspondences by gating hole pairs/triplets on the
**3-D inter-hole distance** (a rigid invariant, independent of the estimated up-axis), then
keeps the transform that makes the most *other* holes line up (**max consensus**). The key
lever is **`--min-inliers`** (default 3): a stray pairing can't recruit a third hole to
agree, so the outliers are voted out. The ellipse shape (size/axis-ratio/orientation) is
used only to break ties between geometrically equivalent solutions — never to reject a
match.

`auto` tries `4dof` (up-assisted: only 2 inliers needed, handles a straight/collinear
tube) and falls back to `6dof` (full Kabsch, ≥3 non-collinear) if the up-axis proves
unreliable. The command prints the matched count, RMS, the **uniqueness margin** (how many
more holes the best solution explains than the next distinct one — `0` means ambiguous),
the ellipse **shape score**, and warnings. Lower `--min-inliers` to 2 only if you trust a
2-hole match.

**The vertical (Z) caveat.** The skylight openings pin down XY and yaw well, but the
*vertical* offset is only a best fit of the tube's **ceiling** level to the aerial
**surface** level — two different elevations — so it tends to lift the tube toward the
surface. `register` prints the **per-skylight vertical residual** so you can see whether
the openings even agree on depth (all near 0 = consistent; a big spread = a bad up-axis or
genuinely different roof depths). Correct the depth directly with `merge --z-offset` (you
can judge it by eye in the viewer); the constellation XY/yaw stays put.

`--show` plots the constellations **before and after** registration, side by side, looking
down the aerial up-axis: *left* is the raw, unaligned state (tube holes offset/rotated from
the aerial, long green correspondence lines); *right* is after the transform (tube mapped
into the aerial frame, matched ellipses snapping on top of each other). In both, **aerial**
holes are blue, **tube** holes orange, **matched** holes solid and joined by green lines,
rejected outliers faded/dashed — so the leftover singletons are the outliers the consensus
threw out.

### `merge` — apply and combine
```
#   -t / --transform   transform .json from `register` (required)
#        --refine       GICP (small_gicp) on the matched rims before merging
#        --dof          [--refine] 4 (yaw + translation; default) or 6 (full rigid)
#        --rim-inflate  [--refine] grow each skylight's ellipse by this factor to gather its rim
#        --rim-radius   [--refine] fallback rim radius for holes with no ellipse; also the GICP clamp
#        --rim-height   [--refine] vertical band around each opening for GICP
#        --show-rims    plot the rim points GICP operates on (aerial vs tube, before/after)
#   -z / --z-offset    slide the tube vertically (along aerial up) to set its depth
#        --color / --no-color  keep aerial RGB + elevation-colour the tube (default on)
#        --cmap         matplotlib colormap for the tube's elevation colour (default viridis)
#   -s / --source-field  also add a 'source' channel (0=aerial, 1=tube)
#   -c / --chunk-size / -q / --quiet  as elsewhere
```
The tube is transformed into the **aerial** frame and the two clouds are written as a
single `.pcd`. By default (`--color`) the merged cloud has a packed `rgb` field: the
**aerial** points keep their own RGB (grey if the aerial cloud has none), and the **tube**
points are shaded by **elevation** (output-frame Z) with `--cmap` — so the photographic
surface and the depth-coloured tube read distinctly in `pcl_viewer`. `--no-color` writes a
plain `x y z` cloud; `--source-field` adds a `0/1` origin channel either way.

`--refine` runs **GICP** (plane-to-plane, via `small_gicp`) on the matched skylight rims
(the two clouds barely overlap, so a global registration would just flatten the tube onto
the ground). By default it is **4-DOF** (`--dof 4`) — yaw about the aerial up-axis plus
translation, so the tube's *tilt is locked* and it can't be laid flat; pass `--dof 6` for a
full rigid fit. Each rim is gathered **per hole**: points inside that skylight's fitted
ellipse (from `register`) grown by `--rim-inflate`, and within `--rim-height` of the opening
level — so the patch adapts to each hole's size/shape and excludes the deep tube body and
far ground. Holes with no ellipse data fall back to a `--rim-radius` circle (that radius is
also the GICP correspondence clamp). `--rim-height` is the key knob: shrink it until the
refine stops being pulled toward the ground — use `--show-rims` to see exactly which points
each opening's patch is capturing (aerial rim in blue, tube rim before/after refine in
orange/green), so you can tell whether the band is grabbing the deep tube body before you
trust the fit. As a safety net the refinement is **rejected** (the landmark alignment kept,
with a warning) if it would move the matched skylights by more than `--rim-radius` or make
the rim fit worse — so it can never make things dramatically worse. There is also a thin
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
pcl_viewer merged.pcd      # aerial keeps its RGB, the tube is elevation-coloured (rgb field)
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
