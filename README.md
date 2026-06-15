# lava-pcd

Tools for processing large point cloud files — convert, downsample, voxelize, crop, merge.

The first capability is converting `.laz`/`.las` LiDAR clouds into **binary `.pcd`** files
(fields `x y z intensity`) that work both in Python and with `pcl_viewer`.

## Install

```bash
pip install -e .
```

This pulls `laspy[lazrs]` (a pure-Rust LAZ backend, so no external `laszip` is needed).

## Convert a file

```bash
lava-pcd convert input.laz output.pcd
# options:
#   -f / --fields      attributes to export: auto, xyz, intensity, rgb, all (default auto)
#   -v / --voxel       voxel-downsample resolution in coord units (0 = off).
#                      Prompted interactively if not given.
#   -c / --chunk-size  points read per chunk (lower = less memory; default 5,000,000)
#   -o / --origin      local-origin shift: 'header' (default), 'none', or 'x,y,z'
#       --parallel     use the multi-threaded LAZ backend (faster; panics on some files)
#   -q / --quiet       suppress the progress bar
```

`--fields auto` exports **RGB** when the file is colourised (e.g. photogrammetry /
aerial products, where LiDAR intensity is often empty), otherwise **intensity**. RGB is
written as PCL's packed `rgb` float field, which `pcl_viewer` and Open3D colour by
automatically. Use `intensity`, `rgb`, `xyz`, or `all` to force a choice.

`--voxel`/`-v` downsamples to one centroid per cubic voxel of the given edge length
(in coordinate units, e.g. metres). `0` disables it. Downsampling is streamed, so it
stays memory-bounded even for hundred-million-point clouds.

From Python:

```python
from lava_pcd import laz_to_pcd

res = laz_to_pcd("input.laz", "output.pcd", voxel_size=0.1, chunk_size=2_000_000)
print(f"{res.source_count} -> {res.point_count} points, origin {res.origin}")
```

## Viewing

```bash
pcl_viewer output.pcd      # press 2 to colour by the intensity field
```

## Notes

- **Reading is chunked**, so memory stays bounded for very large clouds.
- Coordinates are stored as **float32**. This is fine for typical LiDAR ranges but loses
  precision for large survey-grade (e.g. raw UTM) coordinates. If you need full precision,
  a future option can subtract a local origin before writing.
- Only `intensity` is preserved alongside XYZ for now; RGB and other LAS dimensions are
  ignored.

## Development

```bash
pip install -e ".[dev]"
pytest -q
```
