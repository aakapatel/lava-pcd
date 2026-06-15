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
#   -c / --chunk-size  points read per chunk (lower = less memory; default 5,000,000)
#   -q / --quiet       suppress the progress bar
```

From Python:

```python
from lava_pcd import laz_to_pcd

n = laz_to_pcd("input.laz", "output.pcd", chunk_size=2_000_000)
print(f"{n} points written")
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
