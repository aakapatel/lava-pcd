"""Round-trip test: synthetic .las -> .pcd -> parsed back."""

from __future__ import annotations

from pathlib import Path

import laspy
import numpy as np
import pytest

from lava_pcd.convert import laz_to_pcd
from lava_pcd.io.laz_reader import LazChunkReader, in_bounds_mask


def _make_las(
    path: Path,
    xyz: np.ndarray,
    intensity: np.ndarray,
    rgb: np.ndarray | None = None,
    offsets=(0.0, 0.0, 0.0),
) -> None:
    # Format 0 = xyz + intensity (no RGB); format 2 adds RGB.
    point_format = 2 if rgb is not None else 0
    header = laspy.LasHeader(point_format=point_format, version="1.2")
    # Fine scales so float32 round-trip stays well within tolerance.
    header.scales = [0.001, 0.001, 0.001]
    header.offsets = list(offsets)
    las = laspy.LasData(header)
    las.x, las.y, las.z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    las.intensity = intensity
    if rgb is not None:
        las.red, las.green, las.blue = rgb[:, 0], rgb[:, 1], rgb[:, 2]
    las.write(path)


def _read_pcd(path: Path) -> tuple[dict, np.ndarray]:
    """Parse a binary PCD into (header dict, (N, ncols) float32 array)."""
    with open(path, "rb") as fh:
        header = {}
        while True:
            line = fh.readline().decode("ascii").strip()
            if not line or line.startswith("#"):
                continue
            key, _, value = line.partition(" ")
            header[key] = value
            if key == "DATA":
                break
        assert header["DATA"] == "binary"
        n = int(header["POINTS"])
        ncols = len(header["FIELDS"].split())
        data = np.frombuffer(fh.read(n * ncols * 4), dtype=np.float32).reshape(n, ncols)
    return header, data


def test_laz_to_pcd_roundtrip(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    n = 1234
    xyz = rng.uniform(-50.0, 50.0, size=(n, 3))
    intensity = rng.integers(0, 65535, size=n, dtype=np.uint16)

    las_path = tmp_path / "in.las"
    pcd_path = tmp_path / "out.pcd"
    _make_las(las_path, xyz, intensity)

    # Use a small chunk size to exercise the chunked path. Header offset is 0,
    # so origin="header" is a no-op here.
    result = laz_to_pcd(las_path, pcd_path, chunk_size=500, show_progress=False)
    assert result.point_count == n
    assert result.origin == (0.0, 0.0, 0.0)
    assert result.sidecar_path is not None and result.sidecar_path.exists()

    header, data = _read_pcd(pcd_path)
    assert header["FIELDS"] == "x y z intensity"
    assert int(header["POINTS"]) == n
    assert data.shape == (n, 4)

    np.testing.assert_allclose(data[:, :3], xyz, atol=1e-2)
    np.testing.assert_array_equal(data[:, 3].astype(np.uint16), intensity)


def test_origin_shift_preserves_precision(tmp_path: Path) -> None:
    """Large UTM-like coords must survive via the local-origin shift."""
    rng = np.random.default_rng(2)
    n = 1000
    base = np.array([479157.0, 7089825.0, 100.0])
    local = rng.uniform(0.0, 500.0, size=(n, 3))  # within float32 mm range
    xyz = base + local
    intensity = rng.integers(0, 65535, size=n, dtype=np.uint16)

    las_path = tmp_path / "utm.las"
    pcd_path = tmp_path / "utm.pcd"
    # LAS header offset near the data -> used as origin by default.
    _make_las(las_path, xyz, intensity, offsets=tuple(base))

    result = laz_to_pcd(las_path, pcd_path, show_progress=False)
    np.testing.assert_allclose(result.origin, base, atol=1e-6)

    _, data = _read_pcd(pcd_path)
    # Reconstruct global coords; mm-level accuracy must hold.
    recovered = data[:, :3].astype(np.float64) + np.array(result.origin)
    np.testing.assert_allclose(recovered, xyz, atol=1e-2)


def test_voxel_downsample(tmp_path: Path) -> None:
    # Two tight clusters > 1 voxel apart -> downsamples to ~2 points.
    rng = np.random.default_rng(3)
    a = np.array([1.0, 1.0, 1.0]) + rng.uniform(0, 0.05, size=(500, 3))
    b = np.array([20.0, 20.0, 20.0]) + rng.uniform(0, 0.05, size=(500, 3))
    xyz = np.vstack([a, b])
    intensity = rng.integers(0, 65535, size=len(xyz), dtype=np.uint16)

    las_path = tmp_path / "clusters.las"
    pcd_path = tmp_path / "clusters.pcd"
    _make_las(las_path, xyz, intensity)

    result = laz_to_pcd(las_path, pcd_path, voxel_size=1.0, show_progress=False)
    assert result.source_count == 1000
    assert result.point_count == 2  # one centroid per cluster
    assert result.voxel_size == 1.0

    header, data = _read_pcd(pcd_path)
    assert int(header["POINTS"]) == 2
    centroids = np.sort(data[:, :3], axis=0)
    expected = np.sort(np.vstack([a.mean(0), b.mean(0)]), axis=0)
    np.testing.assert_allclose(centroids, expected, atol=0.05)


def test_no_downsample_when_zero(tmp_path: Path) -> None:
    rng = np.random.default_rng(4)
    xyz = rng.uniform(-10, 10, size=(300, 3))
    intensity = rng.integers(0, 65535, size=300, dtype=np.uint16)
    las_path = tmp_path / "in.las"
    pcd_path = tmp_path / "out.pcd"
    _make_las(las_path, xyz, intensity)

    result = laz_to_pcd(las_path, pcd_path, voxel_size=0.0, show_progress=False)
    assert result.point_count == 300


def _unpack_rgb(rgb_float: np.ndarray) -> np.ndarray:
    packed = rgb_float.view(np.uint32)
    r = (packed >> 16) & 0xFF
    g = (packed >> 8) & 0xFF
    b = packed & 0xFF
    return np.column_stack([r, g, b]).astype(np.uint8)


def test_rgb_auto_export(tmp_path: Path) -> None:
    """A colourised cloud should export packed RGB by default ('auto')."""
    rng = np.random.default_rng(5)
    n = 800
    xyz = rng.uniform(-10, 10, size=(n, 3))
    intensity = np.zeros(n, dtype=np.uint16)  # empty intensity, like the real file
    rgb8 = rng.integers(0, 256, size=(n, 3), dtype=np.uint16)
    rgb16 = (rgb8 * 257).astype(np.uint16)  # 8-bit -> 16-bit, exact round-trip

    las_path = tmp_path / "color.las"
    pcd_path = tmp_path / "color.pcd"
    _make_las(las_path, xyz, intensity, rgb=rgb16)

    result = laz_to_pcd(las_path, pcd_path, show_progress=False)
    assert result.fields == ("x", "y", "z", "rgb")

    header, data = _read_pcd(pcd_path)
    assert header["FIELDS"] == "x y z rgb"
    np.testing.assert_allclose(data[:, :3], xyz, atol=1e-2)
    np.testing.assert_array_equal(_unpack_rgb(data[:, 3]), rgb8.astype(np.uint8))


def test_rgb_with_voxel_downsample(tmp_path: Path) -> None:
    """RGB must survive downsampling (averaged per voxel, then re-packed)."""
    rng = np.random.default_rng(6)
    # One tight cluster, all the same colour -> single point keeps that colour.
    xyz = np.array([5.0, 5.0, 5.0]) + rng.uniform(0, 0.05, size=(400, 3))
    intensity = np.zeros(len(xyz), dtype=np.uint16)
    rgb16 = np.full((len(xyz), 3), [100 * 257, 150 * 257, 200 * 257], dtype=np.uint16)

    las_path = tmp_path / "c.las"
    pcd_path = tmp_path / "c.pcd"
    _make_las(las_path, xyz, intensity, rgb=rgb16)

    result = laz_to_pcd(las_path, pcd_path, voxel_size=1.0, show_progress=False)
    assert result.point_count == 1
    _, data = _read_pcd(pcd_path)
    np.testing.assert_array_equal(_unpack_rgb(data[:, 3])[0], [100, 150, 200])


def test_missing_input(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        laz_to_pcd(tmp_path / "nope.laz", tmp_path / "out.pcd")


def test_bad_suffix(tmp_path: Path) -> None:
    bogus = tmp_path / "in.txt"
    bogus.write_text("not a point cloud")
    with pytest.raises(ValueError):
        laz_to_pcd(bogus, tmp_path / "out.pcd")


def test_in_bounds_mask() -> None:
    gx = np.array([0.0, 5.0, 1e6, -1e6])
    gy = np.array([0.0, 5.0, 0.0, 0.0])
    gz = np.array([0.0, 5.0, 0.0, 0.0])
    mask = in_bounds_mask(gx, gy, gz, (0, 0, 0), (10, 10, 10))
    np.testing.assert_array_equal(mask, [True, True, False, False])


def test_reader_drops_out_of_bounds(tmp_path: Path) -> None:
    # 100 good points in [0,10] plus 5 sentinel points far outside.
    good = np.random.default_rng(8).uniform(0, 10, size=(100, 3))
    bad = np.full((5, 3), 2_147_483.0)  # INT32-limit-style sentinels
    xyz = np.vstack([good, bad])
    intensity = np.zeros(len(xyz), dtype=np.uint16)
    las_path = tmp_path / "withbad.las"
    _make_las(las_path, xyz, intensity)

    with LazChunkReader(las_path, filter_bounds=True) as r:
        # laspy recomputes header bounds to include the sentinels on write, so
        # tighten them to the real data range for the test.
        r._reader.header.maxs = [10.0, 10.0, 10.0]
        r._reader.header.mins = [0.0, 0.0, 0.0]
        kept = sum(len(c) for c in r.chunks(("x", "y", "z")))
    assert kept == 100
    assert r.dropped == 5


def test_writer_count_rewrite_when_fewer(tmp_path: Path) -> None:
    """Header POINTS must reflect actual count even if < reserved max_points."""
    from lava_pcd.io.pcd_writer import BinaryPcdWriter

    pcd = tmp_path / "few.pcd"
    pts = np.arange(12, dtype=np.float32).reshape(3, 4)  # 3 points, x y z intensity
    with BinaryPcdWriter(pcd, max_points=1_000_000, fields=("x", "y", "z", "intensity")) as w:
        w.write_chunk(pts)

    header, data = _read_pcd(pcd)
    assert int(header["POINTS"]) == 3  # may be space-padded to reserved width
    assert int(header["WIDTH"]) == 3
    assert data.shape == (3, 4)
    np.testing.assert_array_equal(data, pts)


def test_reproject_to_utm(tmp_path: Path) -> None:
    """Lon/lat (EPSG:4326) input reprojects to metric UTM with recoverable origin."""
    import pyproj

    # A few points near Iceland (lon ~ -21.4, lat ~ 63.9), Z in metres.
    lon = np.array([-21.40, -21.39, -21.38])
    lat = np.array([63.930, 63.931, 63.932])
    z = np.array([214.0, 250.0, 300.0])
    xyz = np.column_stack([lon, lat, z])
    intensity = np.zeros(len(xyz), dtype=np.uint16)

    header = laspy.LasHeader(point_format=0, version="1.4")
    header.scales = [1e-7, 1e-7, 1e-3]  # fine enough for degrees
    header.offsets = [-21.0, 63.0, 0.0]
    header.add_crs(pyproj.CRS.from_epsg(4326))
    las = laspy.LasData(header)
    las.x, las.y, las.z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    las.intensity = intensity
    las_path = tmp_path / "geo.las"
    las.write(las_path)

    pcd_path = tmp_path / "utm.pcd"
    result = laz_to_pcd(
        las_path, pcd_path, fields="xyz", reproject="EPSG:32627", show_progress=False
    )
    assert result.reproject == "EPSG:32627"

    _, data = _read_pcd(pcd_path)
    # Reconstruct global UTM coords and compare to a direct pyproj transform.
    t = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32627", always_xy=True)
    ex, ny = t.transform(lon, lat)
    recovered = data[:, :3].astype(np.float64) + np.array(result.origin)
    np.testing.assert_allclose(recovered[:, 0], ex, atol=0.05)
    np.testing.assert_allclose(recovered[:, 1], ny, atol=0.05)
    np.testing.assert_allclose(recovered[:, 2], z, atol=0.05)
    # Origin should be near the data (UTM 27N easting ~480k, northing ~7.09M).
    assert 470_000 < result.origin[0] < 490_000
    assert 7_080_000 < result.origin[1] < 7_100_000
