"""Round-trip test: synthetic .las -> .pcd -> parsed back."""

from __future__ import annotations

import struct
from pathlib import Path

import laspy
import numpy as np
import pytest

from lava_pcd.convert import laz_to_pcd


def _make_las(path: Path, xyz: np.ndarray, intensity: np.ndarray) -> None:
    header = laspy.LasHeader(point_format=3, version="1.2")
    # Fine scales so float32 round-trip stays well within tolerance.
    header.scales = [0.001, 0.001, 0.001]
    header.offsets = [0.0, 0.0, 0.0]
    las = laspy.LasData(header)
    las.x = xyz[:, 0]
    las.y = xyz[:, 1]
    las.z = xyz[:, 2]
    las.intensity = intensity
    las.write(path)


def _read_pcd(path: Path) -> tuple[dict, np.ndarray]:
    """Parse a binary PCD into (header dict, (N, 4) float32 array)."""
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
        data = np.frombuffer(fh.read(n * 4 * 4), dtype=np.float32).reshape(n, 4)
    return header, data


def test_laz_to_pcd_roundtrip(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    n = 1234
    xyz = rng.uniform(-50.0, 50.0, size=(n, 3))
    intensity = rng.integers(0, 65535, size=n, dtype=np.uint16)

    las_path = tmp_path / "in.las"
    pcd_path = tmp_path / "out.pcd"
    _make_las(las_path, xyz, intensity)

    # Use a small chunk size to exercise the chunked path.
    count = laz_to_pcd(las_path, pcd_path, chunk_size=500, show_progress=False)
    assert count == n

    header, data = _read_pcd(pcd_path)
    assert header["FIELDS"] == "x y z intensity"
    assert int(header["POINTS"]) == n
    assert data.shape == (n, 4)

    np.testing.assert_allclose(data[:, :3], xyz, atol=1e-2)
    np.testing.assert_array_equal(data[:, 3].astype(np.uint16), intensity)


def test_missing_input(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        laz_to_pcd(tmp_path / "nope.laz", tmp_path / "out.pcd")


def test_bad_suffix(tmp_path: Path) -> None:
    bogus = tmp_path / "in.txt"
    bogus.write_text("not a point cloud")
    with pytest.raises(ValueError):
        laz_to_pcd(bogus, tmp_path / "out.pcd")
