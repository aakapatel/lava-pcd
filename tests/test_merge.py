"""Tests for the skylight-based merge pipeline: holes -> register -> merge."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from lava_pcd.geometry import (
    as_matrix,
    rotation_about_axis,
    transform_points,
)
from lava_pcd.holes import Hole, HoleSet, detect_holes, estimate_up
from lava_pcd.io.pcd_writer import BinaryPcdWriter
from lava_pcd.merge import merge_clouds
from lava_pcd.register import match_constellations


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _write_pcd(path: Path, xyz: np.ndarray, origin=(0.0, 0.0, 0.0)) -> None:
    xyz = np.ascontiguousarray(xyz, dtype=np.float32)
    with BinaryPcdWriter(path, max_points=len(xyz), fields=("x", "y", "z"),
                         origin=origin) as w:
        if len(xyz):
            w.write_chunk(xyz)


def _read_pcd_xyz(path: Path) -> tuple[dict, np.ndarray]:
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
        n = int(header["POINTS"])
        ncols = len(header["FIELDS"].split())
        data = np.frombuffer(fh.read(n * ncols * 4), dtype=np.float32).reshape(n, ncols)
    return header, data


def _ground_with_voids(centers, radius=4.0, extent=100.0, spacing=1.0, z=0.0):
    """A dense flat ground (z=const) with circular voids removed at ``centers``."""
    g = np.arange(0.0, extent, spacing)
    gx, gy = np.meshgrid(g, g)
    pts = np.column_stack([gx.ravel(), gy.ravel(), np.full(gx.size, z)])
    keep = np.ones(len(pts), dtype=bool)
    for cx, cy in centers:
        keep &= (pts[:, 0] - cx) ** 2 + (pts[:, 1] - cy) ** 2 > radius ** 2
    return pts[keep]


# skylight centres in the world / aerial frame
_CENTERS = [(25.0, 30.0), (62.0, 35.0), (40.0, 72.0), (78.0, 76.0)]


# --------------------------------------------------------------------------- #
# up-axis estimation
# --------------------------------------------------------------------------- #
def test_estimate_up_recovers_plane_normal() -> None:
    rng = np.random.default_rng(0)
    # A tilted plane with normal n (positive z); estimate_up should recover ~n.
    n = np.array([0.2, -0.3, 0.93])
    n /= np.linalg.norm(n)
    # two in-plane basis vectors
    e1 = np.cross(n, [1.0, 0.0, 0.0]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    uv = rng.uniform(-50, 50, size=(8000, 2))
    pts = uv[:, :1] * e1 + uv[:, 1:] * e2
    pts += rng.normal(0, 0.02, pts.shape)  # thin slab
    up = np.asarray(estimate_up(pts))
    assert abs(float(np.dot(up, n))) > 0.9


# --------------------------------------------------------------------------- #
# constellation matching (unit, no I/O)
# --------------------------------------------------------------------------- #
def _holeset(centroids, up=(0.0, 0.0, 1.0), origin=(0.0, 0.0, 0.0)):
    holes = [Hole(i, tuple(c), 3.0, 28.0, 28) for i, c in enumerate(centroids)]
    return HoleSet(holes, origin, tuple(up), "test", 1.0, "")


def test_match_4dof_recovers_transform_with_outliers() -> None:
    rng = np.random.default_rng(3)
    A = np.column_stack([rng.uniform(0, 100, 6), rng.uniform(0, 100, 6),
                         rng.uniform(-1, 1, 6)])  # roughly planar
    R_g = rotation_about_axis(np.array([0.0, 0.0, 1.0]), 0.9) @ \
        rotation_about_axis(np.array([1.0, 0.0, 0.0]), 0.15)
    t_g = np.array([12.0, -7.0, 3.0])
    # tube = R_g^T (A - t_g); recovering M should give (R_g, t_g)
    B = (A - t_g) @ R_g
    up_tube = R_g.T @ np.array([0.0, 0.0, 1.0])
    B += rng.normal(0, 0.05, B.shape)
    # add two spurious tube holes with no aerial match
    B = np.vstack([B, [[200.0, 5.0, 0.0], [-150.0, 40.0, 2.0]]])

    aerial = _holeset(A, up=(0, 0, 1))
    tube = _holeset(B, up=up_tube)
    tf = match_constellations(aerial, tube, mode="4dof", tolerance=2.0)

    assert len(tf.inliers) == 6
    pred = transform_points(tube.centroids()[[p for p, _ in tf.inliers]], tf.array)
    truth = aerial.centroids()[[a for _, a in tf.inliers]]
    np.testing.assert_allclose(pred, truth, atol=0.5)


def test_match_6dof_recovers_full_rotation() -> None:
    rng = np.random.default_rng(7)
    A = rng.uniform(0, 80, size=(7, 3))  # non-planar
    R_g = rotation_about_axis(np.array([0.4, 0.7, 0.5]), 1.1)
    t_g = np.array([-20.0, 15.0, 8.0])
    B = (A - t_g) @ R_g
    aerial = _holeset(A)
    tube = _holeset(B)
    tf = match_constellations(aerial, tube, mode="6dof", tolerance=1.0)
    assert len(tf.inliers) == 7
    pred = transform_points(tube.centroids()[[p for p, _ in tf.inliers]], tf.array)
    truth = aerial.centroids()[[a for _, a in tf.inliers]]
    np.testing.assert_allclose(pred, truth, atol=0.5)
    # recovered transform maps tube -> aerial: A = R_g @ B + t_g
    np.testing.assert_allclose(tf.array, as_matrix(R_g, t_g), atol=1e-2)


def test_match_too_few_holes_raises() -> None:
    aerial = _holeset([(0, 0, 0)])
    tube = _holeset([(1, 1, 1)])
    with pytest.raises(ValueError):
        match_constellations(aerial, tube)


# --------------------------------------------------------------------------- #
# detection
# --------------------------------------------------------------------------- #
def test_detect_aerial_voids(tmp_path: Path) -> None:
    pcd = tmp_path / "aerial.pcd"
    _write_pcd(pcd, _ground_with_voids(_CENTERS))
    hs = detect_holes(pcd, mode="aerial", resolution=1.0, min_area=10.0)
    assert len(hs.holes) == len(_CENTERS)
    found = sorted((h.centroid[0], h.centroid[1]) for h in hs.holes)
    for (cx, cy), (fx, fy) in zip(sorted(_CENTERS), found):
        assert abs(fx - cx) < 1.5 and abs(fy - cy) < 1.5


# --------------------------------------------------------------------------- #
# end-to-end: holes -> register -> merge
# --------------------------------------------------------------------------- #
def test_pipeline_holes_register_merge(tmp_path: Path) -> None:
    # ground-truth transform mapping tube-local -> aerial(world)
    R_g = rotation_about_axis(np.array([0.0, 0.0, 1.0]), 1.3) @ \
        rotation_about_axis(np.array([1.0, 0.0, 0.0]), 0.2)
    t_g = np.array([15.0, -25.0, 5.0])
    up_tube = tuple(R_g.T @ np.array([0.0, 0.0, 1.0]))

    world = _ground_with_voids(_CENTERS)
    tube_local = (world - t_g) @ R_g

    aerial_pcd = tmp_path / "aerial.pcd"
    tube_pcd = tmp_path / "tube.pcd"
    _write_pcd(aerial_pcd, world)
    _write_pcd(tube_pcd, tube_local)

    aerial_holes = detect_holes(aerial_pcd, mode="aerial", resolution=1.0, min_area=10.0)
    tube_holes = detect_holes(tube_pcd, mode="ceiling", up=up_tube, resolution=1.0,
                              min_area=10.0)
    assert len(aerial_holes.holes) == len(_CENTERS)
    assert len(tube_holes.holes) == len(_CENTERS)

    tf = match_constellations(aerial_holes, tube_holes, mode="4dof", tolerance=2.0)
    assert len(tf.inliers) >= 3

    # matched tube skylights land on the aerial ones after the transform
    pred = transform_points(
        tube_holes.centroids()[[p for p, _ in tf.inliers]], tf.array
    )
    truth = aerial_holes.centroids()[[a for _, a in tf.inliers]]
    np.testing.assert_allclose(pred, truth, atol=1.5)

    out = tmp_path / "merged.pcd"
    res = merge_clouds(aerial_pcd, tube_pcd, out, tf, source_field=True,
                       show_progress=False)
    assert res.point_count == res.aerial_count + res.tube_count
    header, data = _read_pcd_xyz(out)
    assert header["FIELDS"] == "x y z source"
    assert data.shape[0] == res.point_count
    # the source channel splits the two clouds
    assert set(np.unique(data[:, 3])) == {0.0, 1.0}


def test_merge_rejects_output_equal_input(tmp_path: Path) -> None:
    pcd = tmp_path / "a.pcd"
    _write_pcd(pcd, _ground_with_voids(_CENTERS))
    tf = match_constellations(
        _holeset([(0, 0, 0), (10, 0, 0), (0, 10, 0)]),
        _holeset([(0, 0, 0), (10, 0, 0), (0, 10, 0)]),
    )
    with pytest.raises(ValueError):
        merge_clouds(pcd, tmp_path / "b.pcd", pcd, tf, show_progress=False)
