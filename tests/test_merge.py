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
import matplotlib
matplotlib.use("Agg")  # headless: plt.show() is a no-op, so the viewer is testable

from lava_pcd.holes import (
    Hole,
    HoleSet,
    build_occupancy,
    detect_holes,
    detect_skylights,
    estimate_up,
    show_occupancy,
)
from lava_pcd.io.pcd_writer import BinaryPcdWriter
from lava_pcd.merge import icp_refine, merge_clouds, visualize_rims
from lava_pcd.register import (
    Transform,
    match_constellations,
    vertical_residuals,
    visualize_match,
)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _write_pcd(path: Path, xyz: np.ndarray, origin=(0.0, 0.0, 0.0)) -> None:
    xyz = np.ascontiguousarray(xyz, dtype=np.float32)
    with BinaryPcdWriter(path, max_points=len(xyz), fields=("x", "y", "z"),
                         origin=origin) as w:
        if len(xyz):
            w.write_chunk(xyz)


def _write_pcd_rgb(path: Path, xyz: np.ndarray, rgb_packed: np.ndarray) -> None:
    arr = np.column_stack([xyz.astype(np.float32), rgb_packed.astype(np.float32)])
    with BinaryPcdWriter(path, max_points=len(arr), fields=("x", "y", "z", "rgb")) as w:
        if len(arr):
            w.write_chunk(np.ascontiguousarray(arr, dtype=np.float32))


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


def _holeset_shaped(centroids, areas, up=(0.0, 0.0, 1.0)):
    holes = []
    for i, (c, ar) in enumerate(zip(centroids, areas)):
        r = float(np.sqrt(ar / np.pi))
        holes.append(Hole(i, tuple(c), r, float(ar), int(ar),
                          semi_major=r, semi_minor=r, orientation=0.0))
    return HoleSet(holes, (0.0, 0.0, 0.0), tuple(up), "test", 1.0, "")


def test_match_robust_to_many_outliers() -> None:
    rng = np.random.default_rng(11)
    # 3 true holes (non-collinear); the rest are outliers in each set
    A_true = np.array([[10.0, 10.0, 0.0], [40.0, 15.0, 1.0], [25.0, 45.0, -1.0]])
    R_g = rotation_about_axis(np.array([0.0, 0.0, 1.0]), 0.8) @ \
        rotation_about_axis(np.array([1.0, 0.0, 0.0]), 0.1)
    t_g = np.array([12.0, -8.0, 4.0])
    up_tube = tuple(R_g.T @ np.array([0.0, 0.0, 1.0]))
    B_true = (A_true - t_g) @ R_g

    A_out = rng.uniform(-60.0, 120.0, size=(8, 3))
    B_out = rng.uniform(-120.0, 120.0, size=(9, 3))
    A = np.vstack([A_true, A_out])            # true holes at indices 0,1,2
    B = np.vstack([B_true, B_out])
    aerial = _holeset(A, up=(0, 0, 1))
    tube = _holeset(B, up=up_tube)

    tf = match_constellations(aerial, tube, mode="4dof", tolerance=1.0, min_inliers=3)
    assert {(0, 0), (1, 1), (2, 2)} <= set(tf.inliers)   # true correspondences found
    assert tf.n_inliers == 3 and tf.margin >= 1          # outliers excluded, unique
    pred = transform_points(tube.centroids()[[0, 1, 2]], tf.array)
    np.testing.assert_allclose(pred, A_true, atol=0.5)


def test_match_6dof_robust_to_outliers() -> None:
    rng = np.random.default_rng(5)
    A_true = np.array([[0.0, 0.0, 0.0], [30.0, 5.0, 10.0], [12.0, 28.0, -6.0]])
    R_g = rotation_about_axis(np.array([0.3, 0.6, 0.5]), 1.0)
    t_g = np.array([-15.0, 20.0, 7.0])
    B_true = (A_true - t_g) @ R_g
    A = np.vstack([A_true, rng.uniform(-80, 120, size=(7, 3))])
    B = np.vstack([B_true, rng.uniform(-120, 120, size=(8, 3))])
    tf = match_constellations(_holeset(A), _holeset(B), mode="6dof",
                              tolerance=1.0, min_inliers=3)
    assert {(0, 0), (1, 1), (2, 2)} <= set(tf.inliers)
    np.testing.assert_allclose(
        transform_points(_holeset(B).centroids()[[0, 1, 2]], tf.array), A_true, atol=0.5)


def test_match_too_few_inliers_raises() -> None:
    # only 2 holes truly correspond; with min_inliers=3 this must fail
    A_true = np.array([[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]])
    R_g = rotation_about_axis(np.array([0.0, 0.0, 1.0]), 0.5)
    B_true = A_true @ R_g
    A = np.vstack([A_true, [[100.0, 100.0, 0.0]]])
    B = np.vstack([B_true, [[-90.0, 60.0, 0.0]]])
    with pytest.raises(ValueError):
        match_constellations(_holeset(A), _holeset(B), mode="6dof", min_inliers=3)


def test_match_ambiguous_constellation_warns() -> None:
    # equilateral triangle in both -> several correspondences explain all 3 equally
    tri = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [5.0, 8.66, 0.0]])
    tf = match_constellations(_holeset(tri), _holeset(tri), mode="6dof",
                              tolerance=0.5, min_inliers=3)
    assert tf.n_inliers == 3
    assert tf.margin <= 0
    assert any("ambiguous" in w for w in tf.warnings)


def test_register_visualize_headless() -> None:
    A_true = np.array([[10.0, 10.0, 0.0], [40.0, 15.0, 1.0], [25.0, 45.0, -1.0]])
    R_g = rotation_about_axis(np.array([0.0, 0.0, 1.0]), 0.6) @ \
        rotation_about_axis(np.array([1.0, 0.0, 0.0]), 0.1)
    t_g = np.array([5.0, -5.0, 2.0])
    up_tube = tuple(R_g.T @ np.array([0.0, 0.0, 1.0]))
    B_true = (A_true - t_g) @ R_g
    aerial = _holeset(np.vstack([A_true, [[100.0, 100.0, 0.0]]]))   # + 1 outlier
    tube = _holeset(np.vstack([B_true, [[-80.0, 50.0, 0.0]]]), up=up_tube)
    tf = match_constellations(aerial, tube, tolerance=1.0)
    visualize_match(aerial, tube, tf)  # must not raise under the Agg backend


def test_match_shape_breaks_geometric_tie() -> None:
    # same equilateral geometry (ambiguous), but distinct hole sizes pick identity
    tri = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [5.0, 8.66, 0.0]])
    areas = [10.0, 50.0, 200.0]
    tf = match_constellations(_holeset_shaped(tri, areas), _holeset_shaped(tri, areas),
                              mode="6dof", tolerance=0.5, min_inliers=3)
    # size agreement makes the identity correspondence the best of the ties
    assert set(tf.inliers) == {(0, 0), (1, 1), (2, 2)}


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


def _ground_with_sparse_hole(center, radius=6.0, extent=60.0, spacing=0.5):
    """Dense ground (4 pts/cell) except a disc thinned to exactly 1 pt/cell.

    Inside the disc only the integer-coordinate points survive, so every 1.0-unit
    cell there keeps exactly one point: occupied at min_density=1 (hole hidden),
    empty at min_density=2 (hole revealed) -- with no fully empty cells to give it
    away.
    """
    g = np.arange(0.0, extent, spacing)
    gx, gy = np.meshgrid(g, g)
    pts = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
    inside = (pts[:, 0] - center[0]) ** 2 + (pts[:, 1] - center[1]) ** 2 < radius ** 2
    on_coarse = np.isclose(pts[:, 0] % 1.0, 0.0) & np.isclose(pts[:, 1] % 1.0, 0.0)
    keep = ~inside | on_coarse
    return pts[keep]


def test_density_threshold_catches_less_occupied_hole(tmp_path: Path) -> None:
    # 4 points/cell on the ground; exactly 1 point/cell inside the thinned disc.
    pcd = tmp_path / "sparse.pcd"
    _write_pcd(pcd, _ground_with_sparse_hole((30.0, 30.0)))

    # min_density=1: the disc still has returns, so it stays "occupied" -> missed.
    miss = detect_holes(pcd, mode="aerial", resolution=1.0, min_area=15.0, min_density=1)
    assert len(miss.holes) == 0
    # min_density=2: the thinned disc drops below threshold -> detected as a void.
    hit = detect_holes(pcd, mode="aerial", resolution=1.0, min_area=15.0, min_density=2)
    assert len(hit.holes) == 1
    cx, cy, _ = hit.holes[0].centroid
    assert abs(cx - 30.0) < 2.0 and abs(cy - 30.0) < 2.0


def _ground_density_gradient(hole_center=(60.0, 40.0), hole_r=5.0,
                             extent=80.0, spacing=0.5):
    """Dense ground for x<40 (4 pts/cell), sparse for x>=40 (1 pt/cell), with a
    hole in the *sparse* half. No single absolute threshold separates the hole
    from the sparse ground; a relative one does."""
    g = np.arange(0.0, extent, spacing)
    gx, gy = np.meshgrid(g, g)
    pts = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
    dense = pts[:, 0] < 40.0
    coarse = np.isclose(pts[:, 0] % 1.0, 0.0) & np.isclose(pts[:, 1] % 1.0, 0.0)
    keep = dense | coarse
    inside = ((pts[:, 0] - hole_center[0]) ** 2
              + (pts[:, 1] - hole_center[1]) ** 2 < hole_r ** 2)
    keep &= ~inside
    return pts[keep]


def test_relative_threshold_handles_varying_density(tmp_path: Path) -> None:
    pcd = tmp_path / "gradient.pcd"
    _write_pcd(pcd, _ground_density_gradient())

    # Absolute threshold tuned for the dense half marks the whole sparse half
    # "empty" (border-connected), so the real hole isn't isolated -> missed.
    absolute = detect_holes(pcd, mode="aerial", resolution=1.0, min_area=10.0,
                            min_density=2)
    assert all(abs(h.centroid[0] - 60.0) > 5 or abs(h.centroid[1] - 40.0) > 5
               for h in absolute.holes)

    # Relative threshold adapts to the local density and finds the hole.
    rel = detect_holes(pcd, mode="aerial", resolution=1.0, min_area=10.0,
                       relative=0.4, relative_window=21)
    near = [h for h in rel.holes
            if abs(h.centroid[0] - 60.0) < 2 and abs(h.centroid[1] - 40.0) < 2]
    assert len(near) == 1


def _slab_with_voids(centers_r, extent=40.0, spacing=0.5):
    """A solid square slab with circular voids removed at ``centers_r``."""
    g = np.arange(0.0, extent, spacing)
    gx, gy = np.meshgrid(g, g)
    pts = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
    keep = np.ones(len(pts), dtype=bool)
    for (cx, cy), r in centers_r:
        keep &= (pts[:, 0] - cx) ** 2 + (pts[:, 1] - cy) ** 2 > r ** 2
    return pts[keep]


def test_edge_margin_drops_boundary_holes(tmp_path: Path) -> None:
    # one genuine interior hole + one nick near the bottom edge (both enclosed)
    pcd = tmp_path / "slab.pcd"
    _write_pcd(pcd, _slab_with_voids([((20.0, 20.0), 3.0), ((20.0, 2.5), 1.6)]))

    both = detect_holes(pcd, mode="aerial", resolution=1.0, min_area=5.0, edge_margin=0)
    assert len(both.holes) == 2

    inner = detect_holes(pcd, mode="aerial", resolution=1.0, min_area=5.0, edge_margin=6)
    assert len(inner.holes) == 1
    cx, cy, _ = inner.holes[0].centroid
    assert abs(cx - 20.0) < 2 and abs(cy - 20.0) < 2


def _slab_with_ellipse_void(center, a, b, angle_deg, extent=60.0, spacing=0.5):
    """Solid slab with one rotated elliptical void (semi-axes a>b at angle_deg)."""
    g = np.arange(0.0, extent, spacing)
    gx, gy = np.meshgrid(g, g)
    pts = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
    th = np.radians(angle_deg)
    dx, dy = pts[:, 0] - center[0], pts[:, 1] - center[1]
    xr = dx * np.cos(th) + dy * np.sin(th)     # major axis
    yr = -dx * np.sin(th) + dy * np.cos(th)    # minor axis
    inside = (xr / a) ** 2 + (yr / b) ** 2 < 1.0
    return pts[~inside]


def test_ellipse_fit_recovers_shape_and_orientation(tmp_path: Path) -> None:
    pcd = tmp_path / "ellipse.pcd"
    _write_pcd(pcd, _slab_with_ellipse_void((30.0, 30.0), a=10.0, b=4.0, angle_deg=30.0))
    hs = detect_holes(pcd, mode="aerial", resolution=0.5, min_area=20.0)
    assert len(hs.holes) == 1
    h = hs.holes[0]
    # equivalent ellipse of a uniform fill: semi-axis ~ the true semi-axis
    assert 8.0 < h.semi_major < 12.0
    assert 2.5 < h.semi_minor < 5.5
    assert h.semi_major > h.semi_minor
    # orientation ~30 deg, modulo 180
    deg = np.degrees(h.orientation) % 180.0
    assert min(abs(deg - 30.0), abs(deg - 210.0), abs(deg + 150.0)) < 12.0


def test_merge_overlapping_ellipses(tmp_path: Path) -> None:
    # two separate voids (a band of occupied ground between them), centres 8 apart
    pcd = tmp_path / "pair.pcd"
    _write_pcd(pcd, _slab_with_voids([((18.0, 20.0), 3.0), ((26.0, 20.0), 3.0)]))

    # distinct connected components, so without merging there are two holes
    split = detect_holes(pcd, mode="aerial", resolution=1.0, min_area=10.0,
                        merge_overlap=False)
    assert len(split.holes) == 2

    # ellipse radii ~3 each, 8 apart: factor 1 doesn't reach, a larger factor fuses
    near = detect_holes(pcd, mode="aerial", resolution=1.0, min_area=10.0,
                       merge_overlap=True, merge_factor=1.0)
    assert len(near.holes) == 2
    fused = detect_holes(pcd, mode="aerial", resolution=1.0, min_area=10.0,
                        merge_overlap=True, merge_factor=1.6)
    assert len(fused.holes) == 1
    cx, cy, _ = fused.holes[0].centroid
    assert abs(cx - 22.0) < 2.0 and abs(cy - 20.0) < 2.0       # midway between the two
    assert fused.holes[0].semi_major > split.holes[0].semi_major  # elongated union


def test_occupancy_viewer_runs_headless(tmp_path: Path) -> None:
    pcd = tmp_path / "aerial.pcd"
    _write_pcd(pcd, _ground_with_voids(_CENTERS))
    grid = build_occupancy(pcd, mode="aerial", resolution=1.0)
    assert grid.counts.shape[0] > 0 and grid.extent[1] > grid.extent[0]
    holes, mask = detect_skylights(grid, min_area=10.0)
    assert len(holes) == len(_CENTERS)
    hs = HoleSet(holes, grid.origin, grid.up, "aerial", 1.0, str(pcd))
    # non-interactive viewer returns the holeset unchanged and doesn't raise
    out = show_occupancy(grid, hs, mask, interactive=False)
    assert out is hs


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
    # colour by default: aerial keeps rgb (grey here, no source rgb), tube by elevation
    assert header["FIELDS"] == "x y z rgb source"
    assert data.shape[0] == res.point_count
    # the source channel (now col 4) splits the two clouds
    assert set(np.unique(data[:, 4])) == {0.0, 1.0}
    # tube points (source 1) carry a spread of elevation colours, not one constant
    tube_rgb = data[data[:, 4] == 1.0][:, 3]
    assert len(np.unique(tube_rgb)) > 1


def test_merge_preserves_aerial_rgb_and_colors_tube(tmp_path: Path) -> None:
    from lava_pcd.register import Transform

    rng = np.random.default_rng(0)
    A = rng.uniform(0.0, 20.0, size=(200, 3)); A[:, 2] = 0.0
    red = np.array([(255 << 16)], dtype=np.uint32).view(np.float32)[0]   # packed (255,0,0)
    aerial_pcd = tmp_path / "aerial.pcd"
    _write_pcd_rgb(aerial_pcd, A, np.full(len(A), red, dtype=np.float32))

    B = rng.uniform(0.0, 20.0, size=(200, 3)); B[:, 2] = rng.uniform(-10.0, -2.0, 200)
    tube_pcd = tmp_path / "tube.pcd"
    _write_pcd(tube_pcd, B)

    tf = Transform(np.eye(4).tolist(), [], 0.0, "4dof", (0, 0, 0), (0, 0, 0))
    out = tmp_path / "merged.pcd"
    res = merge_clouds(aerial_pcd, tube_pcd, out, tf, color=True, show_progress=False)

    header, data = _read_pcd_xyz(out)
    assert header["FIELDS"] == "x y z rgb"
    a_n = res.aerial_count
    aerial_bits = np.ascontiguousarray(data[:a_n, 3]).view(np.uint32)
    # aerial RGB preserved exactly: pure red
    assert np.all((aerial_bits >> 16) & 0xFF == 255)
    assert np.all((aerial_bits >> 8) & 0xFF == 0)
    assert np.all(aerial_bits & 0xFF == 0)
    # tube coloured by elevation -> a spread of colours, not constant
    tube_bits = np.ascontiguousarray(data[a_n:, 3]).view(np.uint32)
    assert len(np.unique(tube_bits)) > 5


def _ring(radius=5.0, n=240, center=(0.0, 0.0), z=0.0):
    th = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([center[0] + radius * np.cos(th),
                            center[1] + radius * np.sin(th), np.full(n, z)])


def _identity_tf():
    return Transform(np.eye(4).tolist(), [(0, 0)], 0.0, "4dof",
                     (0, 0, 0), (0, 0, 0), aerial_up=(0, 0, 1), anchors=[[0.0, 0.0, 0.0]])


def test_icp_refine_corrects_inplane_without_collapse(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    aerial = _ring(5.0, 240, (0.0, 0.0), 0.0)            # opening rim at z=0
    tube_ring = _ring(5.0, 240, (1.5, 0.0), 0.0)         # same rim, shifted +1.5 in x
    # deep tube body well below the opening -- must be excluded by the height band,
    # else point-to-point ICP would flatten the tube onto the ground.
    ang = rng.uniform(0, 2 * np.pi, 600); rad = rng.uniform(0, 4, 600)
    deep = np.column_stack([rad * np.cos(ang), rad * np.sin(ang),
                            rng.uniform(-15.0, -6.0, 600)])
    ap, tp = tmp_path / "a.pcd", tmp_path / "t.pcd"
    _write_pcd(ap, aerial)
    _write_pcd(tp, np.vstack([tube_ring, deep]))

    ref = icp_refine(ap, tp, _identity_tf(), radius=8.0, rim_height=3.0)
    M = ref.array
    assert ref.mode.endswith("+icp")                    # accepted
    np.testing.assert_allclose(M[:3, :3], np.eye(3), atol=0.1)    # no tilt
    np.testing.assert_allclose(M[:3, 3], [-1.5, 0.0, 0.0], atol=0.4)  # in-plane fix, no Z collapse


def test_visualize_rims_headless(tmp_path: Path) -> None:
    ap, tp = tmp_path / "a.pcd", tmp_path / "t.pcd"
    _write_pcd(ap, _ring(5.0, 240, (0.0, 0.0), 0.0))
    _write_pcd(tp, _ring(5.0, 240, (1.5, 0.0), 0.0))
    tf = _identity_tf()
    ref = icp_refine(ap, tp, tf, radius=8.0, rim_height=3.0)
    visualize_rims(ap, tp, tf, refined=ref, radius=8.0, rim_height=3.0)  # must not raise


def test_icp_refine_rejects_large_correction(tmp_path: Path) -> None:
    ap, tp = tmp_path / "a.pcd", tmp_path / "t.pcd"
    _write_pcd(ap, _ring(5.0, 240, (0.0, 0.0), 0.0))
    _write_pcd(tp, _ring(5.0, 240, (1.5, 0.0), 0.0))
    # the fit wants a 1.5 shift, but max_shift caps it -> reject, keep landmark
    ref = icp_refine(ap, tp, _identity_tf(), radius=8.0, rim_height=3.0, max_shift=0.5)
    assert ref.mode == "4dof"
    np.testing.assert_allclose(ref.array, np.eye(4), atol=1e-9)
    assert any("rejected" in w for w in ref.warnings)


def test_merge_z_offset_shifts_tube_vertically(tmp_path: Path) -> None:
    A = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    B = np.array([[2.0, 2.0, -10.0], [3.0, 3.0, -12.0]])
    ap, tp = tmp_path / "a.pcd", tmp_path / "t.pcd"
    _write_pcd(ap, A)
    _write_pcd(tp, B)
    tf = Transform(np.eye(4).tolist(), [], 0.0, "4dof", (0, 0, 0), (0, 0, 0),
                   aerial_up=(0, 0, 1))
    out = tmp_path / "m.pcd"
    res = merge_clouds(ap, tp, out, tf, color=False, z_offset=5.0, show_progress=False)
    _, data = _read_pcd_xyz(out)
    # aerial unchanged, tube lifted by +5 along z
    np.testing.assert_allclose(np.sort(data[res.aerial_count:, 2]),
                               np.sort(B[:, 2] + 5.0), atol=1e-4)


def test_vertical_residuals(tmp_path: Path) -> None:
    # aligned holes -> residuals near 0; offset one tube hole in Z -> shows up
    A = np.array([[0.0, 0.0, 0.0], [20.0, 0.0, 0.0], [10.0, 18.0, 0.0]])
    aerial = _holeset(A)
    tube = _holeset(A.copy())  # identity-aligned
    tf = match_constellations(aerial, tube, mode="6dof", tolerance=0.5)
    res = vertical_residuals(aerial, tube, tf)
    assert max(abs(r) for r in res) < 0.2


def test_merge_rejects_output_equal_input(tmp_path: Path) -> None:
    pcd = tmp_path / "a.pcd"
    _write_pcd(pcd, _ground_with_voids(_CENTERS))
    tf = match_constellations(
        _holeset([(0, 0, 0), (10, 0, 0), (0, 10, 0)]),
        _holeset([(0, 0, 0), (10, 0, 0), (0, 10, 0)]),
    )
    with pytest.raises(ValueError):
        merge_clouds(pcd, tmp_path / "b.pcd", pcd, tf, show_progress=False)
