"""Streaming voxel-grid downsampling.

Open3D's ``voxel_down_sample`` needs the whole cloud in memory at once, which is
impractical for hundred-million-point clouds. This accumulator instead consumes
``(k, 4)`` chunks (``x y z intensity``) incrementally: each chunk is reduced to
its occupied voxels, buffered, and periodically compacted so memory stays
bounded by the *downsampled* size (plus a flush buffer) rather than the input.

Each output point is the centroid of the points in its voxel, with intensity
averaged the same way -- matching the usual voxel-downsample semantics.
"""

from __future__ import annotations

import numpy as np

DEFAULT_FLUSH_POINTS = 20_000_000


def _reduce(vox: np.ndarray, acc: np.ndarray, cnt: np.ndarray):
    """Sum ``acc``/``cnt`` over identical voxel rows. Returns (vox, acc, cnt)."""
    uniq, inv = np.unique(vox, axis=0, return_inverse=True)
    inv = inv.ravel()
    m = len(uniq)
    acc_out = np.zeros((m, acc.shape[1]), dtype=np.float64)
    np.add.at(acc_out, inv, acc)
    cnt_out = np.bincount(inv, weights=cnt, minlength=m)
    return uniq, acc_out, cnt_out


class VoxelDownsampler:
    """Accumulate ``(k, 4)`` chunks into voxel centroids.

    ``voxel_size`` is the edge length of each cubic voxel (same units as the
    coordinates). ``flush_points`` caps how many buffered voxel rows accumulate
    before an intermediate compaction.
    """

    def __init__(self, voxel_size: float, flush_points: int = DEFAULT_FLUSH_POINTS) -> None:
        if voxel_size <= 0:
            raise ValueError(f"voxel_size must be positive, got {voxel_size}")
        self.voxel_size = float(voxel_size)
        self.flush_points = flush_points
        self._vox: list[np.ndarray] = []
        self._acc: list[np.ndarray] = []  # summed x, y, z, intensity (float64)
        self._cnt: list[np.ndarray] = []  # point count per voxel (float64)
        self._buffered = 0

    def add(self, points: np.ndarray) -> None:
        """Add a ``(k, 4)`` chunk of ``x y z intensity``."""
        if points.shape[0] == 0:
            return
        pts = points.astype(np.float64, copy=False)
        vidx = np.floor(pts[:, :3] / self.voxel_size).astype(np.int64)
        uniq, inv = np.unique(vidx, axis=0, return_inverse=True)
        inv = inv.ravel()
        m = len(uniq)
        acc = np.zeros((m, 4), dtype=np.float64)
        np.add.at(acc, inv, pts)
        cnt = np.bincount(inv, minlength=m).astype(np.float64)

        self._vox.append(uniq)
        self._acc.append(acc)
        self._cnt.append(cnt)
        self._buffered += m
        if self._buffered >= self.flush_points:
            self._compact()

    def _compact(self) -> None:
        if len(self._vox) <= 1:
            return
        vox, acc, cnt = _reduce(
            np.concatenate(self._vox),
            np.concatenate(self._acc),
            np.concatenate(self._cnt),
        )
        self._vox, self._acc, self._cnt = [vox], [acc], [cnt]
        self._buffered = len(vox)

    def result(self) -> np.ndarray:
        """Return the downsampled cloud as a ``(N, 4)`` float32 array."""
        if not self._vox:
            return np.empty((0, 4), dtype=np.float32)
        self._compact()
        acc, cnt = self._acc[0], self._cnt[0]
        out = (acc / cnt[:, None]).astype(np.float32)
        return out
