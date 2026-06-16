"""Small rigid-geometry helpers shared by the skylight-merge pipeline.

Pure numpy, no external dependency. Conventions:

* A rigid transform is a rotation ``R`` (3x3) and translation ``t`` (3,) acting
  as ``y = R @ x + t``.
* :func:`as_matrix` / :func:`from_matrix` convert to/from a 4x4 homogeneous
  matrix (the form stored in :class:`lava_pcd.register.Transform`).
"""

from __future__ import annotations

import numpy as np

_EPS = 1e-9


def normalize(v: np.ndarray) -> np.ndarray:
    """Return ``v`` scaled to unit length (unchanged if it is ~zero)."""
    v = np.asarray(v, dtype=np.float64)
    n = np.linalg.norm(v)
    return v / n if n > _EPS else v


def rotation_align(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation matrix taking unit vector ``a`` onto unit vector ``b``.

    Uses Rodrigues' formula for the shortest-arc rotation. Handles the
    anti-parallel case (180 degrees) by rotating about an arbitrary
    perpendicular axis.
    """
    a = normalize(a)
    b = normalize(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    s = float(np.linalg.norm(v))
    if s < _EPS:
        if c > 0:  # already aligned
            return np.eye(3)
        # Anti-parallel: rotate 180 deg about any axis perpendicular to a.
        axis = np.cross(a, [1.0, 0.0, 0.0])
        if np.linalg.norm(axis) < _EPS:
            axis = np.cross(a, [0.0, 1.0, 0.0])
        return rotation_about_axis(normalize(axis), np.pi)
    vx = np.array([
        [0.0, -v[2], v[1]],
        [v[2], 0.0, -v[0]],
        [-v[1], v[0], 0.0],
    ])
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))


def rotation_about_axis(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotation matrix of ``angle`` radians about unit ``axis`` (Rodrigues)."""
    axis = normalize(axis)
    x, y, z = axis
    c, s = np.cos(angle), np.sin(angle)
    C = 1.0 - c
    return np.array([
        [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def kabsch(src: np.ndarray, dst: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Best-fit rigid transform mapping ``src`` onto ``dst`` (no scaling).

    ``src`` and ``dst`` are ``(N, 3)`` corresponding point sets. Returns
    ``(R, t)`` minimising ``sum ||R @ src_i + t - dst_i||^2`` via SVD, with a
    reflection guard so ``det(R) = +1``.
    """
    src = np.asarray(src, dtype=np.float64)
    dst = np.asarray(dst, dtype=np.float64)
    cs = src.mean(axis=0)
    cd = dst.mean(axis=0)
    H = (src - cs).T @ (dst - cd)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = cd - R @ cs
    return R, t


def as_matrix(R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Pack rotation ``R`` and translation ``t`` into a 4x4 homogeneous matrix."""
    M = np.eye(4)
    M[:3, :3] = R
    M[:3, 3] = t
    return M


def from_matrix(M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split a 4x4 homogeneous matrix into ``(R, t)``."""
    M = np.asarray(M, dtype=np.float64)
    return M[:3, :3].copy(), M[:3, 3].copy()


def transform_points(xyz: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Apply a 4x4 homogeneous transform to an ``(N, 3)`` point array."""
    R, t = from_matrix(M)
    return xyz @ R.T + t
