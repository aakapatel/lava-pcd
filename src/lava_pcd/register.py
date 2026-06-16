"""Match two skylight constellations and recover the rigid transform.

Given the :class:`~lava_pcd.holes.HoleSet` from the aerial cloud and the one from
the lava-tube cloud, find the correspondence between their skylights and the
rigid transform mapping **tube-local -> aerial-local** coordinates. The clouds
barely overlap, so this landmark/constellation match is what aligns them.

Two strategies:

* ``4dof`` (default, robust) -- the up-axis recorded in each :class:`HoleSet`
  reduces the problem to yaw-about-up plus translation. Correspondences are found
  by a brute-force search over hole-pair to hole-pair matches in the horizontal
  plane (counts are tiny), scored by how many other holes then line up. Needs >=2
  matched holes and is not degenerate for collinear skylights.
* ``6dof`` (fallback) -- no reliable up-axis: a RANSAC over correspondence
  triplets (consistent pairwise distances) with a Kabsch fit. Needs >=3
  non-collinear matched holes.

The result is refined later by a local ICP on the matched rims
(:func:`lava_pcd.merge.icp_refine`).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

from lava_pcd.geometry import (
    as_matrix,
    kabsch,
    normalize,
    rotation_about_axis,
    rotation_align,
    transform_points,
)
from lava_pcd.holes import HoleSet

DEFAULT_TOLERANCE = 2.0   # max landmark mismatch (coordinate units) to count an inlier
_MAX_HOLES_6DOF = 25      # triplet RANSAC blows up beyond this; suggest 4dof/manual


@dataclass
class Transform:
    """Rigid transform (tube-local -> aerial-local) plus its support."""

    matrix: list[list[float]]                 # 4x4 homogeneous
    inliers: list[tuple[int, int]]            # (tube_id, aerial_id) correspondences
    rms: float
    mode: str
    aerial_origin: tuple[float, float, float]
    tube_origin: tuple[float, float, float]
    anchors: list[list[float]] = field(default_factory=list)  # matched rim centres (aerial frame)
    warnings: list[str] = field(default_factory=list)

    @property
    def array(self) -> np.ndarray:
        return np.asarray(self.matrix, dtype=np.float64)

    def to_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "matrix": [list(map(float, row)) for row in self.matrix],
            "inliers": [list(map(int, c)) for c in self.inliers],
            "rms": float(self.rms),
            "mode": self.mode,
            "aerial_origin": list(self.aerial_origin),
            "tube_origin": list(self.tube_origin),
            "anchors": [list(map(float, a)) for a in self.anchors],
            "warnings": list(self.warnings),
        }
        path.write_text(json.dumps(payload, indent=2))
        return path

    @classmethod
    def from_json(cls, path: str | Path) -> "Transform":
        d = json.loads(Path(path).read_text())
        return cls(
            matrix=[list(map(float, row)) for row in d["matrix"]],
            inliers=[tuple(map(int, c)) for c in d["inliers"]],
            rms=float(d["rms"]),
            mode=d["mode"],
            aerial_origin=tuple(d["aerial_origin"]),
            tube_origin=tuple(d["tube_origin"]),
            anchors=[list(map(float, a)) for a in d.get("anchors", [])],
            warnings=list(d.get("warnings", [])),
        )


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _plane_basis(up: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Right-handed in-plane basis ``(e1, e2)`` with ``e1 x e2 = up``."""
    up = normalize(up)
    seed = np.array([1.0, 0.0, 0.0]) if abs(up[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = normalize(seed - up * np.dot(seed, up))
    e2 = np.cross(up, e1)
    return e1, e2


def _kabsch2d(src: np.ndarray, dst: np.ndarray) -> tuple[float, np.ndarray]:
    """2-D rigid fit ``dst ~ R(angle) @ src + t``. Returns ``(angle, t)``."""
    cs, cd = src.mean(axis=0), dst.mean(axis=0)
    H = (src - cs).T @ (dst - cd)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, d]) @ U.T
    angle = float(np.arctan2(R[1, 0], R[0, 0]))
    t = cd - R @ cs
    return angle, t


def _assign_inliers(pred: np.ndarray, target: np.ndarray, tol: float) -> list[tuple[int, int]]:
    """Greedy nearest-neighbour assignment of ``pred`` rows to ``target`` rows.

    Returns one-to-one ``(pred_idx, target_idx)`` pairs within ``tol``.
    """
    pairs: list[tuple[int, int, float]] = []
    for i in range(len(pred)):
        d = np.linalg.norm(target - pred[i], axis=1)
        j = int(np.argmin(d))
        if d[j] <= tol:
            pairs.append((i, j, float(d[j])))
    pairs.sort(key=lambda p: p[2])
    used_p, used_t, out = set(), set(), []
    for i, j, _ in pairs:
        if i in used_p or j in used_t:
            continue
        used_p.add(i); used_t.add(j); out.append((i, j))
    return out


def _collinear(points: np.ndarray, tol: float) -> bool:
    """True if all points lie ~on a line (degenerate for a 3-D rotation fit)."""
    if len(points) < 3:
        return True
    c = points - points.mean(axis=0)
    s = np.linalg.svd(c, compute_uv=False)
    return s[1] < tol  # second singular value tiny -> rank 1


# --------------------------------------------------------------------------- #
# 4-DOF (up-assisted) matching
# --------------------------------------------------------------------------- #
def _match_4dof(aerial: HoleSet, tube: HoleSet, tol: float):
    a_up = normalize(np.asarray(aerial.up))
    e1, e2 = _plane_basis(a_up)
    R_up = rotation_align(np.asarray(tube.up), a_up)

    A = aerial.centroids()
    B = tube.centroids() @ R_up.T  # tube holes with up aligned to aerial up
    A2 = np.column_stack([A @ e1, A @ e2])
    B2 = np.column_stack([B @ e1, B @ e2])
    Ah, Bh = A @ a_up, B @ a_up

    best_inliers: list[tuple[int, int]] = []
    # Each ordered aerial pair vs ordered tube pair fixes a 2-D rotation+translation.
    for i0, i1 in permutations(range(len(A2)), 2):
        dA = A2[i1] - A2[i0]
        lA = np.hypot(*dA)
        angA = np.arctan2(dA[1], dA[0])
        for p0, p1 in permutations(range(len(B2)), 2):
            dB = B2[p1] - B2[p0]
            if abs(np.hypot(*dB) - lA) > tol:
                continue
            ang = angA - np.arctan2(dB[1], dB[0])
            c, s = np.cos(ang), np.sin(ang)
            R2 = np.array([[c, -s], [s, c]])
            t2 = A2[i0] - R2 @ B2[p0]
            pred = B2 @ R2.T + t2
            inl = _assign_inliers(pred, A2, tol)
            if len(inl) > len(best_inliers):
                best_inliers = inl

    if len(best_inliers) < 2:
        raise ValueError(
            "4-DOF match failed: fewer than 2 skylights line up. Check detections "
            "(use --show/--manual), loosen --tolerance, or supply a better --up."
        )

    tube_idx = [p for p, _ in best_inliers]
    aer_idx = [a for _, a in best_inliers]
    ang, t2 = _kabsch2d(B2[tube_idx], A2[aer_idx])
    dh = float(np.mean(Ah[aer_idx] - Bh[tube_idx]))
    R = rotation_about_axis(a_up, ang) @ R_up
    t = e1 * t2[0] + e2 * t2[1] + a_up * dh
    return as_matrix(R, t), best_inliers


# --------------------------------------------------------------------------- #
# 6-DOF (Kabsch RANSAC) matching
# --------------------------------------------------------------------------- #
def _match_6dof(aerial: HoleSet, tube: HoleSet, tol: float):
    A = aerial.centroids()
    B = tube.centroids()
    if max(len(A), len(B)) > _MAX_HOLES_6DOF:
        raise ValueError(
            f"6-DOF matching enumerates triplets and is impractical beyond "
            f"{_MAX_HOLES_6DOF} holes; prune detections or use --mode 4dof."
        )

    def pdist(p, i, j):
        return float(np.linalg.norm(p[i] - p[j]))

    best_inliers: list[tuple[int, int]] = []
    for p0, p1, p2 in combinations(range(len(B)), 3):
        db = (pdist(B, p0, p1), pdist(B, p0, p2), pdist(B, p1, p2))
        for i0, i1, i2 in permutations(range(len(A)), 3):
            if (abs(pdist(A, i0, i1) - db[0]) > tol
                    or abs(pdist(A, i0, i2) - db[1]) > tol
                    or abs(pdist(A, i1, i2) - db[2]) > tol):
                continue
            R, t = kabsch(B[[p0, p1, p2]], A[[i0, i1, i2]])
            pred = B @ R.T + t
            inl = _assign_inliers(pred, A, tol)
            if len(inl) > len(best_inliers):
                best_inliers = inl

    if len(best_inliers) < 3:
        raise ValueError(
            "6-DOF match failed: fewer than 3 skylights line up. Need >=3 "
            "non-collinear matches; check detections or try --mode 4dof."
        )
    tube_idx = [p for p, _ in best_inliers]
    aer_idx = [a for _, a in best_inliers]
    R, t = kabsch(B[tube_idx], A[aer_idx])
    return as_matrix(R, t), best_inliers


# --------------------------------------------------------------------------- #
# public API
# --------------------------------------------------------------------------- #
def match_constellations(
    aerial: HoleSet,
    tube: HoleSet,
    mode: str = "4dof",
    tolerance: float = DEFAULT_TOLERANCE,
) -> Transform:
    """Match ``tube`` skylights to ``aerial`` skylights; return the rigid transform.

    ``mode`` is ``"4dof"`` (up-assisted; default, robust) or ``"6dof"`` (full
    Kabsch-RANSAC fallback). ``tolerance`` is the max landmark mismatch (in
    coordinate units) for a hole to count as an inlier. The returned
    :class:`Transform` maps tube-local -> aerial-local coordinates.
    """
    if mode not in ("4dof", "6dof"):
        raise ValueError(f"mode must be '4dof' or '6dof', got {mode!r}")
    if len(aerial.holes) < 2 or len(tube.holes) < 2:
        raise ValueError(
            f"need >=2 holes in each cloud to register "
            f"(aerial={len(aerial.holes)}, tube={len(tube.holes)})"
        )

    if mode == "4dof":
        matrix, inliers = _match_4dof(aerial, tube, tolerance)
    else:
        matrix, inliers = _match_6dof(aerial, tube, tolerance)

    # Residual RMS over the matched landmarks.
    A = aerial.centroids()
    B = tube.centroids()
    tube_idx = [p for p, _ in inliers]
    aer_idx = [a for _, a in inliers]
    pred = transform_points(B[tube_idx], np.asarray(matrix))
    rms = float(np.sqrt(np.mean(np.sum((pred - A[aer_idx]) ** 2, axis=1))))

    warnings: list[str] = []
    if len(inliers) == 2:
        warnings.append(
            "only 2 matched skylights: alignment has no redundancy; verify with "
            "--refine and visual inspection."
        )
    if _collinear(A[aer_idx], tolerance):
        warnings.append(
            "matched skylights are nearly collinear: rotation about that line is "
            "weakly constrained; --refine (ICP on the rims) is recommended."
        )
    if rms > tolerance:
        warnings.append(f"landmark RMS {rms:.2f} exceeds tolerance {tolerance:.2f}.")

    return Transform(
        matrix=[list(map(float, row)) for row in np.asarray(matrix)],
        inliers=[(int(p), int(a)) for p, a in inliers],
        rms=rms,
        mode=mode,
        aerial_origin=tuple(map(float, aerial.origin)),
        tube_origin=tuple(map(float, tube.origin)),
        anchors=[list(map(float, A[a])) for a in aer_idx],
        warnings=warnings,
    )
