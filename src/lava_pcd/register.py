"""Match two skylight constellations and recover the rigid transform.

Given the :class:`~lava_pcd.holes.HoleSet` from the aerial cloud and the one from
the lava-tube cloud, find the correspondence between their skylights and the
rigid transform mapping **tube-local -> aerial-local** coordinates. The clouds
barely overlap, so this landmark/constellation match is what aligns them.

It is built to tolerate a **large fraction of outlier holes** (small noisy ones,
or holes present in only one map). Correspondence-finding is decoupled from the
fit: candidate matches are gated by the **3-D inter-hole distance** (a rigid
invariant, independent of the estimated up-axis), then verified by **max
consensus** -- the transform that makes the most *other* holes line up wins.
Requiring >= ``min_inliers`` (default 3) mutually consistent skylights is what
rejects the outliers: a stray pairing can't recruit a third hole to agree.

Two fitting paths (chosen automatically by ``mode="auto"``):

* ``4dof`` -- the up-axis in each :class:`HoleSet` reduces the fit to yaw-about-up
  plus translation; needs only 2 inliers and handles a near-collinear (straight)
  tube. Candidates are hole *pairs*.
* ``6dof`` -- full Kabsch from hole *triplets*; the fallback when no up-axis is
  trustworthy. Needs >=3 non-collinear inliers.

Ellipse shape (size, axis-ratio, orientation) is used only as a **soft tie-break**
between geometrically equivalent solutions, never as a gate. The result carries an
inlier count, a uniqueness ``margin`` and a ``shape_score``, and is refined later
by a local ICP on the matched rims (:func:`lava_pcd.merge.icp_refine`).
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
    n_inliers: int = 0                        # number of matched skylights
    margin: int = 0                           # best inlier count minus next distinct solution
    shape_score: float = 0.0                  # ellipse-shape penalty over inliers (lower = better)
    aerial_up: tuple[float, float, float] = (0.0, 0.0, 1.0)  # aerial up-axis (for constrained ICP)
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
            "n_inliers": int(self.n_inliers),
            "margin": int(self.margin),
            "shape_score": float(self.shape_score),
            "aerial_up": list(self.aerial_up),
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
            n_inliers=int(d.get("n_inliers", 0)),
            margin=int(d.get("margin", 0)),
            shape_score=float(d.get("shape_score", 0.0)),
            aerial_up=tuple(d.get("aerial_up", (0.0, 0.0, 1.0))),
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


def _pairwise(P: np.ndarray) -> np.ndarray:
    """``(n, n)`` matrix of 3-D distances (the rigid invariant)."""
    return np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)


def _axis_ratio(h) -> float:
    return float(h.semi_major / max(h.semi_minor, 1e-9))


def _shape_penalty(aerial: HoleSet, tube: HoleSet, inliers, matrix) -> float:
    """Soft ellipse disagreement over the matched holes (0 = identical shapes).

    Mean over inlier matches of relative size difference, axis-ratio difference,
    and major-axis misalignment *after* the transform. Used only to break ties
    between geometrically equivalent solutions -- never to reject a match.
    """
    z = np.array([0.0, 0.0, 1.0])
    Ra = rotation_align(np.asarray(aerial.up), z)
    Rb = rotation_align(np.asarray(tube.up), z)
    R = np.asarray(matrix)[:3, :3]
    terms: list[float] = []
    for t_idx, a_idx in inliers:
        ha, hb = aerial.holes[a_idx], tube.holes[t_idx]
        sa, sb = max(ha.area, 1e-9), max(hb.area, 1e-9)
        size_t = abs(sb - sa) / (sb + sa)
        ara, arb = _axis_ratio(ha), _axis_ratio(hb)
        ratio_t = abs(arb - ara) / (arb + ara) if (arb + ara) > 1e-9 else 0.0
        # major-axis directions in each cloud's local frame, tube one transformed
        va = Ra.T @ np.array([np.cos(ha.orientation), np.sin(ha.orientation), 0.0])
        vb = R @ (Rb.T @ np.array([np.cos(hb.orientation), np.sin(hb.orientation), 0.0]))
        na, nb = np.linalg.norm(va), np.linalg.norm(vb)
        if na < 1e-9 or nb < 1e-9:
            orient_t = 0.0
        else:
            cosang = abs(float(np.dot(va, vb)) / (na * nb))  # |cos| -> undirected axes
            orient_t = float(np.arccos(np.clip(cosang, 0.0, 1.0))) / (np.pi / 2)
        terms.append((size_t + ratio_t + orient_t) / 3.0)
    return float(np.mean(terms)) if terms else 0.0


# --------------------------------------------------------------------------- #
# candidate generation (minimal-sample hypotheses, gated by 3-D distance)
# --------------------------------------------------------------------------- #
def _candidates_4dof(aerial: HoleSet, tube: HoleSet, tol: float):
    """Up-constrained *pair* hypotheses. Yields ``(matrix, inliers)``.

    The candidate pair is gated by the 3-D inter-hole distance (up-independent);
    the transform itself aligns the up-axes, then solves yaw + translation.
    """
    a_up = normalize(np.asarray(aerial.up))
    e1, e2 = _plane_basis(a_up)
    R_up = rotation_align(np.asarray(tube.up), a_up)
    A = aerial.centroids()
    B = tube.centroids()
    B3 = B @ R_up.T  # tube holes with up aligned to aerial up
    A2 = np.column_stack([A @ e1, A @ e2]); Ah = A @ a_up
    B2 = np.column_stack([B3 @ e1, B3 @ e2]); Bh = B3 @ a_up
    DA, DB = _pairwise(A), _pairwise(B)

    out, seen = [], set()
    for i0, i1 in permutations(range(len(A)), 2):
        for p0, p1 in permutations(range(len(B)), 2):
            if abs(DA[i0, i1] - DB[p0, p1]) > tol:
                continue
            dA2, dB2 = A2[i1] - A2[i0], B2[p1] - B2[p0]
            ang = np.arctan2(dA2[1], dA2[0]) - np.arctan2(dB2[1], dB2[0])
            c, s = np.cos(ang), np.sin(ang)
            R2 = np.array([[c, -s], [s, c]])
            t2 = A2[i0] - R2 @ B2[p0]
            dh = float(Ah[i0] - Bh[p0])
            R = rotation_about_axis(a_up, ang) @ R_up
            t = e1 * t2[0] + e2 * t2[1] + a_up * dh
            M = as_matrix(R, t)
            inl = _assign_inliers(transform_points(B, M), A, tol)
            if len(inl) >= 2:
                key = frozenset(inl)
                if key not in seen:
                    seen.add(key); out.append((M, inl))
    return out


def _candidates_6dof(aerial: HoleSet, tube: HoleSet, tol: float):
    """Full Kabsch *triplet* hypotheses. Yields ``(matrix, inliers)``."""
    A = aerial.centroids()
    B = tube.centroids()
    if max(len(A), len(B)) > _MAX_HOLES_6DOF:
        raise ValueError(
            f"6-DOF matching enumerates triplets and is impractical beyond "
            f"{_MAX_HOLES_6DOF} holes; prune detections or use --mode 4dof."
        )
    DA, DB = _pairwise(A), _pairwise(B)
    out, seen = [], set()
    for p0, p1, p2 in combinations(range(len(B)), 3):
        db = (DB[p0, p1], DB[p0, p2], DB[p1, p2])
        for i0, i1, i2 in permutations(range(len(A)), 3):
            if (abs(DA[i0, i1] - db[0]) > tol or abs(DA[i0, i2] - db[1]) > tol
                    or abs(DA[i1, i2] - db[2]) > tol):
                continue
            R, t = kabsch(B[[p0, p1, p2]], A[[i0, i1, i2]])
            M = as_matrix(R, t)
            inl = _assign_inliers(transform_points(B, M), A, tol)
            if len(inl) >= 3:
                key = frozenset(inl)
                if key not in seen:
                    seen.add(key); out.append((M, inl))
    return out


def _score_candidates(aerial: HoleSet, tube: HoleSet, cands):
    """Dedupe candidates by inlier set and rank by (inliers desc, RMS asc, shape asc).

    Returns a list of ``(count, rms, shape, matrix, inliers)`` groups, best first.
    """
    A, B = aerial.centroids(), tube.centroids()
    by_set: dict[frozenset, tuple] = {}
    for M, inl in cands:
        ti = [t for t, _ in inl]; ai = [a for _, a in inl]
        pred = transform_points(B[ti], M)
        rms = float(np.sqrt(np.mean(np.sum((pred - A[ai]) ** 2, axis=1))))
        sp = _shape_penalty(aerial, tube, inl, M)
        key = frozenset(inl)
        rec = (len(inl), rms, sp, M, inl)
        if key not in by_set or (rms, sp) < (by_set[key][1], by_set[key][2]):
            by_set[key] = rec
    return sorted(by_set.values(), key=lambda r: (-r[0], r[1], r[2]))


def _refit(aerial: HoleSet, tube: HoleSet, inliers, mode: str) -> np.ndarray:
    """Re-estimate the transform from all inlier correspondences."""
    A, B = aerial.centroids(), tube.centroids()
    ti = [t for t, _ in inliers]; ai = [a for _, a in inliers]
    if mode == "6dof":
        R, t = kabsch(B[ti], A[ai])
        return as_matrix(R, t)
    # 4-DOF up-constrained refit
    a_up = normalize(np.asarray(aerial.up))
    e1, e2 = _plane_basis(a_up)
    R_up = rotation_align(np.asarray(tube.up), a_up)
    B3 = B @ R_up.T
    A2 = np.column_stack([A @ e1, A @ e2]); Ah = A @ a_up
    B2 = np.column_stack([B3 @ e1, B3 @ e2]); Bh = B3 @ a_up
    ang, t2 = _kabsch2d(B2[ti], A2[ai])
    dh = float(np.mean(Ah[ai] - Bh[ti]))
    R = rotation_about_axis(a_up, ang) @ R_up
    t = e1 * t2[0] + e2 * t2[1] + a_up * dh
    return as_matrix(R, t)


# --------------------------------------------------------------------------- #
# public API
# --------------------------------------------------------------------------- #
def match_constellations(
    aerial: HoleSet,
    tube: HoleSet,
    mode: str = "auto",
    tolerance: float = DEFAULT_TOLERANCE,
    min_inliers: int = 3,
) -> Transform:
    """Match ``tube`` skylights to ``aerial`` skylights; return the rigid transform.

    Robust to a large fraction of outlier holes: correspondences are found by a
    deterministic max-consensus search over hole pairs (``4dof``) or triplets
    (``6dof``) gated by the 3-D inter-hole distance, then verified by how many
    *other* holes line up. ``min_inliers`` (default 3) mutually consistent
    skylights are required -- a stray pairing can't recruit a third, which is what
    rejects the outliers. Ellipse shape only breaks ties between equally-good
    geometric solutions.

    ``mode`` is ``"auto"`` (try 4-DOF, fall back to 6-DOF; default), ``"4dof"``
    (up-assisted; also handles a near-collinear/straight tube) or ``"6dof"`` (full
    Kabsch). ``tolerance`` is the max landmark mismatch (coordinate units). The
    returned :class:`Transform` maps tube-local -> aerial-local and carries the
    inlier count, a uniqueness ``margin`` and the ellipse ``shape_score``.
    """
    if mode not in ("auto", "4dof", "6dof"):
        raise ValueError(f"mode must be 'auto', '4dof' or '6dof', got {mode!r}")
    if len(aerial.holes) < 2 or len(tube.holes) < 2:
        raise ValueError(
            f"need >=2 holes in each cloud to register "
            f"(aerial={len(aerial.holes)}, tube={len(tube.holes)})"
        )

    if mode in ("auto", "4dof"):
        groups4 = _score_candidates(aerial, tube, _candidates_4dof(aerial, tube, tolerance))
    else:
        groups4 = []
    if mode in ("auto", "6dof"):
        try:
            groups6 = _score_candidates(aerial, tube, _candidates_6dof(aerial, tube, tolerance))
        except ValueError:
            if mode == "6dof":
                raise
            groups6 = []  # too many holes for triplets; 4-DOF still applies in auto
    else:
        groups6 = []

    best4 = groups4[0][0] if groups4 else 0
    best6 = groups6[0][0] if groups6 else 0
    if best4 >= best6:
        groups, chosen = groups4, "4dof"
    else:
        groups, chosen = groups6, "6dof"

    if not groups:
        raise ValueError(
            "no consistent correspondence found: no two skylights have matching "
            "inter-hole distances. Loosen --tolerance or recheck detections."
        )
    best = groups[0]
    best_count = best[0]
    if best_count < min_inliers:
        raise ValueError(
            f"only {best_count} mutually-consistent skylight(s) found "
            f"(need >= {min_inliers}). Check detections (use --show/--manual), "
            f"loosen --tolerance, or lower --min-inliers if you trust fewer."
        )

    inliers = best[4]
    matrix = _refit(aerial, tube, inliers, chosen)
    A, B = aerial.centroids(), tube.centroids()
    ti = [t for t, _ in inliers]; ai = [a for _, a in inliers]
    pred = transform_points(B[ti], matrix)
    rms = float(np.sqrt(np.mean(np.sum((pred - A[ai]) ** 2, axis=1))))
    shape_score = _shape_penalty(aerial, tube, inliers, matrix)
    second_count = groups[1][0] if len(groups) > 1 else 0
    margin = int(best_count - second_count)

    warnings: list[str] = []
    if margin <= 0:
        warnings.append(
            f"ambiguous: a different correspondence set explains as many skylights "
            f"({best_count}); verify with --refine and visual inspection."
        )
    if best_count == 2:
        warnings.append(
            "only 2 matched skylights: alignment has no redundancy; verify with "
            "--refine and visual inspection."
        )
    if _collinear(A[ai], tolerance):
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
        mode=chosen,
        aerial_origin=tuple(map(float, aerial.origin)),
        tube_origin=tuple(map(float, tube.origin)),
        n_inliers=int(best_count),
        margin=margin,
        shape_score=float(shape_score),
        aerial_up=tuple(map(float, aerial.up)),
        anchors=[list(map(float, A[a])) for a in ai],
        warnings=warnings,
    )


# --------------------------------------------------------------------------- #
# visualization
# --------------------------------------------------------------------------- #
def _major_axis_3d(hole, R_to_z: np.ndarray) -> np.ndarray:
    """Hole major-axis unit vector in the cloud's local frame."""
    v = np.array([np.cos(hole.orientation), np.sin(hole.orientation), 0.0])
    return R_to_z.T @ v


def _draw_ellipses(ax, holes, xs, ys, axis_dirs, matched, color, e1, e2) -> None:
    from matplotlib.patches import Ellipse
    for i, h in enumerate(holes):
        ang = np.degrees(np.arctan2(axis_dirs[i] @ e2, axis_dirs[i] @ e1))
        a = h.semi_major if h.semi_major > 0 else h.radius
        b = h.semi_minor if h.semi_minor > 0 else h.radius
        m = i in matched
        ax.add_patch(Ellipse(
            (xs[i], ys[i]), 2 * max(a, 0.5), 2 * max(b, 0.5), angle=ang,
            fill=False, color=color, lw=2.0 if m else 1.0,
            ls="-" if m else "--", alpha=1.0 if m else 0.4,
        ))
        ax.text(xs[i], ys[i], str(h.id), color=color, fontsize=7,
                ha="center", va="center")


def visualize_match(
    aerial: HoleSet, tube: HoleSet, transform: Transform, title: str | None = None
) -> None:
    """Plot the constellations **before and after** registration, side by side.

    Both panels look down the aerial up-axis (same projection, so they're directly
    comparable). **Left = before**: the tube holes in their raw frame, generally
    offset/rotated from the aerial ones (long green correspondence lines). **Right
    = after**: the tube holes transformed into the aerial frame — matched ellipses
    snap on top of each other (green lines collapse to ~0). In both: **aerial**
    holes blue, **tube** holes orange, **matched** (inlier) holes solid and joined
    by green lines, rejected outliers faded/dashed.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    z = np.array([0.0, 0.0, 1.0])
    e1, e2 = _plane_basis(normalize(np.asarray(aerial.up)))
    Ra = rotation_align(np.asarray(aerial.up), z)
    Rb = rotation_align(np.asarray(tube.up), z)
    M = transform.array
    Rm = M[:3, :3]

    A3 = aerial.centroids()
    Ax, Ay = A3 @ e1, A3 @ e2
    A_axis = [_major_axis_3d(h, Ra) for h in aerial.holes]
    in_aer = {a for _, a in transform.inliers}
    in_tube = {t for t, _ in transform.inliers}

    # tube before (raw, tube-local) and after (mapped into the aerial frame)
    B_before = tube.centroids()
    B_after = transform_points(B_before, M)
    axis_before = [_major_axis_3d(h, Rb) for h in tube.holes]
    axis_after = [Rm @ v for v in axis_before]

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    def panel(ax, B3, B_axis, subtitle):
        Bx, By = B3 @ e1, B3 @ e2
        _draw_ellipses(ax, aerial.holes, Ax, Ay, A_axis, in_aer, "tab:blue", e1, e2)
        _draw_ellipses(ax, tube.holes, Bx, By, B_axis, in_tube, "tab:orange", e1, e2)
        for t, a in transform.inliers:
            ax.plot([Ax[a], Bx[t]], [Ay[a], By[t]], "-", color="green", lw=1.0, zorder=1)
        allx = np.concatenate([Ax, Bx]); ally = np.concatenate([Ay, By])
        pad = 10.0
        ax.set_xlim(allx.min() - pad, allx.max() + pad)
        ax.set_ylim(ally.min() - pad, ally.max() + pad)
        ax.set_aspect("equal")
        ax.set_title(subtitle)
        ax.set_xlabel("aerial up-plane a"); ax.set_ylabel("aerial up-plane b")

    panel(axes[0], B_before, axis_before, "before registration")
    panel(axes[1], B_after, axis_after, f"after registration  (RMS {transform.rms:.2f})")
    fig.suptitle(title or (
        f"register: {transform.n_inliers} matches   margin {transform.margin}   "
        f"mode {transform.mode}"
    ))
    fig.legend(handles=[
        Line2D([0], [0], color="tab:blue", lw=2, label="aerial holes"),
        Line2D([0], [0], color="tab:orange", lw=2, label="tube holes"),
        Line2D([0], [0], color="green", lw=1.2, label="matches"),
    ], loc="lower center", ncol=3)
    plt.show()
