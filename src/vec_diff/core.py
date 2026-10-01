from __future__ import annotations

import hashlib
import math
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.spatial.distance import pdist
from scipy.stats import wasserstein_distance


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_float32(x) -> np.ndarray:
    arr = np.asarray(x, dtype=np.float32)
    return np.ascontiguousarray(arr)


def row_hashes_dense(x: np.ndarray, coords: np.ndarray | None = None) -> list[str]:
    x = canonical_float32(x)
    c = canonical_float32(coords) if coords is not None else None
    if c is not None and len(c) != len(x):
        raise ValueError("coords must have the same number of rows as X")
    out = []
    for i in range(len(x)):
        h = hashlib.blake2b(digest_size=16)
        h.update(x[i].tobytes(order="C"))
        if c is not None:
            h.update(c[i, :3].tobytes(order="C"))
        out.append(h.hexdigest())
    return out


def multiset_overlap(a: Iterable[str], b: Iterable[str]) -> tuple[int, int, int, float]:
    ca, cb = Counter(a), Counter(b)
    common = sum((ca & cb).values())
    na, nb = sum(ca.values()), sum(cb.values())
    denom = max(na, nb, 1)
    return common, na, nb, common / denom


def expression_differences(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    aa, bb = canonical_float32(a), canonical_float32(b)
    if aa.shape != bb.shape:
        raise ValueError("expression arrays must have the same shape")
    d = aa.astype(np.float64) - bb.astype(np.float64)
    ad = np.abs(d)
    return {
        "mean_abs_diff": float(ad.mean()) if ad.size else 0.0,
        "rmse": float(np.sqrt(np.mean(d * d))) if d.size else 0.0,
        "max_abs_diff": float(ad.max()) if ad.size else 0.0,
        "exact_fraction": float(np.mean(aa == bb)) if aa.size else 1.0,
    }


def pseudobulk(x: np.ndarray) -> np.ndarray:
    xx = canonical_float32(x)
    if xx.shape[0] == 0:
        return np.full(xx.shape[1], np.nan)
    return xx.astype(np.float64).mean(axis=0)


def pseudobulk_stats(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    pa, pb = pseudobulk(a), pseudobulk(b)
    d = pb - pa
    finite = np.isfinite(pa) & np.isfinite(pb)
    if finite.sum() >= 2 and np.std(pa[finite]) > 0 and np.std(pb[finite]) > 0:
        corr = float(np.corrcoef(pa[finite], pb[finite])[0, 1])
    elif np.array_equal(pa[finite], pb[finite]):
        corr = 1.0
    else:
        corr = float("nan")
    return {
        "pearson": corr,
        "mean_abs_delta": float(np.nanmean(np.abs(d))) if len(d) else 0.0,
        "max_abs_delta": float(np.nanmax(np.abs(d))) if len(d) else 0.0,
    }


def rms_radius(coords: np.ndarray) -> float:
    c = canonical_float32(coords)[:, :3].astype(np.float64)
    if len(c) == 0:
        return float("nan")
    centered = c - c.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum(centered * centered, axis=1))))


def kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    """Best proper-rotation + translation RMSD mapping b onto a, same row correspondence."""
    aa = canonical_float32(a)[:, :3].astype(np.float64)
    bb = canonical_float32(b)[:, :3].astype(np.float64)
    if aa.shape != bb.shape:
        raise ValueError("coordinate arrays must have equal shape")
    if len(aa) == 0:
        return {"rmsd": 0.0, "relative_rmsd": 0.0, "det_rotation": 1.0}
    ac = aa - aa.mean(axis=0)
    bc = bb - bb.mean(axis=0)
    h = bc.T @ ac
    u, _, vt = np.linalg.svd(h)
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    aligned = bc @ r
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - ac) ** 2, axis=1))))
    scale = max(rms_radius(aa), 1e-12)
    return {"rmsd": rmsd, "relative_rmsd": rmsd / scale, "det_rotation": float(np.linalg.det(r))}


def spatial_stats(a: np.ndarray, b: np.ndarray, sample_size: int = 512, seed: int = 0) -> dict[str, object]:
    aa = canonical_float32(a)[:, :3].astype(np.float64)
    bb = canonical_float32(b)[:, :3].astype(np.float64)
    ra, rb = rms_radius(aa), rms_radius(bb)
    result: dict[str, object] = {
        "centroid_a": aa.mean(axis=0).tolist() if len(aa) else [float("nan")] * 3,
        "centroid_b": bb.mean(axis=0).tolist() if len(bb) else [float("nan")] * 3,
        "centroid_shift": float(np.linalg.norm(aa.mean(axis=0) - bb.mean(axis=0))) if len(aa) and len(bb) else float("nan"),
        "rms_radius_a": ra,
        "rms_radius_b": rb,
        "log_scale_b_over_a": float(math.log(rb / ra)) if ra > 0 and rb > 0 else float("nan"),
    }
    if aa.shape == bb.shape:
        result.update({f"kabsch_{k}": v for k, v in kabsch_rmsd(aa, bb).items()})
        row_delta = np.linalg.norm(aa - bb, axis=1)
        result["same_row_coord_rmse"] = float(np.sqrt(np.mean(row_delta**2))) if len(row_delta) else 0.0
        result["same_row_coord_max"] = float(row_delta.max()) if len(row_delta) else 0.0
    def choose(c):
        if len(c) <= sample_size:
            return c
        rng = np.random.default_rng(seed)
        return c[np.sort(rng.choice(len(c), size=sample_size, replace=False))]
    sa, sb = choose(aa), choose(bb)
    da = pdist(sa) if len(sa) >= 2 else np.array([], dtype=float)
    db = pdist(sb) if len(sb) >= 2 else np.array([], dtype=float)
    if len(da) and len(db):
        wd = float(wasserstein_distance(da, db))
        result["pairwise_distance_wasserstein"] = wd
        result["pairwise_distance_wasserstein_relative"] = wd / max(float(np.median(da)), 1e-12)
    else:
        result["pairwise_distance_wasserstein"] = float("nan")
        result["pairwise_distance_wasserstein_relative"] = float("nan")
    return result


def classify(
    *,
    byte_identical: bool,
    genes_compatible: bool,
    task: str,
    expression_rowwise_identical: bool,
    expression_multiset_overlap: float,
    paired_multiset_overlap: float | None,
    kabsch_relative_rmsd: float | None,
    coord_rowwise_identical: bool | None,
    tolerance: float = 1e-6,
) -> tuple[str, str]:
    if byte_identical:
        return "BYTE_IDENTICAL", "The two files are byte-for-byte identical."
    if not genes_compatible:
        return "INCOMPATIBLE", "The files do not contain the same gene set, so scorer-aware equivalence cannot be established."
    if task == "T1":
        if expression_rowwise_identical:
            return "SCORER_CONTENT_IDENTICAL", "Expression content is identical after scorer-style float32 casting; differences are outside T1 scorer-relevant X/gene content."
        if expression_multiset_overlap >= 1.0 - 1e-12:
            return "ROW_ORDER_ONLY", "The exact predicted expression-row multiset is identical; only cell row order and/or non-scorer metadata differ. T1 metrics are row-order invariant."
        return "MATERIAL", "The predicted expression distribution differs in scorer-relevant content."
    if expression_rowwise_identical and coord_rowwise_identical:
        return "SCORER_CONTENT_IDENTICAL", "Expression and the first three spatial coordinates are identical after scorer-style float32 casting."
    if paired_multiset_overlap is not None and paired_multiset_overlap >= 1.0 - 1e-12:
        return "ROW_ORDER_ONLY", "The exact (expression, spatial_3D[:3]) paired-cell multiset is identical; only row order and/or non-scorer metadata differ."
    if expression_rowwise_identical and kabsch_relative_rmsd is not None and kabsch_relative_rmsd <= tolerance:
        return "RIGID_FRAME_ONLY", "Expression is rowwise identical and coordinates differ only by a proper rigid rotation/translation within tolerance; current VEC spatial metrics are rotation/translation invariant."
    return "MATERIAL", "Expression and/or expression-to-position coupling differs in scorer-relevant content."
