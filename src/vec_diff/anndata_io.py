from __future__ import annotations

import hashlib
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import sparse

from .core import canonical_float32


def load_anndata_backed(path: Path):
    try:
        import anndata as ad
    except ImportError as exc:
        raise RuntimeError(
            "anndata is required. Install with `pip install anndata`."
        ) from exc
    return ad.read_h5ad(path, backed="r")


def chunk_rows(
    n_vars: int,
    max_bytes: int = 64 * 1024 * 1024,
) -> int:
    return max(1, min(1024, max_bytes // max(n_vars * 4, 1)))


def dense_chunk(
    adata,
    start: int,
    stop: int,
    col_order: np.ndarray | None = None,
) -> np.ndarray:
    x = adata.X[start:stop]
    if sparse.issparse(x):
        x = x.toarray()
    x = np.asarray(x)
    if col_order is not None:
        x = x[:, col_order]
    return canonical_float32(x)


def gene_alignment(a, b) -> dict[str, object]:
    ga = [str(x) for x in a.var_names]
    gb = [str(x) for x in b.var_names]
    set_a, set_b = set(ga), set(gb)
    duplicates_a = [g for g, n in Counter(ga).items() if n > 1]
    duplicates_b = [g for g, n in Counter(gb).items() if n > 1]
    same_set = (
        set_a == set_b
        and len(ga) == len(gb)
        and not duplicates_a
        and not duplicates_b
    )
    same_order = ga == gb
    order_b = None
    if same_set:
        pos = {g: i for i, g in enumerate(gb)}
        order_b = np.asarray([pos[g] for g in ga], dtype=int)

    return {
        "genes_a": ga,
        "genes_b": gb,
        "same_set": same_set,
        "same_order": same_order,
        "duplicates_a": duplicates_a,
        "duplicates_b": duplicates_b,
        "missing_from_b": sorted(set_a - set_b)[:20],
        "extra_in_b": sorted(set_b - set_a)[:20],
        "order_b_to_a": order_b,
    }


def expression_summary(
    a,
    b,
    col_order_b: np.ndarray | None,
    compute_hashes: bool = True,
) -> dict[str, object]:
    if a.n_vars != b.n_vars and col_order_b is None:
        return {"comparable": False}

    nvars = a.n_vars
    step = chunk_rows(nvars)
    na, nb = a.n_obs, b.n_obs

    sum_a = np.zeros(nvars, dtype=np.float64)
    sum_b = np.zeros(nvars, dtype=np.float64)
    min_a = math.inf
    max_a = -math.inf
    min_b = math.inf
    max_b = -math.inf
    hashes_a: list[str] = []
    hashes_b: list[str] = []

    same_shape = na == nb and a.n_vars == b.n_vars
    abs_sum = 0.0
    sq_sum = 0.0
    max_abs = 0.0
    exact = 0
    count = 0

    def hash_rows(x):
        out = []
        for row in x:
            digest = hashlib.blake2b(
                np.ascontiguousarray(row).tobytes(),
                digest_size=16,
            )
            out.append(digest.hexdigest())
        return out

    for start in range(0, max(na, nb), step):
        if start < na:
            xa = dense_chunk(a, start, min(start + step, na))
            sum_a += xa.astype(np.float64).sum(axis=0)
            if xa.size:
                min_a = min(min_a, float(xa.min()))
                max_a = max(max_a, float(xa.max()))
            if compute_hashes:
                hashes_a.extend(hash_rows(xa))
        else:
            xa = None

        if start < nb:
            xb = dense_chunk(
                b,
                start,
                min(start + step, nb),
                col_order_b,
            )
            sum_b += xb.astype(np.float64).sum(axis=0)
            if xb.size:
                min_b = min(min_b, float(xb.min()))
                max_b = max(max_b, float(xb.max()))
            if compute_hashes:
                hashes_b.extend(hash_rows(xb))
        else:
            xb = None

        if same_shape and xa is not None and xb is not None:
            delta = xa.astype(np.float64) - xb.astype(np.float64)
            abs_delta = np.abs(delta)
            abs_sum += float(abs_delta.sum())
            sq_sum += float((delta * delta).sum())
            max_abs = max(
                max_abs,
                float(abs_delta.max()) if abs_delta.size else 0.0,
            )
            exact += int(np.count_nonzero(xa == xb))
            count += int(xa.size)

    pa = sum_a / max(na, 1)
    pb = sum_b / max(nb, 1)
    pd = pb - pa
    finite = np.isfinite(pa) & np.isfinite(pb)
    corr = float("nan")

    if finite.sum() >= 2:
        if np.std(pa[finite]) > 0 and np.std(pb[finite]) > 0:
            corr = float(np.corrcoef(pa[finite], pb[finite])[0, 1])
        elif np.array_equal(pa[finite], pb[finite]):
            corr = 1.0

    out: dict[str, object] = {
        "comparable": True,
        "n_cells_a": na,
        "n_cells_b": nb,
        "n_genes": nvars,
        "min_a": None if math.isinf(min_a) else min_a,
        "max_a": None if math.isinf(max_a) else max_a,
        "min_b": None if math.isinf(min_b) else min_b,
        "max_b": None if math.isinf(max_b) else max_b,
        "pseudobulk_a": pa,
        "pseudobulk_b": pb,
        "pseudobulk_pearson": corr,
        "pseudobulk_mean_abs_delta": (
            float(np.mean(np.abs(pd))) if len(pd) else 0.0
        ),
        "pseudobulk_max_abs_delta": (
            float(np.max(np.abs(pd))) if len(pd) else 0.0
        ),
        "row_hashes_a": hashes_a,
        "row_hashes_b": hashes_b,
    }

    if same_shape:
        out.update(
            {
                "rowwise_mean_abs_diff": abs_sum / max(count, 1),
                "rowwise_rmse": math.sqrt(sq_sum / max(count, 1)),
                "rowwise_max_abs_diff": max_abs,
                "rowwise_exact_fraction": exact / max(count, 1),
                "rowwise_identical": exact == count,
            }
        )
    else:
        out.update({"rowwise_identical": False})

    return out


def read_coords(adata) -> np.ndarray | None:
    if "spatial_3D" not in adata.obsm:
        return None
    coords = np.asarray(adata.obsm["spatial_3D"], dtype=np.float32)
    if (
        coords.ndim != 2
        or coords.shape[0] != adata.n_obs
        or coords.shape[1] < 3
    ):
        return None
    coords = np.ascontiguousarray(coords[:, :3])
    if not np.isfinite(coords).all():
        return None
    return coords


def paired_hashes(
    a,
    col_order: np.ndarray | None,
    coords: np.ndarray,
) -> list[str]:
    out = []
    step = chunk_rows(a.n_vars)
    for start in range(0, a.n_obs, step):
        stop = min(start + step, a.n_obs)
        x = dense_chunk(a, start, stop, col_order)
        chunk_coords = coords[start:stop]
        for row, coord in zip(x, chunk_coords):
            digest = hashlib.blake2b(digest_size=16)
            digest.update(np.ascontiguousarray(row).tobytes())
            digest.update(
                np.ascontiguousarray(
                    coord,
                    dtype=np.float32,
                ).tobytes()
            )
            out.append(digest.hexdigest())
    return out
