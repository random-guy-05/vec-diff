from __future__ import annotations

from pathlib import Path

import numpy as np

from .anndata_io import expression_summary, gene_alignment, load_anndata_backed, paired_hashes, read_coords
from .core import classify, multiset_overlap, sha256_file, spatial_stats


def _key_diff(a, b):
    sa, sb = set(a), set(b)
    return sorted(sa - sb), sorted(sb - sa)


def compare_files(path_a: Path, path_b: Path, task: str, top_genes: int = 15, spatial_sample: int = 512) -> dict:
    path_a, path_b = path_a.resolve(), path_b.resolve()
    sha_a, sha_b = sha256_file(path_a), sha256_file(path_b)
    a = load_anndata_backed(path_a)
    b = load_anndata_backed(path_b)
    try:
        genes = gene_alignment(a, b)
        col_order_b = genes["order_b_to_a"] if genes["same_set"] else None
        expr = expression_summary(a, b, col_order_b, compute_hashes=bool(genes["same_set"]))
        expr_overlap = 0.0
        if genes["same_set"]:
            _, _, _, expr_overlap = multiset_overlap(expr["row_hashes_a"], expr["row_hashes_b"])
        expr["expression_multiset_overlap"] = expr_overlap
        pa, pb = expr.get("pseudobulk_a"), expr.get("pseudobulk_b")
        top = []
        if genes["same_set"] and pa is not None and pb is not None:
            d = np.asarray(pb) - np.asarray(pa)
            idx = np.argsort(np.abs(d))[::-1][:top_genes]
            top = [{"gene": genes["genes_a"][int(i)], "a": float(pa[i]), "b": float(pb[i]), "delta": float(d[i])} for i in idx]

        spatial = None
        paired_overlap = None
        coord_identical = None
        kabsch_rel = None
        if task in {"T2", "T3"}:
            ca, cb = read_coords(a), read_coords(b)
            spatial = {"present_a": ca is not None, "present_b": cb is not None}
            if ca is not None and cb is not None:
                coord_identical = bool(ca.shape == cb.shape and np.array_equal(ca, cb))
                spatial["rowwise_identical"] = coord_identical
                if genes["same_set"]:
                    ha = paired_hashes(a, None, ca)
                    hb = paired_hashes(b, col_order_b, cb)
                    _, _, _, paired_overlap = multiset_overlap(ha, hb)
                    spatial["paired_multiset_overlap"] = paired_overlap
                stats = spatial_stats(ca, cb, sample_size=spatial_sample)
                spatial.update(stats)
                kabsch_rel = stats.get("kabsch_relative_rmsd")
            else:
                spatial["paired_multiset_overlap"] = None

        obs_a, obs_b = _key_diff(a.obs.columns, b.obs.columns)
        layers_a, layers_b = _key_diff(a.layers.keys(), b.layers.keys())
        obsm_a, obsm_b = _key_diff(a.obsm.keys(), b.obsm.keys())
        uns_a, uns_b = _key_diff(a.uns.keys(), b.uns.keys())
        var_a, var_b = _key_diff(a.var.columns, b.var.columns)
        metadata = {
            "obs_only_a": obs_a, "obs_only_b": obs_b,
            "layers_only_a": layers_a, "layers_only_b": layers_b,
            "obsm_only_a": obsm_a, "obsm_only_b": obsm_b,
            "uns_only_a": uns_a, "uns_only_b": uns_b,
            "var_columns_only_a": var_a, "var_columns_only_b": var_b,
        }
        if task in {"T2", "T3"} and spatial is not None and (not spatial["present_a"] or not spatial["present_b"]):
            classification = "INCOMPATIBLE"
            interpretation = "At least one file is missing a finite cells×>=3 obsm['spatial_3D']; the current spatial-task scorer contract cannot compare it as a valid T2/T3 prediction."
        else:
            classification, interpretation = classify(
                byte_identical=sha_a == sha_b,
                genes_compatible=bool(genes["same_set"]),
                task=task,
                expression_rowwise_identical=bool(expr.get("rowwise_identical", False)),
                expression_multiset_overlap=expr_overlap,
                paired_multiset_overlap=paired_overlap,
                kabsch_relative_rmsd=float(kabsch_rel) if kabsch_rel is not None else None,
                coord_rowwise_identical=coord_identical,
            )
        # Drop bulky internal arrays/hashes from report.
        for k in ["pseudobulk_a", "pseudobulk_b", "row_hashes_a", "row_hashes_b"]:
            expr.pop(k, None)
        genes.pop("order_b_to_a", None)
        return {
            "schema_version": 1,
            "task": task,
            "classification": classification,
            "interpretation": interpretation,
            "files": {
                "a": {"path": str(path_a), "sha256": sha_a, "size_bytes": path_a.stat().st_size, "shape": [int(a.n_obs), int(a.n_vars)]},
                "b": {"path": str(path_b), "sha256": sha_b, "size_bytes": path_b.stat().st_size, "shape": [int(b.n_obs), int(b.n_vars)]},
            },
            "genes": genes,
            "expression": expr,
            "top_genes": top,
            "spatial": spatial,
            "metadata": metadata,
        }
    finally:
        try: a.file.close()
        except Exception: pass
        try: b.file.close()
        except Exception: pass
