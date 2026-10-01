from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


def _safe(obj):
    if isinstance(obj, np.ndarray):
        return _safe(obj.tolist())
    if isinstance(obj, (np.floating, np.integer)):
        return _safe(obj.item())
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    if isinstance(obj, dict):
        return {k: _safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_safe(v) for v in obj]
    return obj


def write_json(report: dict, path: Path) -> None:
    path.write_text(json.dumps(_safe(report), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _fmt(v, digits=6):
    if v is None:
        return "—"
    try:
        x = float(v)
        if not math.isfinite(x):
            return "—"
        return f"{x:.{digits}g}"
    except (TypeError, ValueError, OverflowError):
        return str(v)


def markdown(report: dict) -> str:
    f = report["files"]
    g = report["genes"]
    e = report["expression"]
    lines = [
        "# VEC Diff report",
        "",
        f"**Classification: `{report['classification']}`**",
        "",
        report["interpretation"],
        "",
        f"- Task: **{report['task']}**",
        f"- A: `{f['a']['path']}`",
        f"- B: `{f['b']['path']}`",
        f"- SHA256 A: `{f['a']['sha256']}`",
        f"- SHA256 B: `{f['b']['sha256']}`",
        "",
        "## Gene contract",
        "",
        f"- Same gene set: **{g['same_set']}**",
        f"- Same gene order: **{g['same_order']}**",
        f"- Genes: {len(g['genes_a'])} vs {len(g['genes_b'])}",
    ]
    if g["missing_from_b"]:
        lines.append(f"- Missing from B (first 20): `{g['missing_from_b']}`")
    if g["extra_in_b"]:
        lines.append(f"- Extra in B (first 20): `{g['extra_in_b']}`")
    lines.extend([
        "",
        "## Expression",
        "",
        f"- Cells: {e.get('n_cells_a', '—')} vs {e.get('n_cells_b', '—')}",
        f"- Rowwise identical after float32 cast/alignment: **{e.get('rowwise_identical', False)}**",
        f"- Exact expression-row multiset overlap: **{_fmt(e.get('expression_multiset_overlap'))}**",
        f"- Pseudobulk Pearson: {_fmt(e.get('pseudobulk_pearson'))}",
        f"- Pseudobulk mean absolute delta: {_fmt(e.get('pseudobulk_mean_abs_delta'))}",
        f"- Pseudobulk max absolute delta: {_fmt(e.get('pseudobulk_max_abs_delta'))}",
        f"- Same-row mean absolute difference: {_fmt(e.get('rowwise_mean_abs_diff'))}",
        f"- Same-row RMSE: {_fmt(e.get('rowwise_rmse'))}",
        "",
        "### Largest pseudobulk changes",
        "",
        "| gene | A | B | B - A |",
        "|---|---:|---:|---:|",
    ])
    for r in report.get("top_genes", []):
        lines.append(f"| {r['gene']} | {_fmt(r['a'])} | {_fmt(r['b'])} | {_fmt(r['delta'])} |")
    if report.get("spatial") is not None:
        s = report["spatial"]
        lines.extend([
            "",
            "## Spatial content",
            "",
            f"- `spatial_3D` present: A={s['present_a']}, B={s['present_b']}",
            f"- First-three-coordinate rowwise identical: **{s.get('rowwise_identical', False)}**",
            f"- Exact paired (expression + coordinates) multiset overlap: **{_fmt(s.get('paired_multiset_overlap'))}**",
            f"- RMS radius A / B: {_fmt(s.get('rms_radius_a'))} / {_fmt(s.get('rms_radius_b'))}",
            f"- log(scale B/A): {_fmt(s.get('log_scale_b_over_a'))}",
            f"- Centroid shift: {_fmt(s.get('centroid_shift'))} *(frame-dependent; current VEC spatial metrics are translation-invariant)*",
            f"- Proper-Kabsch relative RMSD: {_fmt(s.get('kabsch_relative_rmsd'))} *(same-row correspondence diagnostic)*",
            f"- Pairwise-distance Wasserstein / median(A): {_fmt(s.get('pairwise_distance_wasserstein_relative'))}",
        ])
    m = report["metadata"]
    lines.extend([
        "",
        "## Non-scorer / container metadata",
        "",
        f"- obs columns only in A: `{m['obs_only_a']}`",
        f"- obs columns only in B: `{m['obs_only_b']}`",
        f"- layers only in A: `{m['layers_only_a']}`",
        f"- layers only in B: `{m['layers_only_b']}`",
        f"- obsm keys only in A: `{m['obsm_only_a']}`",
        f"- obsm keys only in B: `{m['obsm_only_b']}`",
        f"- uns keys only in A: `{m['uns_only_a']}`",
        f"- uns keys only in B: `{m['uns_only_b']}`",
        "",
        "## Interpretation notes",
        "",
        "- VEC Diff compares expression after casting to `float32`, matching the scorer contract.",
        "- Cell row order is not itself meaningful to the current VEC metrics. For T2/T3, however, each expression row remains paired with its own `spatial_3D` coordinate; the paired-cell multiset therefore matters.",
        "- Current VEC spatial metrics are invariant to translation and rotation. A rigid-frame-only classification is only emitted when expression is rowwise identical and a proper Kabsch alignment leaves negligible residual.",
        "- This tool is a comparison aid, not an official submission validator or scorer.",
    ])
    return "\n".join(lines) + "\n"
