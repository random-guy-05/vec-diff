# VEC Diff

**Scorer-aware diff for Virtual Embryo Challenge `.h5ad` predictions.**

`sha256sum` tells you two submission files differ. VEC Diff tells you **whether they differ in content the VEC scorer actually uses, and how**.

> Independent community tool. Not an official Virtual Embryo Challenge validator or scorer and not affiliated with the organizers.

## Why

Late in model development, directories often contain variants such as:

```text
best.h5ad
best_fixed.h5ad
best_v2.h5ad
best_v2_final.h5ad
```

A byte diff is not enough. Two files may differ only in `uns`, compression, cell row order, or coordinate frame, while another pair may have nearly identical pseudobulk expression but a different single-cell distribution or expression-to-position coupling.

VEC Diff follows the current VEC submission semantics:

- Expression is compared after `float32` casting, matching the scorer contract.
- Gene sets/order are audited explicitly.
- T1 cell row order is treated as irrelevant to distributional scoring.
- T2/T3 preserve the **pairing between each expression row and its own first three `spatial_3D` coordinates**.
- Exact row multisets are compared independent of row order.
- Spatial reports include tissue scale and a rotation/translation-invariant shape sketch.
- A proper Kabsch alignment detects coordinate changes that are only a global rigid frame transform, which current VEC spatial metrics are designed to ignore.
- Changes in `obs`, `layers`, extra `obsm`, `uns`, etc. are listed separately from scorer-relevant content.

## Install

```bash
git clone https://github.com/random-guy-05/vec-diff.git
cd vec-diff
python -m pip install -e .
```

## Quick start

```bash
vec-diff model_A.h5ad model_B.h5ad --task T2
```

Write reproducible reports:

```bash
vec-diff A.h5ad B.h5ad \
  --task T2 \
  --json diff.json \
  --markdown diff.md
```

## Classifications

VEC Diff emits one high-level classification:

- **`BYTE_IDENTICAL`** — the files themselves have the same SHA256.
- **`SCORER_CONTENT_IDENTICAL`** — scorer-relevant values are identical after current scorer-style casting/alignment; container metadata may differ.
- **`ROW_ORDER_ONLY`** — the same exact predicted cell multiset exists in a different row order. For T2/T3 this is checked on the paired `(expression, spatial_3D[:3])` cell, so expression-coordinate coupling is preserved.
- **`RIGID_FRAME_ONLY`** — spatial task only: expression is rowwise identical and coordinates differ only by a proper rotation/translation within a strict numerical tolerance. Current VEC spatial metrics are rotation/translation invariant.
- **`MATERIAL`** — expression and/or spatial coupling materially differs.
- **`INCOMPATIBLE`** — the gene sets differ, so a scorer-aware equivalence claim cannot be made.

These labels are intentionally conservative. VEC Diff does **not** claim that two merely similar distributions will receive identical scores.

## Example

```text
ROW_ORDER_ONLY
The exact (expression, spatial_3D[:3]) paired-cell multiset is identical;
only row order and/or non-scorer metadata differ.

Expression
  expression_multiset_overlap       1
  pseudobulk_pearson                 1
  pseudobulk_mean_abs_delta          0

Spatial
  paired_multiset_overlap            1
  log_scale_b_over_a                 0
```

Or:

```text
MATERIAL
Expression and/or expression-to-position coupling differs in scorer-relevant content.

Expression
  expression_multiset_overlap       0.744
  pseudobulk_pearson                 0.9986

Spatial
  paired_multiset_overlap            0.120
  log_scale_b_over_a                 0.087
```

The second example makes an important point: high pseudobulk correlation does not mean two generative submissions are the same.

## Machine-readable / CI use

```bash
# Return exit code 3 for recognized scorer-equivalent variants.
vec-diff A.h5ad B.h5ad --task T1 --fail-if-equivalent

# Return exit code 4 when scorer-relevant content materially differs.
vec-diff A.h5ad B.h5ad --task T2 --fail-if-material
```

This makes VEC Diff usable in checkpoint/export pipelines.

## What is compared

### Gene contract

- counts, exact names, duplicates, order
- same-set/different-order detection
- logical alignment for the comparison when the same unique gene set is present

### Expression

- cell count
- `float32`-canonical exact row hashes
- exact row-multiset overlap
- same-row MAE / RMSE / max difference when shapes match
- pseudobulk Pearson
- pseudobulk mean/max absolute change
- genes with the largest pseudobulk delta

### Spatial (T2/T3)

Only the first three columns of `obsm["spatial_3D"]` are scorer-relevant under the current file contract.

VEC Diff reports:

- exact coordinate equality
- exact paired `(X row, coordinate row)` multiset overlap
- centroid shift (reported but labeled frame-dependent)
- RMS tissue radius and log scale ratio
- proper-Kabsch residual for same-row rigid-frame changes
- a deterministic pairwise-distance Wasserstein sketch, which is invariant to translation/rotation and helps reveal shape/scale changes

### Container metadata

Differences in `obs`, `var` columns, `layers`, extra `obsm`, and `uns` keys are shown separately so a user can understand why byte hashes differ even when scorer-relevant content does not.

## Memory behavior

Expression comparisons use backed AnnData and bounded row chunks rather than densifying an entire T1 transcriptome in memory. Sparse chunks are densified only one bounded chunk at a time because the official scorer itself ultimately evaluates numeric expression values after `float32` casting.

## Limitations

- This is not an official format validator. Use the Challenge's official checks / scorer for validity.
- Exact multiset hashing deliberately tests exact scorer-cast cell values. Similar-but-not-identical cells are summarized numerically rather than declared equivalent.
- `RIGID_FRAME_ONLY` currently requires same row correspondence for the Kabsch diagnostic. Detecting an arbitrary cell permutation *plus* an arbitrary rigid transform would require a registration/matching problem and is intentionally not guessed.
- Current VEC rules/evaluation are authoritative and may change. This release targets the public contract documented on **2026-09-30**.

## Development

```bash
python -m pip install -e '.[dev]'
pytest
ruff check src tests
```

## Sources

- Challenge submissions/evaluation: https://virtualembryo.ai/challenge/evaluation
- Challenge rules: https://virtualembryo.ai/challenge/rules
- Community Contribution Award: https://virtualembryo.ai/challenge/community

## License

MIT.
