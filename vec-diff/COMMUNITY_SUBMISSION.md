# Community Contribution submission text

## Title
VEC Diff — scorer-aware comparison of `.h5ad` prediction files

## Short description
VEC Diff compares two Virtual Embryo Challenge prediction files in the semantics that matter to the current scorer rather than merely comparing file bytes. It audits gene sets/order, compares expression after scorer-style float32 casting, detects exact cell-multiset equivalence despite row reordering, quantifies pseudobulk and single-cell differences, and for T2/T3 preserves expression-coordinate pairing while reporting tissue scale, geometry changes, and rigid rotation/translation-only differences. It also separates differences in scorer-relevant content from irrelevant AnnData metadata such as extra layers/uns fields and emits Markdown/JSON reports plus CI-friendly exit codes.

## Contribution / community value
During iterative model development it is easy to accumulate many `.h5ad` variants and accidentally treat metadata-only rewrites, reordered cells, or duplicated exports as distinct model candidates. Conversely, two files can have nearly identical pseudobulk expression while differing substantially in cell-state distribution or expression-to-position coupling. VEC Diff makes those distinctions explicit before competitors spend time evaluating or selecting submissions. It is open, installable, memory-bounded for large T1 files, and uses no hidden Challenge data.

## Suggested tags
VEC; tooling; AnnData; debugging; submission workflow; reproducibility
