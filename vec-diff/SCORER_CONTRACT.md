# VEC scorer contract snapshot used by VEC Diff 1.0.0

Snapshot date: **2026-09-30**.

Authoritative source: https://virtualembryo.ai/challenge/evaluation

VEC Diff's scorer-aware comparison is based on these current public properties:

- Submission expression lives in AnnData `.X` and is evaluated as non-negative finite `float32` values.
- Gene names/order are part of the submission contract.
- T1 submits expression only.
- T2/T3 require `obsm["spatial_3D"]`; only the first three coordinate columns are read by the current scorer contract.
- Submitted cell-type labels are ignored by the official scorer.
- Predicted cell count does not have to equal truth cell count.
- Current metrics do not assume predicted cell `i` corresponds to truth cell `i`; row order alone therefore carries no biological meaning.
- For spatial tasks, each expression row is nevertheless coupled to its own coordinate row. VEC Diff hashes that pair together when testing row-order-only equivalence.
- Current spatial terms are designed to be invariant to global translation and rotation. VEC Diff uses a proper Kabsch alignment only as a conservative diagnostic for same-row rigid-frame changes.

The official Challenge documentation remains authoritative if the contract changes after this snapshot.
