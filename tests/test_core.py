import numpy as np

from vec_diff.core import (
    classify,
    expression_differences,
    kabsch_rmsd,
    multiset_overlap,
    row_hashes_dense,
    spatial_stats,
)


def test_row_hash_order_invariant_overlap():
    x = np.array([[1, 2], [3, 4], [5, 6]], dtype=float)
    h1 = row_hashes_dense(x)
    h2 = row_hashes_dense(x[[2, 0, 1]])
    common, na, nb, frac = multiset_overlap(h1, h2)
    assert (common, na, nb) == (3, 3, 3)
    assert frac == 1.0


def test_pair_hash_preserves_coupling():
    x = np.array([[1, 2], [3, 4]], dtype=float)
    c = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    good = row_hashes_dense(x[[1, 0]], c[[1, 0]])
    bad = row_hashes_dense(x[[1, 0]], c)
    _, _, _, fg = multiset_overlap(row_hashes_dense(x, c), good)
    _, _, _, fb = multiset_overlap(row_hashes_dense(x, c), bad)
    assert fg == 1.0
    assert fb < 1.0


def test_kabsch_rotation_translation_is_near_zero():
    a = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]], dtype=float)
    r = np.array([[0,-1,0],[1,0,0],[0,0,1]], dtype=float)
    b = a @ r.T + np.array([8, -3, 2])
    out = kabsch_rmsd(a, b)
    assert out["relative_rmsd"] < 1e-6
    assert out["det_rotation"] > 0.999


def test_expression_float32_semantics():
    a = np.array([[1.0, 2.0]], dtype=np.float64)
    b = a.copy(); b[0, 0] += 1e-10
    d = expression_differences(a, b)
    assert d["exact_fraction"] == 1.0


def test_classify_t1_row_order_only():
    c, _ = classify(
        byte_identical=False,
        genes_compatible=True,
        task="T1",
        expression_rowwise_identical=False,
        expression_multiset_overlap=1.0,
        paired_multiset_overlap=None,
        kabsch_relative_rmsd=None,
        coord_rowwise_identical=None,
    )
    assert c == "ROW_ORDER_ONLY"


def test_classify_spatial_rigid_frame_only():
    c, _ = classify(
        byte_identical=False,
        genes_compatible=True,
        task="T2",
        expression_rowwise_identical=True,
        expression_multiset_overlap=1.0,
        paired_multiset_overlap=0.0,
        kabsch_relative_rmsd=1e-10,
        coord_rowwise_identical=False,
    )
    assert c == "RIGID_FRAME_ONLY"


def test_spatial_scale_detected():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(100, 3))
    b = a * 2
    s = spatial_stats(a, b, sample_size=100)
    assert abs(s["log_scale_b_over_a"] - np.log(2)) < 1e-5


def test_spatial_sampling_identical_cloud_is_zero():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(1000, 3))
    s = spatial_stats(a, a.copy(), sample_size=100, seed=7)
    assert s["pairwise_distance_wasserstein"] == 0.0
