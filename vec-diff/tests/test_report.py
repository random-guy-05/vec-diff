from vec_diff.report import markdown


def test_markdown_smoke():
    r = {
        "classification": "MATERIAL",
        "interpretation": "different",
        "task": "T1",
        "files": {"a": {"path": "a", "sha256": "1"}, "b": {"path": "b", "sha256": "2"}},
        "genes": {"same_set": True, "same_order": True, "genes_a": ["g"], "genes_b": ["g"], "missing_from_b": [], "extra_in_b": []},
        "expression": {"n_cells_a": 1, "n_cells_b": 1, "rowwise_identical": False, "expression_multiset_overlap": 0.0},
        "top_genes": [],
        "spatial": None,
        "metadata": {"obs_only_a": [], "obs_only_b": [], "layers_only_a": [], "layers_only_b": [], "obsm_only_a": [], "obsm_only_b": [], "uns_only_a": [], "uns_only_b": []},
    }
    text = markdown(r)
    assert "VEC Diff report" in text
    assert "MATERIAL" in text
