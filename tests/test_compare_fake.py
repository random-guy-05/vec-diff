import numpy as np
import pandas as pd

import vec_diff.compare as cmp


class FakeFile:
    def close(self):
        pass


class FakeAdata:
    def __init__(self, expression, genes, coords=None):
        self.X = np.asarray(expression, dtype=np.float64)
        self.var_names = pd.Index(genes)
        self.n_obs, self.n_vars = self.X.shape
        self.obs = pd.DataFrame(index=range(self.n_obs))
        self.var = pd.DataFrame(index=genes)
        self.layers = {}
        self.obsm = (
            {}
            if coords is None
            else {"spatial_3D": np.asarray(coords, dtype=float)}
        )
        self.uns = {}
        self.file = FakeFile()


def test_compare_t2_detects_paired_row_reordering(tmp_path, monkeypatch):
    expression = np.array(
        [[1, 2], [3, 4], [5, 6]],
        dtype=float,
    )
    coords = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
        dtype=float,
    )
    a = FakeAdata(expression, ["g1", "g2"], coords)
    perm = [2, 0, 1]
    b = FakeAdata(
        expression[perm],
        ["g1", "g2"],
        coords[perm],
    )

    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    path_a.write_bytes(b"a")
    path_b.write_bytes(b"b")

    monkeypatch.setattr(
        cmp,
        "load_anndata_backed",
        lambda path: a if path.name == "a.h5ad" else b,
    )
    report = cmp.compare_files(path_a, path_b, "T2")
    assert report["classification"] == "ROW_ORDER_ONLY"
    assert report["spatial"]["paired_multiset_overlap"] == 1.0


def test_compare_t2_missing_coords_is_incompatible(tmp_path, monkeypatch):
    expression = np.array([[1, 2], [3, 4]], dtype=float)
    a = FakeAdata(
        expression,
        ["g1", "g2"],
        np.zeros((2, 3)),
    )
    b = FakeAdata(expression, ["g1", "g2"], None)

    path_a = tmp_path / "a.h5ad"
    path_b = tmp_path / "b.h5ad"
    path_a.write_bytes(b"a")
    path_b.write_bytes(b"b")

    monkeypatch.setattr(
        cmp,
        "load_anndata_backed",
        lambda path: a if path.name == "a.h5ad" else b,
    )
    report = cmp.compare_files(path_a, path_b, "T2")
    assert report["classification"] == "INCOMPATIBLE"
