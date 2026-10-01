import numpy as np
import pandas as pd

import vec_diff.compare as cmp


class FakeFile:
    def close(self):
        pass


class FakeAdata:
    def __init__(self, X, genes, coords=None):
        self.X = np.asarray(X, dtype=np.float64)
        self.var_names = pd.Index(genes)
        self.n_obs, self.n_vars = self.X.shape
        self.obs = pd.DataFrame(index=range(self.n_obs))
        self.var = pd.DataFrame(index=genes)
        self.layers = {}
        self.obsm = {} if coords is None else {"spatial_3D": np.asarray(coords, dtype=float)}
        self.uns = {}
        self.file = FakeFile()


def test_compare_t2_detects_paired_row_reordering(tmp_path, monkeypatch):
    x = np.array([[1,2],[3,4],[5,6]], dtype=float)
    c = np.array([[0,0,0],[1,0,0],[0,1,0]], dtype=float)
    a = FakeAdata(x, ["g1","g2"], c)
    perm = [2,0,1]
    b = FakeAdata(x[perm], ["g1","g2"], c[perm])
    pa = tmp_path / "a.h5ad"; pb = tmp_path / "b.h5ad"
    pa.write_bytes(b"a"); pb.write_bytes(b"b")
    monkeypatch.setattr(cmp, "load_anndata_backed", lambda p: a if p.name == "a.h5ad" else b)
    r = cmp.compare_files(pa, pb, "T2")
    assert r["classification"] == "ROW_ORDER_ONLY"
    assert r["spatial"]["paired_multiset_overlap"] == 1.0


def test_compare_t2_missing_coords_is_incompatible(tmp_path, monkeypatch):
    x = np.array([[1,2],[3,4]], dtype=float)
    a = FakeAdata(x, ["g1","g2"], np.zeros((2,3)))
    b = FakeAdata(x, ["g1","g2"], None)
    pa = tmp_path / "a.h5ad"; pb = tmp_path / "b.h5ad"
    pa.write_bytes(b"a"); pb.write_bytes(b"b")
    monkeypatch.setattr(cmp, "load_anndata_backed", lambda p: a if p.name == "a.h5ad" else b)
    r = cmp.compare_files(pa, pb, "T2")
    assert r["classification"] == "INCOMPATIBLE"
