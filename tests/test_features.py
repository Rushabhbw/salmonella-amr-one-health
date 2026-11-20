"""Tests for workflow/scripts/04_build_features.py (run from the repository root: pytest)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "04_build_features.py"
spec = importlib.util.spec_from_file_location("build_features", SCRIPT)
bf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bf)


def meta(strings):
    rows = []
    for i, s in enumerate(strings):
        rows.append({"isolate": f"I{i}", "assembly": f"GCA_{i}.1", "in_primary_set": True,
                     "host_group": "clinical" if i % 2 else "non_clinical", "cluster_group": f"G{i}",
                     "amr_genes_core": s, "y_tetracycline": 1.0, "y_ampicillin": 0.0,
                     "y_ceftriaxone": np.nan, "y_ciprofloxacin": 0.0})
    return pd.DataFrame(rows)


def test_parse_calls():
    assert bf.parse_calls("blaTEM-1=COMPLETE,tet(A)=PARTIAL") == [("blaTEM-1", "COMPLETE"), ("tet(A)", "PARTIAL")]
    assert bf.parse_calls("") == [] and bf.parse_calls(None) == [] and bf.parse_calls(float("nan")) == []


def test_presence_uncertain_and_empty():
    m = meta(["blaTEM-1=COMPLETE,gyrA_S83F=POINT,tet(A)=PARTIAL", "", None,
              "sul1=COMPLETE,sul1=PARTIAL_END_OF_CONTIG", "aac(6')-Iaa=HMM"])
    present, uncertain, calls, unknown = bf.build_matrices(m)
    assert present.loc[0, "blaTEM-1"] == 1 and present.loc[0, "gyrA_S83F"] == 1
    assert "tet(A)" not in present.columns           # only weak calls, so not a feature
    assert uncertain.loc[0, "tet(A)"] == 1                       # weak call only
    assert present.loc[[1, 2]].sum().sum() == 0                  # empty and missing strings -> all zeros
    assert present.loc[3, "sul1"] == 1 and uncertain.loc[3, "sul1"] == 0   # strong call wins over weak call
    assert uncertain.loc[4, "aac(6')-Iaa"] == 1 and not unknown
    assert set(present.dtypes.astype(str)) == {"int8"}


def test_unknown_call_type_is_reported_and_weak():
    m = meta(["geneX=NEWCALL"])
    present, uncertain, _, unknown = bf.build_matrices(m)
    assert unknown == {"NEWCALL": 1} and uncertain.loc[0, "geneX"] == 1 and present.shape[1] == 0   # weak-only elements are not features


def test_dictionary_counts_and_kind():
    m = meta(["blaTEM-1=COMPLETE,gyrA_S83F=POINT", "blaTEM-1=COMPLETE", "gyrA_S83F=POINT"])
    present, uncertain, calls, _ = bf.build_matrices(m)
    d = bf.feature_dictionary(m, present, uncertain, calls).set_index("element")
    assert d.loc["blaTEM-1", "n_isolates_all"] == 2 and d.loc["blaTEM-1", "n_groups_all"] == 2
    assert d.loc["gyrA_S83F", "kind"] == "point_mutation" and d.loc["blaTEM-1", "kind"] == "gene"
    assert d.index[0] in ("blaTEM-1", "gyrA_S83F")


def test_end_to_end(tmp_path):
    m = meta(["blaTEM-1=COMPLETE", "", "sul1=COMPLETE,tet(A)=PARTIAL"])
    mp = tmp_path / "m.tsv"; m.to_csv(mp, sep="\t", index=False)
    code = bf.main(["--metadata", str(mp), "--out-dir", str(tmp_path), "--report", str(tmp_path / "r.md")])
    assert code == 0
    f = pd.read_csv(tmp_path / "features_presence.tsv", sep="\t")
    assert len(f) == 3 and {"isolate", "blaTEM-1", "sul1"} <= set(f.columns)
    u = pd.read_csv(tmp_path / "features_uncertain.tsv", sep="\t")
    assert "tet(A)" in u.columns and u.loc[2, "tet(A)"] == 1
    assert "Isolates with no element at all: 1" in (tmp_path / "r.md").read_text(encoding="utf-8")
