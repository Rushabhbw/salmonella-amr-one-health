"""Tests for workflow/scripts/07_make_splits.py (run from the repository root: pytest)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "07_make_splits.py"
spec = importlib.util.spec_from_file_location("splits", SCRIPT)
sp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sp)


def fake_meta():
    rng = np.random.default_rng(0)
    rows = []
    for g in range(60):                                   # 60 groups of varying size
        size = int(rng.integers(1, 12))
        for j in range(size):
            clin = (g % 6 == 0 and j == 0) or (g % 15 == 0)   # some groups mix clinical and non-clinical
            rows.append({"isolate": f"G{g}_{j}", "in_primary_set": g != 59, "cluster_group": f"PDS{g}",
                         "host_group": "clinical" if clin else "non_clinical",
                         "collection_year": int(rng.integers(2005, 2021)) if not (g == 3 and j == 0) else pd.NA,
                         "bioproject": f"PRJ{g % 4}"})
    return pd.DataFrame(rows)


def test_group_folds_keep_groups_together_and_balance():
    m = fake_meta()
    out, shared = sp.make_splits(m)
    p = out[out["fold5"].notna()]
    assert p.groupby("cluster_group")["fold5"].nunique().max() == 1          # a group never spans two folds
    sizes = p["fold5"].value_counts()
    assert sizes.max() - sizes.min() <= p["cluster_group"].value_counts().max()   # balanced to within the largest group
    assert out.loc[~out["in_primary_set"], "fold5"].isna().all()               # non-primary isolates are not split


def test_transfer_disjoint_has_no_shared_groups_and_naive_does():
    m = fake_meta()
    out, shared = sp.make_splits(m)
    assert len(shared) > 0
    _, leaked_disjoint = sp.leakage(out, "transfer_disjoint", (["train"], ["test"]))
    _, leaked_naive = sp.leakage(out, "transfer_naive", (["train"], ["test"]))
    assert leaked_disjoint == 0 and leaked_naive > 0
    assert set(out.loc[out["transfer_disjoint"] == "test", "host_group"]) == {"clinical"}
    ex = out[out["transfer_disjoint"] == "excluded"]
    assert ex["cluster_group"].isin(shared).all() and (ex["host_group"] == "non_clinical").all()


def test_time_split_and_missing_year():
    m = fake_meta()
    out, _ = sp.make_splits(m, cutoff=2012)
    y = pd.to_numeric(out["collection_year"], errors="coerce")
    assert (y[out["time"] == "train"] <= 2012).all() and (y[out["time"] == "test"] > 2012).all()
    assert out.loc[y.isna(), "time"].isna().all()


def test_deterministic_and_end_to_end(tmp_path):
    m = fake_meta()
    a, _ = sp.make_splits(m); b, _ = sp.make_splits(m)
    assert a["fold5"].tolist() == b["fold5"].tolist()
    for d in sp.DRUGS:
        m[f"ast_{d}"] = np.where(np.arange(len(m)) % 5 == 0, "R", "S")
    m.to_csv(tmp_path / "m.tsv", sep="\t", index=False)
    code = sp.main(["--metadata", str(tmp_path / "m.tsv"), "--out", str(tmp_path / "s.tsv"), "--report", str(tmp_path / "r.md")])
    assert code == 0
    text = (tmp_path / "r.md").read_text(encoding="utf-8")
    assert "Leakage check" in text and "transfer_disjoint" in text


def test_random5_is_balanced_but_leaky():
    m = fake_meta()
    out, _ = sp.make_splits(m)
    p = out[out["random5"].notna()]
    sizes = p["random5"].value_counts()
    assert sizes.max() - sizes.min() <= 1                                    # balanced to within one isolate
    assert out.loc[~out["in_primary_set"], "random5"].isna().all()
    n_g, leak_g = sp.fold_leakage(out, "fold5")
    n_r, leak_r = sp.fold_leakage(out, "random5")
    assert leak_g == 0 and leak_r > 0.5 * n_r                                 # grouped folds are clean; random folds are not

