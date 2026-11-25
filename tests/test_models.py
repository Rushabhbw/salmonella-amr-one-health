"""Tests for workflow/scripts/08_train_evaluate_models.py (run from the repository root: pytest)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, ROOT / "workflow" / "scripts" / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ml = load("train_models", "08_train_evaluate_models.py")
sp = load("splits_for_models", "07_make_splits.py")
bl = ml.bl

ELEMENTS = ["tet(A)", "blaTEM-1", "blaCMY-2", "sul2", "aph(6)-Id", "qnrB19", "tet(B)", "rareX"]


def make_world(n_groups=70, seed=0):
    rng = np.random.default_rng(seed)
    rows, feats = [], []
    for g in range(n_groups):
        size = int(rng.integers(2, 9))
        # a lineage carries a fixed backbone of genes; individuals vary a little
        backbone = rng.random(len(ELEMENTS)) < 0.35
        for j in range(size):
            f = backbone.copy()
            flip = rng.random(len(ELEMENTS)) < 0.05
            f = np.where(flip, ~f, f).astype(int)
            clin = (g % 7 == 0 and j == 0) or g % 17 == 0
            noise = rng.random(3) < 0.04
            y_tet = int(f[0]) ^ int(noise[0])
            y_amp = int(f[1] or f[2]) ^ int(noise[1])
            y_cef = int(f[2]) ^ int(noise[2])
            rows.append({"isolate": f"G{g}_{j}", "assembly": f"GCA_{g}_{j}.1", "in_primary_set": g != n_groups - 1,
                         "host_group": "clinical" if clin else "non_clinical", "cluster_group": f"PDS{g}",
                         "collection_year": int(rng.integers(2005, 2021)), "bioproject": f"PRJ{g % 5}",
                         "y_tetracycline": float(y_tet), "y_ampicillin": float(y_amp), "y_ceftriaxone": float(y_cef), "y_ciprofloxacin": 0.0,
                         "ast_tetracycline": "R" if y_tet else "S", "ast_ampicillin": "R" if y_amp else "S",
                         "ast_ceftriaxone": "R" if y_cef else "S", "ast_ciprofloxacin": "S"})
            feats.append(f)
    meta = pd.DataFrame(rows)
    F = pd.DataFrame(np.array(feats), columns=ELEMENTS)
    mapping = pd.DataFrame({"element": ELEMENTS,
                            "tetracycline": [1, 0, 0, 0, 0, 0, 1, 0], "ampicillin": [0, 1, 1, 0, 0, 0, 0, 0],
                            "ceftriaxone": [0, 0, 1, 0, 0, 0, 0, 0], "ciprofloxacin": [0, 0, 0, 0, 0, 1, 0, 0]})
    for d in ["tetracycline", "ampicillin", "ceftriaxone", "ciprofloxacin"]:
        mapping[f"{d}_with_multidrug"] = mapping[d]
    return meta, F, mapping


def prepare():
    meta, F, mapping = make_world()
    splits, _ = sp.make_splits(meta)
    feats_full = pd.concat([meta[["isolate", "assembly", "in_primary_set", "host_group", "cluster_group"]], F], axis=1)
    preds, _ = bl.build(meta, feats_full, mapping)
    prim = meta["in_primary_set"].to_numpy()
    U = pd.DataFrame(np.zeros((len(meta), 1), dtype=int), columns=["tet(A)"])
    return (meta[prim].reset_index(drop=True), F[prim].reset_index(drop=True), U[prim].reset_index(drop=True), mapping,
            splits[prim].reset_index(drop=True), preds[prim].reset_index(drop=True), feats_full, meta, splits, preds)


def test_grouped_folds_never_split_a_cluster_and_baseline_is_aligned():
    meta, F, U, mapping, splits, base, *_ = prepare()
    p, m, imp = ml.run_all(meta, F, U, mapping, splits, base, fast=True, n_boot=60, drugs=["tetracycline"],
                           variants={"lr": ("lr", "presence"), "rf": ("rf", "presence")})
    g = p[(p.split == "fold5") & (p.model == "rf")]
    assert g.groupby("cluster_group")["fold"].nunique().max() == 1          # a group is tested in exactly one fold
    r = p[(p.split == "random5") & (p.model == "rf")]
    assert r.groupby("cluster_group")["fold"].nunique().max() > 1            # random folds split lineages
    # baseline predictions on the test isolates are the baseline's own predictions
    chk = g.merge(base[["isolate", "pred_tetracycline"]], on="isolate")
    assert (chk["baseline_pred"] == chk["pred_tetracycline"]).all()
    # the baseline row and the model rows describe the same test isolates
    row = m[(m.split == "fold5")]
    assert row["n_R"].nunique() == 1 and row["n_S"].nunique() == 1
    assert set(row["model"]) == {"rules (baseline)", "lr", "rf"}


def test_models_learn_an_obvious_signal_and_shortcut_variant_drops_catalog_features():
    meta, F, U, mapping, splits, base, *_ = prepare()
    p, m, imp = ml.run_all(meta, F, U, mapping, splits, base, fast=True, n_boot=60, drugs=["tetracycline"],
                           variants={"lr": ("lr", "presence"), "rf_no_catalog": ("rf", "no_catalog")})
    lr = m[(m.model == "lr") & (m.split == "fold5")].iloc[0]
    assert lr["mcc"] > 0.6
    linked = ["tet(A)", "tet(B)"]
    X = ml.feature_matrix("no_catalog", F, U, linked)
    assert "tet(A)" not in X.columns and "tet(B)" not in X.columns and X.shape[1] == F.shape[1] - 2
    assert ml.feature_matrix("presence+partial", F, U, linked).shape[1] == F.shape[1] + U.shape[1]
    row = imp[(imp.model == "lr") & (imp.element == "tet(A)")].iloc[0]
    assert bool(row["catalog_linked"]) and 0 < row["share"] < 1
    assert abs(imp[(imp.model == "lr")]["share"].sum() - 1) < 1e-5


def test_paired_difference_of_identical_predictions_is_zero():
    counts = np.array([[5, 1, 20, 2]] * 15 + [[1, 3, 8, 1]] * 5, dtype=float)
    lo, hi = ml.paired_delta_mcc(counts, counts.copy(), 200)
    assert lo == 0 and hi == 0
    better = counts.copy(); better[:, [0, 1]] = better[:, [1, 0]]            # swap TP and FN: clearly worse sensitivity
    lo2, hi2 = ml.paired_delta_mcc(better, counts, 200)
    assert hi2 < 0


def test_end_to_end_fast_run_writes_all_outputs(tmp_path):
    meta, F, U, mapping, splits, base, feats_full, meta_all, splits_all, preds_all = prepare()
    meta_all.to_csv(tmp_path / "m.tsv", sep="\t", index=False)
    feats_full.to_csv(tmp_path / "f.tsv", sep="\t", index=False)
    splits_all.to_csv(tmp_path / "s.tsv", sep="\t", index=False)
    preds_all.to_csv(tmp_path / "b.tsv", sep="\t", index=False)
    mapping.to_csv(tmp_path / "map.tsv", sep="\t", index=False)
    code = ml.main(["--metadata", str(tmp_path / "m.tsv"), "--features", str(tmp_path / "f.tsv"), "--uncertain", str(tmp_path / "none.tsv"),
                    "--splits", str(tmp_path / "s.tsv"), "--baseline", str(tmp_path / "b.tsv"), "--mapping", str(tmp_path / "map.tsv"),
                    "--out-dir", str(tmp_path), "--report", str(tmp_path / "ml.md"), "--fast"])
    assert code == 0
    text = (tmp_path / "ml.md").read_text(encoding="utf-8")
    assert "FAST TEST RUN" in text and "## tetracycline" in text and "Shortcut test" in text
    assert (tmp_path / "ml_metrics_fast.tsv").exists() and (tmp_path / "ml_predictions_fast.tsv.gz").exists()
    assert ml.main(["--mapping", str(tmp_path / "missing.tsv")]) == 1


def test_rules_feature_and_hybrid_models():
    meta, F, U, mapping, splits, base, *_ = prepare()
    rule = np.where(F["tet(A)"] + F["tet(B)"] > 0, 1, 0)
    X = ml.feature_matrix("presence+rules", F, U, [], rule)
    assert list(X.columns)[-1] == "rule_evidence" and X.shape[1] == F.shape[1] + 1 and (X["rule_evidence"].to_numpy() == rule).all()
    assert "lr_rules" in ml.VARIANTS and "hgb_rules" in ml.EXPLORATORY
    p, m, imp = ml.run_all(meta, F, U, mapping, splits, base, fast=True, n_boot=60, drugs=["tetracycline"],
                           variants={"lr_rules": ("lr", "presence+rules")})
    assert m[(m.model == "lr_rules") & (m.split == "fold5")].iloc[0]["mcc"] > 0.6

