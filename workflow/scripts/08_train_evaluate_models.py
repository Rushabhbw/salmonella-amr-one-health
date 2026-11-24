#!/usr/bin/env python3
"""
08_train_evaluate_models.py

Train machine-learning models to predict resistance from the resistance-element table, evaluate them on the fixed
splits made by 07_make_splits.py, and compare every model with the rule-based baseline on exactly the same test isolates.

Antibiotics: tetracycline, ampicillin, ceftriaxone (ciprofloxacin is descriptive only, see the charter).
Dataset    : the primary dataset (9,242 isolates); only isolates with an S or R label for the antibiotic are used.

Model variants
  lr             logistic regression (L2, C = 1)                  presence/absence of the 149 elements
  rf             random forest (300 trees)                         same features
  hgb            histogram gradient boosting (200 iterations)      same features
  hgb_partial    histogram gradient boosting                       elements + partial-call indicators
  rf_no_catalog  random forest WITHOUT the elements NCBI links to the antibiotic: a shortcut test. It can only use
                 co-occurring genes, so it shows how much a model can learn from linkage alone.
Hyperparameters are fixed in advance (no tuning on any test data). Threshold 0.5, no class re-weighting.

Splits (from data/processed/splits.tsv)
  fold5              grouped 5-fold CV (no SNP-cluster group in train and test)
  transfer_disjoint  train on non-clinical isolates sharing no group with a clinical isolate, test on clinical isolates
  time               train up to 2016, test later
  random5            LEAKY random 5-fold, shown only to measure how much unfair evaluation inflates results

Comparison with the baseline: sensitivity, specificity, MCC with cluster-bootstrap intervals (whole SNP-cluster groups
are resampled), and the PAIRED difference in MCC (model minus rules) with its own interval, resampling the same groups
for both. ROC-AUC and average precision are pooled over the test isolates (no interval).

Reads : data/processed/{metadata_curated,features_presence,features_uncertain,splits,baseline_predictions}.tsv,
        config/element_classes.tsv
Writes: data/processed/ml_metrics.tsv, data/processed/ml_predictions.tsv.gz, data/processed/ml_feature_importance.tsv,
        reports/ml_report.md

Run from the repository root:
    python workflow/scripts/08_train_evaluate_models.py            # full run (several minutes)
    python workflow/scripts/08_train_evaluate_models.py --fast     # small models and few resamples, for testing only
"""
import argparse
import importlib.util
import os
import warnings
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = Path(__file__).resolve().parent


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


bl = _load("baseline_rules", "06_baseline_rules.py")        # metrics and cluster bootstrap are shared with the baseline

DRUGS = ["tetracycline", "ampicillin", "ceftriaxone"]
SPLITS = ["fold5", "transfer_disjoint", "time", "random5"]
VARIANTS = OrderedDict([
    ("lr", ("lr", "presence")),
    ("rf", ("rf", "presence")),
    ("hgb", ("hgb", "presence")),
    ("hgb_partial", ("hgb", "presence+partial")),
    ("rf_no_catalog", ("rf", "no_catalog")),
])
SEED = 42
N_BOOT = 1000
N_JOBS = int(os.environ.get("N_JOBS", "-1"))
ID_COLS = ("isolate", "assembly", "in_primary_set", "host_group", "cluster_group")


# ----------------------------------------------------------------------------- models and features
def make_model(algo, fast=False, seed=SEED):
    if algo == "lr":
        return LogisticRegression(C=1.0, max_iter=2000)
    if algo == "rf":
        return RandomForestClassifier(n_estimators=30 if fast else 300, max_features="sqrt", n_jobs=N_JOBS, random_state=seed)
    if algo == "hgb":
        return HistGradientBoostingClassifier(max_iter=25 if fast else 200, learning_rate=0.1, random_state=seed)
    raise ValueError(algo)


def feature_matrix(feature_set, F, U, linked):
    """F, U: DataFrames of 0/1 elements (strong, partial). linked: elements to drop for the shortcut test."""
    if feature_set == "presence":
        return F
    if feature_set == "presence+partial":
        return pd.concat([F, U.add_suffix("__partial")], axis=1)
    if feature_set == "no_catalog":
        return F.drop(columns=[c for c in F.columns if c in set(linked)])
    raise ValueError(feature_set)


def split_pairs(split_col, name):
    """Yield (fold label, train mask, test mask) for one split column of splits.tsv."""
    if name in ("fold5", "random5"):
        for k in sorted(split_col.dropna().unique()):
            yield int(k), (split_col.notna() & (split_col != k)).to_numpy(), (split_col == k).to_numpy()
    else:
        yield 0, (split_col == "train").to_numpy(), (split_col == "test").to_numpy()


# ----------------------------------------------------------------------------- evaluation
def paired_delta_mcc(counts_a, counts_b, n_boot, seed=SEED, chunk=250):
    """95% interval for MCC(a) - MCC(b), resampling the same groups for both."""
    g = len(counts_a)
    if g < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    deltas = []
    for start in range(0, n_boot, chunk):
        n = min(chunk, n_boot - start)
        idx = rng.integers(0, g, size=(n, g))
        a = bl.metrics_from_counts(counts_a[idx].sum(axis=1))["mcc"]
        b = bl.metrics_from_counts(counts_b[idx].sum(axis=1))["mcc"]
        deltas.append(a - b)
    d = np.concatenate(deltas)
    d = d[~np.isnan(d)]
    return (float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))) if len(d) > 10 else (np.nan, np.nan)


def counts_by_group(frame, pred_col):
    tmp = pd.DataFrame({"g": frame["cluster_group"], "y": frame["y"].astype(int), "p": frame[pred_col].astype(int)})
    c = pd.DataFrame({"TP": ((tmp.y == 1) & (tmp.p == 1)), "FN": ((tmp.y == 1) & (tmp.p == 0)),
                      "TN": ((tmp.y == 0) & (tmp.p == 0)), "FP": ((tmp.y == 0) & (tmp.p == 1))}).astype(int)
    c["g"] = tmp["g"].values
    return c.groupby("g")[["TP", "FN", "TN", "FP"]].sum().sort_index().to_numpy()


def safe_auc(y, p):
    if len(np.unique(y)) < 2:
        return np.nan, np.nan
    return float(roc_auc_score(y, p)), float(average_precision_score(y, p))


def metric_record(drug, split, model, frame, pred_col, base_counts, n_boot, with_delta):
    counts = counts_by_group(frame, pred_col)
    tot = counts.sum(axis=0)
    point = bl.metrics_from_counts(tot)
    ci = bl.cluster_bootstrap(counts, n_boot=n_boot)
    rec = {"drug": drug, "split": split, "model": model, "n_R": int(tot[0] + tot[1]), "n_S": int(tot[2] + tot[3]),
           "groups": len(counts), "TP": int(tot[0]), "FN": int(tot[1]), "TN": int(tot[2]), "FP": int(tot[3]),
           "sensitivity": float(point["sensitivity"]), "specificity": float(point["specificity"]), "mcc": float(point["mcc"]),
           "sens_lo": ci["sensitivity"][0], "sens_hi": ci["sensitivity"][1], "spec_lo": ci["specificity"][0], "spec_hi": ci["specificity"][1],
           "mcc_lo": ci["mcc"][0], "mcc_hi": ci["mcc"][1], "d_mcc": np.nan, "d_mcc_lo": np.nan, "d_mcc_hi": np.nan, "verdict": "",
           "auc": np.nan, "ap": np.nan}
    if with_delta:
        base_mcc = float(bl.metrics_from_counts(base_counts.sum(axis=0))["mcc"])
        lo, hi = paired_delta_mcc(counts, base_counts, n_boot)
        rec.update({"d_mcc": rec["mcc"] - base_mcc, "d_mcc_lo": lo, "d_mcc_hi": hi})
        rec["verdict"] = "no interval" if np.isnan(lo) else ("better than rules" if lo > 0 else ("worse than rules" if hi < 0 else "no clear difference"))
        if "prob" in frame.columns:
            rec["auc"], rec["ap"] = safe_auc(frame["y"].astype(int), frame["prob"])
    return rec


# ----------------------------------------------------------------------------- main computation
def run_all(meta, F, U, mapping, splits, base, fast=False, n_boot=N_BOOT, drugs=DRUGS, split_names=SPLITS, variants=VARIANTS):
    meta = meta.reset_index(drop=True)
    cg = meta["cluster_group"].to_numpy()
    iso = meta["isolate"].to_numpy()
    preds, records, importance = [], [], []
    for drug in drugs:
        y_all = pd.to_numeric(meta[f"y_{drug}"], errors="coerce")
        linked = set(bl.evidence_elements(mapping, drug, F.columns)) | set(bl.evidence_elements(mapping, f"{drug}_with_multidrug", F.columns))
        base_pred = base[f"pred_{drug}"].to_numpy()
        for vname, (algo, fs) in variants.items():
            X = feature_matrix(fs, F, U, linked).to_numpy(dtype=np.float32)
            for split in split_names:
                col = splits[split]
                chunks = []
                for fold, tr, te in split_pairs(col, split):
                    tr = tr & y_all.notna().to_numpy()
                    te = te & y_all.notna().to_numpy()
                    if tr.sum() < 20 or te.sum() == 0 or len(np.unique(y_all[tr])) < 2:
                        continue
                    model = make_model(algo, fast)
                    model.fit(X[tr], y_all.to_numpy()[tr].astype(int))
                    p = model.predict_proba(X[te])[:, 1]
                    chunks.append(pd.DataFrame({"isolate": iso[te], "cluster_group": cg[te], "drug": drug, "split": split, "model": vname,
                                                "fold": fold, "y": y_all.to_numpy()[te].astype(int), "prob": p,
                                                "pred": (p >= 0.5).astype(int), "baseline_pred": base_pred[te].astype(int)}))
                if not chunks:
                    continue
                frame = pd.concat(chunks, ignore_index=True)
                preds.append(frame)
                bcounts = counts_by_group(frame, "baseline_pred")
                if vname == "lr":       # the baseline is evaluated once per drug and split, on the same test isolates
                    records.append(metric_record(drug, split, "rules (baseline)", frame, "baseline_pred", bcounts, n_boot, with_delta=False))
                records.append(metric_record(drug, split, vname, frame, "pred", bcounts, n_boot, with_delta=True))
        # feature importance on all labelled primary isolates (interpretation only, not evaluation)
        lab = y_all.notna().to_numpy()
        for algo in ("lr", "rf"):
            Xf = F.to_numpy(dtype=np.float32)
            model = make_model(algo, fast)
            model.fit(Xf[lab], y_all.to_numpy()[lab].astype(int))
            w = np.abs(model.coef_[0]) if algo == "lr" else model.feature_importances_
            share = w / w.sum() if w.sum() > 0 else w
            for e, wi, si in zip(F.columns, w, share):
                importance.append({"drug": drug, "model": algo, "element": e, "weight": float(wi), "share": float(si),
                                   "catalog_linked": e in linked, "n_isolates": int(F[e].sum())})
    return (pd.concat(preds, ignore_index=True) if preds else pd.DataFrame()), pd.DataFrame(records), pd.DataFrame(importance)


# ----------------------------------------------------------------------------- report
def f3(v):
    return "NA" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.3f}"


def with_ci(v, lo, hi):
    return f"{f3(v)} ({f3(lo)} to {f3(hi)})" if not any(isinstance(x, float) and np.isnan(x) for x in (lo, hi)) else f3(v)


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def write_report(path, metrics, importance, fast):
    L = ["# Machine-learning report\n", f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}  \nScript: `workflow/scripts/08_train_evaluate_models.py`"
         + ("  \n**FAST TEST RUN: small models, few resamples. Do not interpret these numbers.**" if fast else "") + "\n",
         "## Method\n",
         "- Models, features and hyperparameters are fixed in advance (see the script header). Nothing is tuned on test data. Threshold 0.5.\n"
         "- Same isolates, same splits and the same labels as the rule-based baseline; every model is compared with the rules on exactly the same test isolates.\n"
         "- MCC and sensitivity/specificity have cluster-bootstrap intervals (whole SNP-cluster groups resampled). 'Difference in MCC' is model minus rules, "
         "with a paired interval; 'better' or 'worse' means the interval excludes 0. AUC is pooled over test isolates and has no interval. "
         f"Intervals here use {N_BOOT} resamples (2,000 in the baseline report), so the rules row can differ from the baseline report in the third decimal.\n"
         "- `random5` is a deliberately leaky split (about 92% of test isolates have a close relative in training). It shows how much an unfair evaluation inflates results and is not used to judge any model.\n"
         "- `rf_no_catalog` cannot see the elements NCBI links to the antibiotic. It shows what linkage between genes alone can predict.\n"]
    for drug in DRUGS:
        m = metrics[metrics["drug"] == drug]
        if m.empty:
            continue
        L.append(f"## {drug}\n")
        rows = []
        for split in SPLITS:
            for _, r in m[m["split"] == split].iterrows():
                rows.append([split, r["model"], f"{r['n_R']} / {r['n_S']}", r["groups"], f3(r["sensitivity"]), f3(r["specificity"]),
                             with_ci(r["mcc"], r["mcc_lo"], r["mcc_hi"]),
                             "" if r["model"].startswith("rules") else with_ci(r["d_mcc"], r["d_mcc_lo"], r["d_mcc_hi"]),
                             r["verdict"], "" if np.isnan(r["auc"]) else f3(r["auc"])])
        L.append(md_table(pd.DataFrame(rows, columns=["Split", "Model", "Test R / S", "Groups", "Sensitivity", "Specificity", "MCC (95% CI)",
                                                      "Difference in MCC vs rules (95% CI)", "Verdict", "AUC"])) + "\n")
    # leak and shortcut summaries
    L.append("## How much does an unfair split inflate the result? (MCC, random5 minus fold5)\n")
    rows = []
    for drug in DRUGS:
        for model in list(VARIANTS):
            a = metrics[(metrics.drug == drug) & (metrics.model == model) & (metrics.split == "random5")]
            b = metrics[(metrics.drug == drug) & (metrics.model == model) & (metrics.split == "fold5")]
            if len(a) and len(b):
                rows.append([drug, model, f3(float(b.iloc[0]["mcc"])), f3(float(a.iloc[0]["mcc"])), f3(float(a.iloc[0]["mcc"]) - float(b.iloc[0]["mcc"]))])
    if rows:
        L.append(md_table(pd.DataFrame(rows, columns=["Antibiotic", "Model", "MCC grouped (fair)", "MCC random (leaky)", "Inflation"])) + "\n")
    L.append("## Shortcut test: what can co-occurring genes alone predict? (rf_no_catalog)\n")
    rows = []
    for drug in DRUGS:
        for split in SPLITS:
            r = metrics[(metrics.drug == drug) & (metrics.model == "rf_no_catalog") & (metrics.split == split)]
            if len(r):
                r = r.iloc[0]
                rows.append([drug, split, f"{r['n_R']} / {r['n_S']}", f3(r["sensitivity"]), f3(r["specificity"]), with_ci(r["mcc"], r["mcc_lo"], r["mcc_hi"])])
    if rows:
        L.append(md_table(pd.DataFrame(rows, columns=["Antibiotic", "Split", "Test R / S", "Sensitivity", "Specificity", "MCC (95% CI)"])) + "\n")
        L.append("A high MCC here means the antibiotic's resistance can be guessed from other genes (linkage). Compare `random5` with `fold5`, `transfer_disjoint` and `time`: "
                 "a gap means the shortcut works inside known lineages and does not transfer to new ones.\n")
    # importance
    L.append("## What do the models rely on? (trained on all labelled primary isolates; for interpretation only)\n")
    rows = []
    for drug in DRUGS:
        for model in ("lr", "rf"):
            d = importance[(importance.drug == drug) & (importance.model == model)].sort_values("weight", ascending=False)
            if d.empty:
                continue
            top = d.head(6)
            linked_share = float(d.loc[d.catalog_linked, "share"].sum())
            rows.append([drug, model, f"{100 * linked_share:.0f}%", ", ".join(f"{e}{'' if l else '*'}" for e, l in zip(top["element"], top["catalog_linked"]))])
    if rows:
        L.append(md_table(pd.DataFrame(rows, columns=["Antibiotic", "Model", "Weight on catalog-linked elements", "Top 6 elements (* = not linked by the catalog)"])) + "\n")
        L.append("Logistic regression weight is the absolute coefficient; random forest weight is impurity importance. Both are rough guides, and correlated genes share weight.\n")
    L.append("## How to read this\n")
    L.append("- The rule-based baseline is strong, so a model can only help where the rules miss resistance. 'No clear difference' is an honest result.\n"
             "- Only `fold5`, `transfer_disjoint` and `time` count. `transfer_disjoint` tests only 270 human isolates, mostly from one 2014 study, so its intervals are wide.\n"
             "- No hyperparameter search was done. If a model looks underfitted, tuning would need nested grouped cross-validation.\n"
             "- Labels are submitter-supplied; models can learn label noise as easily as biology.\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(L), encoding="utf-8")


def load_inputs(a):
    meta = pd.read_csv(a.metadata, sep="\t")
    splits = pd.read_csv(a.splits, sep="\t")
    feats = pd.read_csv(a.features, sep="\t")
    base = pd.read_csv(a.baseline, sep="\t")
    mapping = pd.read_csv(a.mapping, sep="\t")
    unc = pd.read_csv(a.uncertain, sep="\t") if Path(a.uncertain).exists() else None
    primary = meta["in_primary_set"].astype(str).str.lower().eq("true").to_numpy()
    order = meta["isolate"]
    f = feats.set_index("isolate").loc[order]
    F = f[[c for c in f.columns if c not in ID_COLS]].astype(int).reset_index(drop=True)
    U = pd.DataFrame(index=F.index)
    if unc is not None:
        u = unc.set_index("isolate").reindex(order)
        U = u[[c for c in u.columns if c not in ID_COLS]].fillna(0).astype(int).reset_index(drop=True)
    splits = splits.set_index("isolate").loc[order].reset_index(drop=True)
    base = base.set_index("isolate").loc[order].reset_index(drop=True)
    sel = np.where(primary)[0]
    return (meta.iloc[sel].reset_index(drop=True), F.iloc[sel].reset_index(drop=True), U.iloc[sel].reset_index(drop=True),
            mapping, splits.iloc[sel].reset_index(drop=True), base.iloc[sel].reset_index(drop=True))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default="data/processed/metadata_curated.tsv")
    ap.add_argument("--features", default="data/processed/features_presence.tsv")
    ap.add_argument("--uncertain", default="data/processed/features_uncertain.tsv")
    ap.add_argument("--splits", default="data/processed/splits.tsv")
    ap.add_argument("--baseline", default="data/processed/baseline_predictions.tsv")
    ap.add_argument("--mapping", default="config/element_classes.tsv")
    ap.add_argument("--out-dir", default="data/processed")
    ap.add_argument("--report", default="reports/ml_report.md")
    ap.add_argument("--fast", action="store_true", help="small models and few resamples (testing only)")
    a = ap.parse_args(argv)
    if not Path(a.mapping).exists():
        print(f"Mapping not found: {a.mapping}. Run 05_map_elements_to_classes.py first.")
        return 1
    meta, F, U, mapping, splits, base = load_inputs(a)
    n_boot = 100 if a.fast else N_BOOT
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        preds, metrics, importance = run_all(meta, F, U, mapping, splits, base, fast=a.fast, n_boot=n_boot)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    suffix = "_fast" if a.fast else ""
    metrics.to_csv(out / f"ml_metrics{suffix}.tsv", sep="\t", index=False, lineterminator="\n")
    importance.to_csv(out / f"ml_feature_importance{suffix}.tsv", sep="\t", index=False, lineterminator="\n")
    preds.to_csv(out / f"ml_predictions{suffix}.tsv.gz", sep="\t", index=False)
    write_report(a.report, metrics, importance, a.fast)
    print(f"Models evaluated: {metrics['model'].nunique() - 1} variants x {metrics['drug'].nunique()} antibiotics x {metrics['split'].nunique()} splits")
    print("Report:", a.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
