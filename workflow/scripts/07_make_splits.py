#!/usr/bin/env python3
"""
07_make_splits.py

Define the evaluation splits used for every model, written ONCE here so that the baseline and all machine-learning
models are judged on exactly the same isolates. Only the primary dataset is split.

  fold5              5-fold cross-validation where whole SNP-cluster groups stay together (no near-identical isolate
                     in both train and test). Measures performance on NEW lineages.
  transfer_naive     train = all non-clinical isolates, test = clinical isolates. LEAKY: many clinical isolates have
                     close relatives among the non-clinical isolates. Kept only to show the size of the leak.
  transfer_disjoint  as above, but non-clinical isolates from any group that also contains a clinical isolate are
                     excluded from training. The clean transfer test.
  time               train = collected up to CUTOFF_YEAR, test = later. Measures forward prediction; shared lineages
                     across years are expected here and are reported.

Reads : data/processed/metadata_curated.tsv
Writes: data/processed/splits.tsv, reports/splits_report.md

Run from the repository root:
    python workflow/scripts/07_make_splits.py
"""
import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

K_FOLDS = 5
SEED = 42
CUTOFF_YEAR = 2016
DRUGS = ["tetracycline", "ampicillin", "ceftriaxone", "ciprofloxacin"]
MIN_R_PER_FOLD = 10


def assign_group_folds(group_sizes, k=K_FOLDS, seed=SEED):
    """Largest groups first, each to the currently smallest fold (ties broken by a seeded shuffle)."""
    order = group_sizes.sample(frac=1.0, random_state=seed).sort_values(ascending=False, kind="stable")
    loads = np.zeros(k)
    assign = {}
    for g, size in order.items():
        f = int(np.argmin(loads))
        assign[g] = f
        loads[f] += size
    return assign


def make_splits(meta, k=K_FOLDS, seed=SEED, cutoff=CUTOFF_YEAR):
    meta = meta.reset_index(drop=True)
    primary = meta["in_primary_set"].astype(str).str.lower().eq("true")
    out = meta[["isolate", "in_primary_set", "host_group", "cluster_group", "collection_year"]].copy()

    p = meta[primary]
    folds = assign_group_folds(p["cluster_group"].value_counts(), k, seed)
    out["fold5"] = pd.array(np.where(primary, meta["cluster_group"].map(folds), np.nan), dtype="Int64")

    clin = meta["host_group"] == "clinical"
    out["transfer_naive"] = np.where(primary, np.where(clin, "test", "train"), None)

    shared = set(p.groupby("cluster_group")["host_group"].nunique().loc[lambda s: s > 1].index)
    excluded = primary & ~clin & meta["cluster_group"].isin(shared)
    td = np.where(clin, "test", np.where(excluded, "excluded", "train"))
    out["transfer_disjoint"] = np.where(primary, td, None)

    year = pd.to_numeric(meta["collection_year"], errors="coerce").astype(float)   # NA -> NaN (nullable Int64 breaks comparisons)
    tm = np.where(year.isna(), None, np.where(year <= cutoff, "train", "test"))
    out["time"] = np.where(primary, tm, None)
    return out, shared


def leakage(out, split_col, labels):
    """Share of test isolates whose cluster group also occurs in train."""
    tr = set(out.loc[out[split_col].isin(labels[0]), "cluster_group"])
    te = out[out[split_col].isin(labels[1])]
    return len(te), int(te["cluster_group"].isin(tr).sum())


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def rs(series):
    return f"{int((series == 'R').sum())} / {int((series == 'S').sum())}"


def write_report(path, meta, out, shared):
    m = meta.reset_index(drop=True)
    L = ["# Evaluation splits report\n", f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}  \nScript: `workflow/scripts/07_make_splits.py`  \n"
         f"Seed {SEED}; folds {K_FOLDS}; time cutoff {CUTOFF_YEAR}. Only the primary dataset is split.\n"]
    # sizes
    rows = []
    for col in ("transfer_naive", "transfer_disjoint", "time"):
        vc = out[col].value_counts()
        rows.append([col, int(vc.get("train", 0)), int(vc.get("test", 0)), int(vc.get("excluded", 0))])
    f = out["fold5"].dropna().astype(int).value_counts().sort_index()
    rows.append(["fold5 (each fold is the test set once)", "see folds table", f"{int(f.min())}-{int(f.max())} per fold", 0])
    L.append("## Sizes\n\n" + md_table(pd.DataFrame(rows, columns=["Split", "Train", "Test", "Excluded"])) + "\n")
    # leakage
    rows = []
    for col, labs in (("transfer_naive", (["train"], ["test"])), ("transfer_disjoint", (["train"], ["test"])), ("time", (["train"], ["test"]))):
        n, k = leakage(out, col, labs)
        rows.append([col, n, k, f"{100 * k / n:.1f}%" if n else "NA"])
    n_t = k_t = 0
    for fold in range(K_FOLDS):
        te = out[out["fold5"] == fold]
        tr = set(out.loc[out["fold5"].notna() & (out["fold5"] != fold), "cluster_group"])
        n_t += len(te); k_t += int(te["cluster_group"].isin(tr).sum())
    rows.append(["fold5 (all folds)", n_t, k_t, f"{100 * k_t / n_t:.1f}%" if n_t else "NA"])
    L.append("## Leakage check: test isolates whose SNP-cluster group also appears in training\n\n"
             + md_table(pd.DataFrame(rows, columns=["Split", "Test isolates", "With a relative in train", "Share"])) + "\n")
    L.append("`fold5` and `transfer_disjoint` must show 0. `transfer_naive` shows how leaky the simple human-vs-other split would be. "
             "`time` is expected to share lineages across years; that is a realistic forward-prediction setting.\n")
    # labels
    prim = m[m["in_primary_set"].astype(str).str.lower() == "true"]
    L.append("## Resistant / susceptible isolates in train and test\n")
    rows = []
    for d in DRUGS:
        for col in ("transfer_disjoint", "time"):
            tr = prim.loc[out.loc[prim.index, col] == "train", f"ast_{d}"]
            te = prim.loc[out.loc[prim.index, col] == "test", f"ast_{d}"]
            rows.append([d, col, rs(tr), rs(te)])
    L.append(md_table(pd.DataFrame(rows, columns=["Antibiotic", "Split", "Train R / S", "Test R / S"])) + "\n")
    # folds
    rows = []
    warn = []
    for fold in range(K_FOLDS):
        idx = out.index[out["fold5"] == fold]
        r = [fold, len(idx), prim.loc[prim.index.intersection(idx), "cluster_group"].nunique()]
        for d in DRUGS:
            n_r = int((m.loc[idx, f"ast_{d}"] == "R").sum())
            r.append(n_r)
            if d != "ciprofloxacin" and n_r < MIN_R_PER_FOLD:
                warn.append(f"fold {fold} has only {n_r} resistant isolates for {d}")
        rows.append(r)
    L.append("## Folds (grouped 5-fold)\n\n" + md_table(pd.DataFrame(rows, columns=["Fold", "Isolates", "Groups"] + [f"R {d}" for d in DRUGS])) + "\n")
    L.append(("WARNING: " + "; ".join(warn) + "\n") if warn else f"Every fold has at least {MIN_R_PER_FOLD} resistant isolates for each core antibiotic.\n")
    # transfer details
    clin = prim[prim["host_group"] == "clinical"]
    L.append("## Why the clean transfer test matters\n")
    L.append(f"- Groups containing both clinical and non-clinical isolates: {len(shared)}\n"
             f"- Clinical isolates in such shared groups: {int(clin['cluster_group'].isin(shared).sum())} of {len(clin)}\n"
             f"- Non-clinical isolates dropped from training in the clean version: {int((out['transfer_disjoint'] == 'excluded').sum())}\n")
    yr = pd.to_numeric(clin["collection_year"], errors="coerce").dropna().astype(int)
    if len(yr):
        top_y = yr.value_counts().head(3)
        L.append("## Concentration of the clinical test set\n")
        L.append("- Clinical isolates by collection year (top 3): " + ", ".join(f"{y}: {n}" for y, n in top_y.items()) + f" (of {len(yr)} with a year)\n")
        if "bioproject" in clin.columns:
            bp = clin["bioproject"].value_counts().head(3)
            L.append("- Clinical isolates by BioProject (top 3): " + ", ".join(f"{b}: {n}" for b, n in bp.items()) + f" (of {len(clin)})\n")
        L.append("- If one year or study dominates, the transfer test mostly measures transfer to that collection, not to human isolates in general.\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(L), encoding="utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default="data/processed/metadata_curated.tsv")
    ap.add_argument("--out", default="data/processed/splits.tsv")
    ap.add_argument("--report", default="reports/splits_report.md")
    a = ap.parse_args(argv)
    meta = pd.read_csv(a.metadata, sep="\t", dtype={"collection_year": "Int64"})
    out, shared = make_splits(meta)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False, lineterminator="\n")
    write_report(a.report, meta, out, shared)
    n_leak = leakage(out, "transfer_disjoint", (["train"], ["test"]))[1]
    print(f"Splits written for {int(out['fold5'].notna().sum())} primary isolates | transfer_disjoint leakage: {n_leak} | time cutoff {CUTOFF_YEAR}")
    print("Report:", a.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
