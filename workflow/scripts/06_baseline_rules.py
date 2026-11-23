#!/usr/bin/env python3
"""
06_baseline_rules.py

Rule-based baseline: predict an isolate RESISTANT to an antibiotic if it carries at least one resistance
element that NCBI links to that antibiotic (config/element_classes.tsv), otherwise SUSCEPTIBLE.
No training is involved, so every labelled isolate can be used for evaluation.

Reads : data/processed/metadata_curated.tsv, data/processed/features_presence.tsv, config/element_classes.tsv,
        data/processed/features_uncertain.tsv (optional: partial calls, used for diagnostics and one sensitivity row)
Writes: data/processed/baseline_predictions.tsv
        reports/baseline_report.md

Metrics per antibiotic and group (primary dataset unless stated):
  sensitivity = TP/(TP+FN), specificity = TN/(TN+FP), MCC,
  very major error rate (VME) = FN/(TP+FN)   (resistant isolates predicted susceptible)
  major error rate (ME)       = FP/(TN+FP)   (susceptible isolates predicted resistant)
Extra rows (primary dataset): the baseline also counting MULTIDRUG-class elements, and also counting partial calls.
Diagnostics: elements behind the false alarms, resistant isolates with no explaining element, and which BioProjects
(studies) carry the most errors, because study-specific label problems would show up there.
Confidence intervals use a CLUSTER BOOTSTRAP: whole SNP-cluster groups are resampled with replacement,
because isolates in the same group are nearly identical and are not independent evidence.

Run from the repository root:
    python workflow/scripts/06_baseline_rules.py
"""
import argparse
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

DRUGS = ["tetracycline", "ampicillin", "ceftriaxone", "ciprofloxacin"]
N_BOOT = 2000
SEED = 42


# ----------------------------------------------------------------------------- metrics
def metrics_from_counts(c):
    """c: array-like with rows or a single vector of [TP, FN, TN, FP]; returns dict of arrays/floats."""
    c = np.asarray(c, dtype=float)
    tp, fn, tn, fp = c[..., 0], c[..., 1], c[..., 2], c[..., 3]
    with np.errstate(divide="ignore", invalid="ignore"):
        sens, spec = tp / (tp + fn), tn / (tn + fp)
        ppv, npv = tp / (tp + fp), tn / (tn + fn)
        acc = (tp + tn) / (tp + fn + tn + fp)
        mcc = (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {"sensitivity": sens, "specificity": spec, "ppv": ppv, "npv": npv, "accuracy": acc,
            "balanced_accuracy": (sens + spec) / 2, "mcc": mcc, "vme": 1 - sens, "me": 1 - spec}


def group_counts(df, pred_col, y_col):
    """Per cluster group [TP, FN, TN, FP] counts, labelled isolates only."""
    d = df[df[y_col].notna()]
    y, p = d[y_col].astype(int), d[pred_col].astype(int)
    tmp = pd.DataFrame({"g": d["cluster_group"], "TP": ((y == 1) & (p == 1)).astype(int), "FN": ((y == 1) & (p == 0)).astype(int),
                        "TN": ((y == 0) & (p == 0)).astype(int), "FP": ((y == 0) & (p == 1)).astype(int)})
    return tmp.groupby("g")[["TP", "FN", "TN", "FP"]].sum().to_numpy()


def cluster_bootstrap(counts, n_boot=N_BOOT, seed=SEED, chunk=250):
    """Resample groups with replacement; returns {metric: (low, high)} (2.5 and 97.5 percentiles)."""
    g = len(counts)
    if g < 2:
        return {k: (np.nan, np.nan) for k in metrics_from_counts(counts.sum(axis=0) if g else [0, 0, 0, 0])}
    rng = np.random.default_rng(seed)
    sums = []
    for start in range(0, n_boot, chunk):
        n = min(chunk, n_boot - start)
        idx = rng.integers(0, g, size=(n, g))
        sums.append(counts[idx].sum(axis=1))
    m = metrics_from_counts(np.vstack(sums))
    out = {}
    for k, v in m.items():
        v = v[~np.isnan(v)]
        out[k] = (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) > 10 else (np.nan, np.nan)
    return out


def evaluate(df, pred_col, y_col):
    counts = group_counts(df, pred_col, y_col)
    tot = counts.sum(axis=0) if len(counts) else np.zeros(4)
    point = metrics_from_counts(tot)
    ci = cluster_bootstrap(counts)
    return {"counts": tot.astype(int), "n_groups": len(counts), "point": point, "ci": ci}


def fmt(v, ci=None, pct=False):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "NA"
    s = f"{100 * v:.1f}%" if pct else f"{v:.3f}"
    if ci is not None and not any(np.isnan(ci)):
        s += f" ({100 * ci[0]:.1f}-{100 * ci[1]:.1f}%)" if pct else f" ({ci[0]:.3f}-{ci[1]:.3f})"
    return s


# ----------------------------------------------------------------------------- prediction
ID_COLS = ("isolate", "assembly", "in_primary_set", "host_group", "cluster_group")


def evidence_elements(mapping, column, available):
    if column not in mapping.columns:
        return []
    return [e for e in mapping.loc[mapping[column].astype(int) == 1, "element"] if e in available]


def any_present(frame, elems, index):
    return frame[elems].max(axis=1).astype(int) if elems else pd.Series(0, index=index)


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def metric_row(name, df, pred_col, y_col):
    ev = evaluate(df, pred_col, y_col)
    tp, fn, tn, fp = ev["counts"]
    p, ci = ev["point"], ev["ci"]
    return [name, f"{tp + fn} / {tn + fp}", ev["n_groups"], f"{tp}/{fn}/{tn}/{fp}",
            fmt(p["sensitivity"], ci["sensitivity"]), fmt(p["specificity"], ci["specificity"]), fmt(p["mcc"], ci["mcc"]),
            fmt(p["vme"], ci["vme"], pct=True), fmt(p["me"], ci["me"], pct=True)]


def bioproject_errors(meta_p, err_mask, lab_mask, top=5):
    """Which BioProjects (studies) carry the most errors, and at what rate."""
    d = pd.DataFrame({"bioproject": meta_p["bioproject"].fillna("unknown"), "lab": lab_mask.values, "err": err_mask.values})
    g = d[d["lab"]].groupby("bioproject").agg(labelled=("lab", "sum"), errors=("err", "sum")).reset_index()
    total = int(g["errors"].sum())
    g["error_rate"] = (100 * g["errors"] / g["labelled"]).round(1).astype(str) + "%"
    g["share_of_all_errors"] = (100 * g["errors"] / total).round(1).astype(str) + "%" if total else "NA"
    g = g.sort_values(["errors", "labelled"], ascending=False).head(top)
    return g[g["errors"] > 0][["bioproject", "labelled", "errors", "error_rate", "share_of_all_errors"]]


def build(meta, feats, mapping, uncertain=None):
    """Returns (predictions DataFrame, list of per-drug section dicts)."""
    meta = meta.reset_index(drop=True)
    ordered = feats.set_index("isolate").loc[meta["isolate"]]
    element_cols = [c for c in ordered.columns if c not in ID_COLS]
    F = ordered[element_cols].astype(int).reset_index(drop=True)
    U = None
    if uncertain is not None:
        u = uncertain.set_index("isolate").reindex(meta["isolate"])
        U = u[[c for c in u.columns if c not in ID_COLS]].fillna(0).astype(int).reset_index(drop=True)
    preds = meta[["isolate", "in_primary_set", "host_group", "cluster_group"]].copy()
    primary = meta["in_primary_set"].astype(str).str.lower().eq("true")
    sections = []
    for d in DRUGS:
        y = pd.to_numeric(meta[f"y_{d}"], errors="coerce")
        preds[f"y_{d}"] = y
        elems = evidence_elements(mapping, d, F.columns)
        pred = any_present(F, elems, F.index)
        preds[f"pred_{d}"] = pred.values
        sub = preds.copy()
        md_elems = evidence_elements(mapping, f"{d}_with_multidrug", F.columns)
        sub[f"pred_{d}_md"] = any_present(F, md_elems, F.index).values
        weak = pd.Series(0, index=F.index)
        if U is not None:
            weak_elems = evidence_elements(mapping, d, U.columns)
            weak = any_present(U, weak_elems, U.index)
        sub[f"pred_{d}_partial"] = ((pred == 1) | (weak == 1)).astype(int).values

        rows = [metric_row("all (primary dataset)", sub[primary], f"pred_{d}", f"y_{d}"),
                metric_row("non-clinical (primary)", sub[primary & (meta["host_group"] == "non_clinical")], f"pred_{d}", f"y_{d}"),
                metric_row("clinical (primary)", sub[primary & (meta["host_group"] == "clinical")], f"pred_{d}", f"y_{d}"),
                metric_row("all isolates incl. non-primary (sensitivity check)", sub, f"pred_{d}", f"y_{d}"),
                metric_row("primary, also counting MULTIDRUG-class elements", sub[primary], f"pred_{d}_md", f"y_{d}")]
        if U is not None:
            rows.append(metric_row("primary, also counting partial calls", sub[primary], f"pred_{d}_partial", f"y_{d}"))
        tab = pd.DataFrame(rows, columns=["Group", "R / S isolates", "Groups", "TP/FN/TN/FP", "Sensitivity", "Specificity", "MCC",
                                          "Very major error (R called S)", "Major error (S called R)"])
        lab = primary & y.notna()
        fp_mask = lab & (y == 0) & (pred == 1)
        fn_mask = lab & (y == 1) & (pred == 0)
        fp_elems = Counter(e for i in np.where(fp_mask)[0] for e in elems if F.iloc[i][e] == 1)
        fn_none = int((fn_mask & (F.sum(axis=1) == 0)).sum())
        fn_weak = int((fn_mask & (weak == 1)).sum()) if U is not None else None
        bp = bioproject_errors(meta[primary], (fp_mask | fn_mask)[primary], lab[primary]) if "bioproject" in meta.columns else None
        sections.append({"drug": d, "elems": elems, "table": tab, "fp_top": fp_elems.most_common(6), "n_fp": int(fp_mask.sum()),
                         "n_fn": int(fn_mask.sum()), "fn_none": fn_none, "fn_other": int(fn_mask.sum()) - fn_none,
                         "fn_weak": fn_weak, "bioproject": bp, "error_rate": float((fp_mask | fn_mask).sum() / max(int(lab.sum()), 1))})
    return preds, sections


def write_report(path, sections):
    L = ["# Rule-based baseline report\n",
         f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}  \nScript: `workflow/scripts/06_baseline_rules.py`\n",
         "## Method\n",
         "- Prediction: RESISTANT if the isolate carries at least one element that NCBI links to the antibiotic (`config/element_classes.tsv`), otherwise SUSCEPTIBLE. No training.\n"
         "- Labels: S = 0, R = 1; intermediate and other calls are excluded for that antibiotic.\n"
         "- Very major error rate = resistant isolates predicted susceptible (FN / all R). Major error rate = susceptible isolates predicted resistant (FP / all S).\n"
         f"- Confidence intervals: cluster bootstrap, {N_BOOT} resamples of whole SNP-cluster groups, seed {SEED}. 'Groups' is the number of independent groups behind each row.\n"
         "- The last two rows of each table are sensitivity analyses on the primary dataset: they change the evidence rule, not the data.\n"
         "- Ciprofloxacin is descriptive only (see the charter): few independent resistant groups.\n"]
    for sec in sections:
        d = sec["drug"]
        L.append(f"## {d}\n")
        L.append(f"Elements used as evidence ({len(sec['elems'])}): {', '.join(sec['elems'][:25])}{' ...' if len(sec['elems']) > 25 else ''}\n")
        L.append(md_table(sec["table"]) + "\n")
        L.append(f"Primary dataset errors: {sec['n_fp']} major (S called R), {sec['n_fn']} very major (R called S).\n")
        if sec["fp_top"]:
            L.append("Elements most often present in susceptible isolates wrongly called resistant: " + ", ".join(f"{e} ({n})" for e, n in sec["fp_top"]) + "\n")
        txt = f"Resistant isolates called susceptible: {sec['fn_none']} carry no element at all; {sec['fn_other']} carry elements, but none linked to this antibiotic."
        if sec["fn_weak"] is not None:
            txt += f" {sec['fn_weak']} of the {sec['n_fn']} have a partial call for an element linked to this antibiotic (the gene may be broken across contigs)."
        L.append(txt + "\n")
        bp = sec["bioproject"]
        if bp is not None and len(bp):
            L.append(f"BioProjects with the most errors (overall error rate {100 * sec['error_rate']:.1f}%). A rate far above that points to a study-specific problem such as label noise or a different test method:\n\n" + md_table(bp) + "\n")
    L.append("## How to read this\n")
    L.append("- A strong baseline leaves little room for machine learning. A model must beat these numbers on the same isolates and splits to be worth reporting.\n"
             "- Wide intervals mean few independent groups, not necessarily a bad result.\n"
             "- NCBI cautions that carrying a resistance gene does not always mean the isolate is resistant, and some resistance has no known gene.\n"
             "- The evidence rules come from NCBI's catalog, not from these labels. Two rule choices were refined after reading a first version of this report (family-level symbols; MULTIDRUG-class elements); both alternatives are shown in the sensitivity rows.\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(L), encoding="utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default="data/processed/metadata_curated.tsv")
    ap.add_argument("--features", default="data/processed/features_presence.tsv")
    ap.add_argument("--uncertain", default="data/processed/features_uncertain.tsv")
    ap.add_argument("--mapping", default="config/element_classes.tsv")
    ap.add_argument("--out", default="data/processed/baseline_predictions.tsv")
    ap.add_argument("--report", default="reports/baseline_report.md")
    a = ap.parse_args(argv)
    if not Path(a.mapping).exists():
        print(f"Mapping not found: {a.mapping}. Run 05_map_elements_to_classes.py first.")
        return 1
    meta = pd.read_csv(a.metadata, sep="\t")
    feats = pd.read_csv(a.features, sep="\t")
    mapping = pd.read_csv(a.mapping, sep="\t")
    uncertain = pd.read_csv(a.uncertain, sep="\t") if Path(a.uncertain).exists() else None
    preds, sections = build(meta, feats, mapping, uncertain)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    preds.to_csv(a.out, sep="\t", index=False, lineterminator="\n")
    write_report(a.report, sections)
    print("Baseline written for:", ", ".join(DRUGS))
    for sec in sections:
        r = sec["table"].iloc[0]
        print(f"  {sec['drug']}: {len(sec['elems'])} elements; primary dataset sensitivity {r['Sensitivity']}, specificity {r['Specificity']}")
    print("Report:", a.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
