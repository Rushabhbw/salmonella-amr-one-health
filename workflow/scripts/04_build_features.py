#!/usr/bin/env python3
"""
04_build_features.py

Turn NCBI's stored AMRFinderPlus results (column amr_genes_core) into a feature table:
one row per isolate, one 0/1 column per resistance gene or point mutation.

Reads : data/processed/metadata_curated.tsv
Writes: data/processed/features_presence.tsv    1 = element present with a COMPLETE or POINT call
        data/processed/features_uncertain.tsv   1 = element seen ONLY with a weak call (PARTIAL, PARTIAL_END_OF_CONTIG,
                                                MISTRANSLATION, HMM) and never with a COMPLETE/POINT call (sensitivity analyses)
        data/processed/feature_dictionary.tsv   one row per element: how many isolates and independent groups carry it
        reports/feature_report.md               counts and sanity checks

Rules:
  * Presence is a yes/no per isolate. Copy number is ignored (NCBI lists each element once and our local run showed duplicates).
  * No feature is dropped here. Rare features are kept so that any filtering can be done inside cross-validation, not on all data.
  * An empty result string means "analysed, no resistance element found" (every curated isolate has an AMRFinderPlus version).

Run from the repository root:
    python workflow/scripts/04_build_features.py
"""
import argparse
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

PAIR = re.compile(r"([^=,]+)=([A-Za-z_]+)")
PRESENT_CALLS = {"COMPLETE", "POINT"}
UNCERTAIN_CALLS = {"PARTIAL", "PARTIAL_END_OF_CONTIG", "MISTRANSLATION", "HMM"}
ID_COLUMNS = ["isolate", "assembly", "in_primary_set", "host_group", "cluster_group"]
DRUGS = ["tetracycline", "ampicillin", "ceftriaxone", "ciprofloxacin"]


def parse_calls(value):
    """'blaTEM-1=COMPLETE,tet(A)=PARTIAL' -> [('blaTEM-1', 'COMPLETE'), ('tet(A)', 'PARTIAL')]"""
    if not isinstance(value, str) or not value.strip():
        return []
    return [(k.strip(), v) for k, v in PAIR.findall(value)]


def build_matrices(meta):
    """Return (present, uncertain, dictionary_rows_info, unknown_calls) for the isolates in meta (row order kept)."""
    present_sets, uncertain_sets, call_types, unknown = [], [], {}, {}
    for value in meta["amr_genes_core"]:
        pres, unc = set(), set()
        for el, call in parse_calls(value):
            call_types.setdefault(el, set()).add(call)
            if call in PRESENT_CALLS:
                pres.add(el)
            elif call in UNCERTAIN_CALLS:
                unc.add(el)
            else:                       # a call type we have not seen before: keep it weak and report it
                unc.add(el)
                unknown[call] = unknown.get(call, 0) + 1
        present_sets.append(pres)
        uncertain_sets.append(unc - pres)   # a weak call does not count if the element also has a strong call
    all_elements = sorted(set().union(*present_sets, *uncertain_sets)) if present_sets else []

    def to_matrix(sets):
        mat = np.zeros((len(sets), len(all_elements)), dtype=np.int8)
        idx = {e: i for i, e in enumerate(all_elements)}
        for r, s in enumerate(sets):
            for e in s:
                mat[r, idx[e]] = 1
        return pd.DataFrame(mat, columns=all_elements, index=meta.index)

    present, uncertain = to_matrix(present_sets), to_matrix(uncertain_sets)
    # keep only elements that have at least one COMPLETE/POINT call; order by how many isolates carry them, then by name
    strong = [e for e in all_elements if present[e].sum() > 0]
    order = sorted(strong, key=lambda e: (-int(present[e].sum()), e))
    return present[order], uncertain, call_types, unknown


def feature_dictionary(meta, present, uncertain, call_types):
    rows = []
    primary = meta["in_primary_set"].astype(str).str.lower().eq("true")
    for el in present.columns:
        has = present[el].astype(bool)
        rows.append({
            "element": el,
            "kind": "point_mutation" if "POINT" in call_types.get(el, set()) else "gene",
            "n_isolates_all": int(has.sum()),
            "n_isolates_primary": int((has & primary).sum()),
            "n_groups_all": int(meta.loc[has, "cluster_group"].nunique()),
            "n_groups_primary": int(meta.loc[has & primary, "cluster_group"].nunique()),
            "n_clinical": int((has & (meta["host_group"] == "clinical")).sum()),
            "n_non_clinical": int((has & (meta["host_group"] == "non_clinical")).sum()),
            "n_uncertain_only": int(uncertain[el].sum()) if el in uncertain.columns else 0,
            "call_types_seen": "/".join(sorted(call_types.get(el, set()))),
        })
    return pd.DataFrame(rows)


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def write_report(path, meta, present, uncertain, dic, call_types, unknown):
    n = len(meta)
    n_el = present.shape[1]
    n_weak_only = int(sum(1 for c in uncertain.columns if uncertain[c].sum() > 0 and c not in present.columns))
    count = present.sum(axis=1)
    zero = count == 0
    L = ["# Feature table report\n", f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}  \nScript: `workflow/scripts/04_build_features.py`\n"]
    L.append("## Rules\n")
    L.append("- Feature = a resistance gene or point mutation from NCBI's stored AMRFinderPlus core results.\n"
             "- 1 = present with a COMPLETE or POINT call. Weak calls (PARTIAL, PARTIAL_END_OF_CONTIG, MISTRANSLATION, HMM) are kept in a separate table and are not used as features.\n"
             "- Copy number is ignored. No feature is removed here; rare-feature filtering must be done inside cross-validation.\n")
    L.append("## Size\n")
    L.append(f"- Isolates: {n}\n- Elements with at least one COMPLETE/POINT call: {n_el}\n"
             f"- Elements seen only with weak calls (never COMPLETE/POINT; not features): {n_weak_only}\n"
             f"- Isolates with no element at all: {int(zero.sum())} ({100 * zero.mean():.1f}%)\n"
             f"- Elements per isolate (median / max): {int(count.median())} / {int(count.max())}\n")
    if unknown:
        L.append(f"- WARNING: unknown call types treated as weak: {unknown}\n")
    L.append("## Elements by how widely they occur (primary dataset)\n")
    rows = []
    for thr in (1, 3, 10, 20):
        rows.append([f">= {thr} independent groups", int((dic["n_groups_primary"] >= thr).sum())])
    for thr in (10, 30, 100):
        rows.append([f">= {thr} isolates", int((dic["n_isolates_primary"] >= thr).sum())])
    L.append(md_table(pd.DataFrame(rows, columns=["Elements present in", "Number of elements"])) + "\n")
    L.append("## Most common elements (primary dataset)\n")
    top = dic.sort_values("n_isolates_primary", ascending=False).head(15)[["element", "kind", "n_isolates_primary", "n_groups_primary", "n_clinical", "n_non_clinical"]]
    L.append(md_table(top) + "\n")
    L.append("## Sanity check: isolates with no element\n")
    L.append("If resistance is explained by known genes, isolates with no element should be mostly susceptible. "
             "Resistant isolates here have no known explanation in the data (unknown mechanism, label noise, or a gap in the database).\n")
    rows = []
    for d in DRUGS:
        lab = meta[f"y_{d}"]
        for name, mask in (("no element", zero), ("at least one element", ~zero)):
            g = lab[mask].dropna()
            rows.append([d, name, int((g == 1).sum()), int((g == 0).sum()), f"{100 * (g == 1).mean():.1f}%" if len(g) else "NA"])
    L.append(md_table(pd.DataFrame(rows, columns=["Antibiotic", "Isolates with", "R", "S", "R share"])) + "\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(L), encoding="utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default="data/processed/metadata_curated.tsv")
    ap.add_argument("--out-dir", default="data/processed")
    ap.add_argument("--report", default="reports/feature_report.md")
    a = ap.parse_args(argv)

    meta = pd.read_csv(a.metadata, sep="\t")
    for d in DRUGS:                      # keep labels numeric (empty = no label)
        meta[f"y_{d}"] = pd.to_numeric(meta[f"y_{d}"], errors="coerce")
    present, uncertain, call_types, unknown = build_matrices(meta)
    dic = feature_dictionary(meta, present, uncertain, call_types)

    ids = meta[ID_COLUMNS].copy()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pd.concat([ids, present], axis=1).to_csv(out / "features_presence.tsv", sep="\t", index=False, lineterminator="\n")
    only_weak = [c for c in uncertain.columns if uncertain[c].sum() > 0]
    pd.concat([ids, uncertain[only_weak]], axis=1).to_csv(out / "features_uncertain.tsv", sep="\t", index=False, lineterminator="\n")
    dic.to_csv(out / "feature_dictionary.tsv", sep="\t", index=False, lineterminator="\n")
    write_report(a.report, meta, present, uncertain, dic, call_types, unknown)
    print(f"Features: {present.shape[1]} elements x {len(meta)} isolates  |  isolates with no element: {int((present.sum(axis=1) == 0).sum())}")
    print("Report:", a.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
