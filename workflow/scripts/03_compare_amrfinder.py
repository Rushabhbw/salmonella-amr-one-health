#!/usr/bin/env python3
"""
03_compare_amrfinder.py

Compare AMRFinderPlus results we ran ourselves with the results NCBI stored for the same assemblies.

Reads : data/interim/amrfinder_local/*.amrfinder.tsv   (made by run_amrfinder_subset.sh)
        data/processed/metadata_curated.tsv            (NCBI's 'AMR genotypes core' string per isolate)
Writes: reports/amrfinder_concordance.md
        data/processed/amrfinder_concordance_per_isolate.tsv

What is compared: the SET of resistance elements (gene names and point mutations) per isolate.
Copy number and position are ignored (NCBI lists each element once). Only core AMR elements are used.

Run from the repository root:
    python workflow/scripts/03_compare_amrfinder.py
"""
import argparse
import re
from datetime import datetime
from pathlib import Path

import pandas as pd

PAIR = re.compile(r"([^=,]+)=([A-Za-z_]+)")


def parse_ncbi_core(value):
    """'blaTEM-1=COMPLETE,gyrA_S83F=POINT' -> {'blaTEM-1': 'COMPLETE', 'gyrA_S83F': 'POINT'}"""
    if not isinstance(value, str):
        return {}
    return {k.strip(): v for k, v in PAIR.findall(value)}


def read_local(path):
    """AMRFinderPlus output table -> {element symbol: set of methods}. Core AMR rows only."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        return {}
    if "Type" in df.columns:
        df = df[df["Type"] == "AMR"]
    if "Scope" in df.columns:
        df = df[df["Scope"] == "core"]
    out = {}
    for sym, meth in zip(df["Element symbol"], df["Method"]):
        out.setdefault(sym, set()).add(meth)
    return out


def compare_sets(local, ncbi):
    a, b = set(local), set(ncbi)
    shared, only_local, only_ncbi = a & b, a - b, b - a
    union = a | b
    return {
        "n_local": len(a), "n_ncbi": len(b), "n_shared": len(shared),
        "only_local": sorted(only_local), "only_ncbi": sorted(only_ncbi),
        "jaccard": (len(shared) / len(union)) if union else 1.0,
        "identical": a == b,
    }


def local_db_version(versions_path):
    p = Path(versions_path)
    if not p.exists():
        return None
    m = re.search(r"Database version:\s*(\S+)", p.read_text(encoding="utf-8", errors="replace"))
    return m.group(1) if m else None


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def run(amr_dir, metadata_path):
    meta = pd.read_csv(metadata_path, sep="\t", dtype=str).drop_duplicates("assembly").set_index("assembly")
    rows, unmatched, discord = [], [], []
    for f in sorted(Path(amr_dir).glob("*.amrfinder.tsv")):
        acc = f.name[: -len(".amrfinder.tsv")]
        if acc not in meta.index:
            unmatched.append(acc)
            continue
        m = meta.loc[acc]
        local = read_local(f)
        ncbi = parse_ncbi_core(m.get("amr_genes_core"))
        c = compare_sets(local, ncbi)
        rows.append({
            "assembly": acc, "isolate": m.get("isolate"), "host_group": m.get("host_group"),
            "serotype": m.get("serotype"),
            "ncbi_software": m.get("amrfinder_version"), "ncbi_catalog": m.get("amrfinder_catalog_version"),
            "n_ncbi": c["n_ncbi"], "n_local": c["n_local"], "n_shared": c["n_shared"],
            "jaccard": round(c["jaccard"], 3), "identical": c["identical"],
            "only_ncbi": ";".join(c["only_ncbi"]), "only_local": ";".join(c["only_local"]),
        })
        for sym in c["only_ncbi"]:
            discord.append({"side": "only in NCBI result", "element": sym, "call_or_method": ncbi[sym]})
        for sym in c["only_local"]:
            discord.append({"side": "only in our result", "element": sym, "call_or_method": "/".join(sorted(local[sym]))})
    return pd.DataFrame(rows), pd.DataFrame(discord), unmatched


def write_report(path, per, disc, unmatched, local_db, versions_text):
    L = ["# AMRFinderPlus concordance: our run vs NCBI's stored result\n"]
    L.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}  \nScript: `workflow/scripts/03_compare_amrfinder.py`\n")
    L.append("## Our run\n")
    L.append("```\n" + (versions_text.strip() if versions_text else "versions.txt not found") + "\n```\n")
    L.append("## Overall\n")
    if per.empty:
        L.append("No result tables could be matched to the curated metadata.\n")
    else:
        n, ident = len(per), int(per["identical"].sum())
        L.append(f"- Isolates compared: {n}\n- Identical sets of resistance elements: {ident} ({100 * ident / n:.1f}%)\n"
                 f"- Mean Jaccard similarity: {per['jaccard'].mean():.3f} (1.000 means identical)\n"
                 f"- Isolates with elements only in NCBI's result: {int((per['only_ncbi'] != '').sum())}\n"
                 f"- Isolates with elements only in our result: {int((per['only_local'] != '').sum())}\n")
        L.append("## By NCBI reference-gene-catalog version\n")
        per = per.copy()
        per["same_db_as_ours"] = per["ncbi_catalog"].apply(lambda v: "yes" if (local_db and v == local_db) else "no")
        g = per.groupby(["ncbi_catalog", "same_db_as_ours"]).agg(isolates=("assembly", "count"), identical=("identical", "sum")).reset_index()
        L.append(md_table(g) + "\n")
        L.append("## Differences\n")
        bad = per[~per["identical"]][["isolate", "host_group", "ncbi_catalog", "only_ncbi", "only_local"]]
        L.append(md_table(bad) + "\n" if len(bad) else "None. Every compared isolate has the same set of elements.\n")
        if len(disc):
            t = (disc.assign(call_or_method=disc["call_or_method"].astype(str))
                     .groupby(["side", "element", "call_or_method"]).size().reset_index(name="isolates")
                     .sort_values("isolates", ascending=False))
            L.append("Differing elements (what the call or method was on the side that reported it):\n\n" + md_table(t) + "\n")
    if unmatched:
        L.append("## Result files that did not match the curated metadata\n\n" + ", ".join(unmatched) + "\n")
    L.append("## How to read this\n")
    L.append("- A difference can come from a different database version, a renamed element, or a PARTIAL call that one run reports and the other does not. "
             "Check the call or method column above before assuming an error.\n"
             "- Only element names are compared, not copy number, position or call quality.\n"
             "- This is a sample, so it supports reproducibility checking only. It does not measure prediction accuracy.\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(L), encoding="utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--amr-dir", default="data/interim/amrfinder_local")
    ap.add_argument("--metadata", default="data/processed/metadata_curated.tsv")
    ap.add_argument("--report", default="reports/amrfinder_concordance.md")
    ap.add_argument("--table", default="data/processed/amrfinder_concordance_per_isolate.tsv")
    a = ap.parse_args(argv)

    per, disc, unmatched = run(a.amr_dir, a.metadata)
    vpath = Path(a.amr_dir) / "versions.txt"
    vtext = vpath.read_text(encoding="utf-8", errors="replace") if vpath.exists() else None
    if vtext:  # keep local folder names (user name, paths) out of a report that is committed to GitHub
        vtext = "\n".join(line for line in vtext.splitlines() if "directory:" not in line.lower())
    write_report(a.report, per, disc, unmatched, local_db_version(vpath), vtext)
    Path(a.table).parent.mkdir(parents=True, exist_ok=True)
    per.to_csv(a.table, sep="\t", index=False, lineterminator="\n")
    if per.empty:
        print("Nothing compared. Check that", a.amr_dir, "contains *.amrfinder.tsv files.")
        return 1
    print(f"Compared {len(per)} isolates | identical: {int(per['identical'].sum())} | mean Jaccard: {per['jaccard'].mean():.3f}")
    print("Report:", a.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
