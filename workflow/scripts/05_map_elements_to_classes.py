#!/usr/bin/env python3
"""
05_map_elements_to_classes.py

Link each resistance element in our feature table to the antibiotic class NCBI assigns to it,
using the AMRFinderPlus Reference Gene Catalog, and derive which elements count as evidence for
resistance to each of our antibiotics (the rules of the rule-based baseline).

Reads : data/external/ReferenceGeneCatalog.txt   (copy it from the AMRFinderPlus database folder; see fetch_catalog.sh)
        data/processed/feature_dictionary.tsv     (made by 04_build_features.py)
Writes: config/element_classes.tsv                (small, commit this file)
        reports/element_class_report.md           (for human review: unmapped elements, class counts, mapped lists)

Matching: an element is looked up in the catalog column `allele` first (e.g. blaTEM-1, gyrA_D87Y), then in `gene_family`.
Drug rules (documented in the report too; they are our reading of NCBI's class and subclass, not facts about each isolate):
    tetracycline  : class or subclass contains TETRACYCLINE
    ampicillin    : class contains BETA-LACTAM (any beta-lactamase)
    ceftriaxone   : class contains BETA-LACTAM and subclass contains CEPHALOSPORIN
    ciprofloxacin : class or subclass contains QUINOLONE
Two refinements (both written in the report):
  * An element counts only if EVERY catalog entry it matches satisfies the rule. A gene-family symbol such as blaTEM stands for
    many alleles (TEM-1 does not break cephalosporins, some other TEM alleles do), so it is not evidence for ceftriaxone.
  * Elements in NCBI's MULTIDRUG class (e.g. ramR mutations, low-level efflux changes) are not drug-specific evidence. They are
    kept in a separate column `<drug>_with_multidrug` for a sensitivity analysis.
NCBI warns that a gene encoding a resistance protein does not always mean the isolate is resistant.

Run from the repository root:
    python workflow/scripts/05_map_elements_to_classes.py
"""
import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd

REQUIRED = ["allele", "gene_family", "class", "subclass", "type", "subtype", "scope"]
DRUG_RULES = {
    "tetracycline": lambda c, s: "TETRACYCLINE" in c or "TETRACYCLINE" in s,
    "ampicillin": lambda c, s: "BETA-LACTAM" in c or "AMPICILLIN" in s,
    "ceftriaxone": lambda c, s: "BETA-LACTAM" in c and "CEPHALOSPORIN" in s,
    "ciprofloxacin": lambda c, s: "QUINOLONE" in c or "QUINOLONE" in s,
}


def read_catalog(path):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise SystemExit(f"Catalog is missing expected columns {missing}. Found: {list(df.columns)}")
    for c in ("class", "subclass", "type", "subtype", "scope"):
        df[c] = df[c].str.strip()
    return df


def build_index(catalog):
    """symbol -> catalog rows, separately for the allele and gene_family columns."""
    by_allele, by_family = {}, {}
    for col, store in (("allele", by_allele), ("gene_family", by_family)):
        for sym, rows in catalog[catalog[col].str.strip() != ""].groupby(catalog[col].str.strip()):
            store[sym] = rows
    return by_allele, by_family


def lookup(element, by_allele, by_family):
    for how, store in (("allele", by_allele), ("gene_family", by_family)):
        if element in store:
            rows = store[element]
            core = rows[rows["scope"].str.lower() == "core"]
            return how, (core if len(core) else rows)
    return "none", None


def drugs_for(pairs):
    """pairs: iterable of (class, subclass) in upper case -> {drug: 0/1, drug_with_multidrug: 0/1}.

    An element counts only if EVERY catalog entry satisfies the rule (a family-level symbol is not evidence if some
    of its alleles do not satisfy it). The plain column excludes MULTIDRUG-class entries; the `_with_multidrug`
    column allows them (sensitivity analysis)."""
    pairs = list(pairs)
    out = {}
    for d, rule in DRUG_RULES.items():
        out[d] = int(bool(pairs) and all(rule(c, s) and "MULTIDRUG" not in c for c, s in pairs))
        out[f"{d}_with_multidrug"] = int(bool(pairs) and all(rule(c, s) for c, s in pairs))
    return out


def map_elements(dictionary, catalog, extra_elements=()):
    """dictionary: feature dictionary (elements with strong calls). extra_elements: elements seen only with weak calls."""
    known = set(dictionary["element"])
    extra = [e for e in extra_elements if e not in known]
    if extra:
        dictionary = pd.concat([dictionary, pd.DataFrame({"element": extra, "kind": "weak_only"})], ignore_index=True)
    by_allele, by_family = build_index(catalog)
    out = []
    for _, r in dictionary.iterrows():
        how, rows = lookup(r["element"], by_allele, by_family)
        if rows is None:
            pairs, cls, sub, typ, subt, scope, dbv = [], "", "", "", "", "", ""
        else:
            pairs = sorted({(c.upper(), s.upper()) for c, s in zip(rows["class"], rows["subclass"])})
            cls = "|".join(sorted({c for c, _ in pairs}))
            sub = "|".join(sorted({s for _, s in pairs}))
            typ = "|".join(sorted(set(rows["type"])))
            subt = "|".join(sorted(set(rows["subtype"])))
            scope = "|".join(sorted(set(rows["scope"])))
            dbv = "|".join(sorted(set(rows["db_version"]))) if "db_version" in rows.columns else ""
        row = {"element": r["element"], "matched_on": how, "class": cls, "subclass": sub, "type": typ,
               "subtype": subt, "scope": scope, "ambiguous": len(pairs) > 1, "catalog_db_version": dbv}
        row.update(drugs_for(pairs))
        for k in ("kind", "n_isolates_primary", "n_groups_primary"):
            if k in dictionary.columns:
                row[k] = r[k]
        out.append(row)
    return pd.DataFrame(out)


def whole_numbers(series, blank=""):
    """Format counts as whole numbers without depending on the pandas version's handling of missing values."""
    return [blank if pd.isna(v) else str(int(v)) for v in pd.to_numeric(series, errors="coerce")]


def clean_counts(df, col="n_isolates_primary", blank="weak calls only"):
    """Show counts as whole numbers; elements that only ever had weak calls have no count."""
    df = df.copy()
    if col in df.columns:
        df[col] = whole_numbers(df[col], blank)
    return df


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def write_report(path, mapping, catalog_path):
    L = ["# Element-to-class mapping report\n",
         f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}  \nScript: `workflow/scripts/05_map_elements_to_classes.py`  \nCatalog file: `{Path(catalog_path).name}`\n"]
    dbv = sorted({v for x in mapping["catalog_db_version"] for v in str(x).split("|") if v})
    L.append(f"Catalog database version(s) seen for matched elements: {', '.join(dbv) if dbv else 'not stated in the file'}\n")
    L.append("## Rules\n")
    L.append("- tetracycline: class or subclass contains TETRACYCLINE\n- ampicillin: class contains BETA-LACTAM (any beta-lactamase)\n"
             "- ceftriaxone: class contains BETA-LACTAM and subclass contains CEPHALOSPORIN (NCBI's definition of CEPHALOSPORIN is the Lahey 2be group; CARBAPENEM subclass is not counted)\n"
             "- ciprofloxacin: class or subclass contains QUINOLONE\n"
             "- An element counts only if every catalog entry it matches satisfies the rule (a family-level symbol such as blaTEM is not evidence for ceftriaxone, because TEM-1 does not break cephalosporins).\n"
             "- MULTIDRUG-class elements are not drug-specific evidence; they are listed separately and used only in a sensitivity analysis.\n"
             "- NCBI cautions that carrying a resistance gene does not always mean the isolate is resistant.\n")
    n = len(mapping)
    unm = mapping[mapping["matched_on"] == "none"]
    L.append("## Matching\n")
    L.append(f"- Elements: {n}\n- Matched on allele: {int((mapping['matched_on'] == 'allele').sum())}\n- Matched on gene family: {int((mapping['matched_on'] == 'gene_family').sum())}\n"
             f"- Not found in the catalog: {len(unm)}\n- With more than one class/subclass pair (ambiguous): {int(mapping['ambiguous'].sum())}\n")
    if len(unm):
        cols = [c for c in ["element", "kind", "n_isolates_primary"] if c in unm.columns]
        L.append("Not found (check these by hand):\n\n" + md_table(clean_counts(unm[cols]).head(30)) + "\n")
    L.append("## Class and subclass of the matched elements\n")
    t = (mapping[mapping["matched_on"] != "none"].groupby(["class", "subclass"]).size().reset_index(name="elements")
         .sort_values("elements", ascending=False).head(25))
    L.append(md_table(t) + "\n")
    L.append("## Elements that count as evidence for each antibiotic\n")
    for d in DRUG_RULES:
        md_only = mapping[(mapping[f"{d}_with_multidrug"] == 1) & (mapping[d] == 0)]
        sel = mapping[mapping[d] == 1]
        L.append(f"### {d} ({len(sel)} elements)\n")
        if len(sel) == 0:
            L.append("WARNING: no element maps to this antibiotic. Check the catalog columns and the rule.\n")
            continue
        cols = [c for c in ["element", "kind", "class", "subclass", "n_isolates_primary", "n_groups_primary"] if c in sel.columns]
        s = sel[cols].copy()
        if "n_isolates_primary" in s.columns:
            s["n_isolates_primary"] = pd.to_numeric(s["n_isolates_primary"], errors="coerce")
            s = s.sort_values("n_isolates_primary", ascending=False)
        if "n_groups_primary" in s.columns:
            s["n_groups_primary"] = whole_numbers(s["n_groups_primary"], "")
        L.append(md_table(clean_counts(s).head(40)) + "\n")
        if len(md_only):
            L.append("Excluded as MULTIDRUG-class or not drug-specific (sensitivity analysis only): " + ", ".join(md_only["element"].astype(str)) + "\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(L), encoding="utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", default="data/external/ReferenceGeneCatalog.txt")
    ap.add_argument("--dictionary", default="data/processed/feature_dictionary.tsv")
    ap.add_argument("--uncertain", default="data/processed/features_uncertain.tsv",
                    help="optional: elements seen only with weak calls are added to the mapping")
    ap.add_argument("--out", default="config/element_classes.tsv")
    ap.add_argument("--report", default="reports/element_class_report.md")
    a = ap.parse_args(argv)
    if not Path(a.catalog).exists():
        print(f"Catalog not found: {a.catalog}. Run (in Ubuntu): bash workflow/scripts/fetch_catalog.sh")
        return 1
    catalog = read_catalog(a.catalog)
    dictionary = pd.read_csv(a.dictionary, sep="\t", dtype=str)
    extra = []
    if Path(a.uncertain).exists():
        extra = [c for c in pd.read_csv(a.uncertain, sep="\t", nrows=0).columns
                 if c not in ("isolate", "assembly", "in_primary_set", "host_group", "cluster_group")]
    mapping = map_elements(dictionary, catalog, extra)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(a.out, sep="\t", index=False, lineterminator="\n")
    write_report(a.report, mapping, a.catalog)
    print(f"Elements: {len(mapping)} | matched: {int((mapping['matched_on'] != 'none').sum())} | not found: {int((mapping['matched_on'] == 'none').sum())}")
    n_added = int((mapping["kind"] == "weak_only").sum()) if "kind" in mapping.columns else 0
    print("Elements per antibiotic:", {d: int(mapping[d].sum()) for d in DRUG_RULES}, "| elements added because they only had weak calls:", n_added)
    print("Report:", a.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
