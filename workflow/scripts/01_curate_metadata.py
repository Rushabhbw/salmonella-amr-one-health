#!/usr/bin/env python3
"""
01_curate_metadata.py

Merge the NCBI Pathogen Detection Isolates Browser downloads (data/raw/*.tsv),
remove duplicates and unusable records, build clean labels, and write:

  data/processed/metadata_curated.tsv   one row per isolate, clean columns and labels
  data/processed/source_mapping.tsv     how free-text isolation sources were grouped
  data/interim/removed_isolates.tsv     every isolate removed, with the reason
  reports/curation_report.md            counts after every step (commit this one)

Run from the repository root:
    python workflow/scripts/01_curate_metadata.py

Rules (all written in the report as well):
  * Label per drug comes from the 'AST phenotypes' column: S -> 0, R -> 1.
    I, NS, ND or missing -> no label for that drug (kept in the table, excluded from that drug's model).
  * Duplicate BioSample: keep the record with an assembly and the highest N50.
  * Isolates without an assembly accession are removed (the local pipeline needs a genome file).
  * Sanity limits on assembly quality (see constants below).
  * Isolates with no SNP cluster get their own single-isolate group, so they are never dropped
    and never share a group with other isolates.
"""
import argparse
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

CORE_DRUGS = ["tetracycline", "ampicillin", "ceftriaxone"]
EXPLORATORY_DRUGS = ["ciprofloxacin"]
ALL_DRUGS = CORE_DRUGS + EXPLORATORY_DRUGS

# Sanity limits for Salmonella assemblies (adjustable; chosen from typical genome size, not from the data)
MIN_LENGTH = 4_300_000
MAX_LENGTH = 5_700_000
MAX_CONTIGS = 500
MIN_N50 = 20_000

# Primary analysis set: isolates whose NCBI AMR results were made with ONE reference-gene-catalog version
# (software 4.2.5 + catalog 2025-12-03.1). All other isolates are kept in the table for sensitivity analyses.
PRIMARY_CATALOG = "2025-12-03.1"


# ----------------------------------------------------------------------------- helpers
def parse_ast(value):
    """'amikacin=S,ampicillin=R,...' -> {'amikacin': 'S', 'ampicillin': 'R', ...}"""
    out = {}
    if not isinstance(value, str):
        return out
    for item in value.split(","):
        if "=" in item:
            drug, call = item.rsplit("=", 1)
            out[drug.strip()] = call.strip()
    return out


def to_label(call):
    """S -> 0, R -> 1, anything else (I, NS, ND, missing) -> NaN."""
    if call == "R":
        return 1.0
    if call == "S":
        return 0.0
    return np.nan


def source_category(text, ifsac, isolation_type):
    """Group free-text isolation source into a small number of documented categories."""
    if isolation_type == "clinical":
        return "human_clinical"
    t = text.lower().strip() if isinstance(text, str) else ""
    i = ifsac.lower().strip() if isinstance(ifsac, str) else ""
    if "turkey" in t:
        return "turkey"
    if re.search(r"chicken|broiler|\bhens?\b|gallus", t):
        return "chicken"
    if re.search(r"pork|swine|\bhogs?\b|\bsows?\b|\bpigs?\b|\bchops?\b", t):
        return "pork_swine"
    if re.search(r"beef|cattle|\bsteer\b|heifer|bovine|\bcows?\b", t):
        return "beef_cattle"
    if re.search(r"poultry|wing|thigh|\blegs?\b|breast|giblet|gizzard|heart|liver|drumstick", t):
        if "turkey" in i:
            return "turkey"
        if "chicken" in i:
            return "chicken"
        return "poultry_unspecified"
    if re.search(r"litter|caecum|cecum|feces|faeces|fecal|stool|intestin|cloaca|carcass", t):
        return "animal_or_farm_sample"
    if "turkey" in i:
        return "turkey"
    if "chicken" in i:
        return "chicken"
    if "pork" in i:
        return "pork_swine"
    if "beef" in i:
        return "beef_cattle"
    if "poultry" in i:
        return "poultry_unspecified"
    if "other animal" in i:
        return "other_animal"
    if t in ("", "missing"):
        return "unknown"
    return "other_or_unclear"


def serotype_from_row(computed_types, serovar):
    """Prefer NCBI's computed serotype; fall back to the submitted serovar."""
    if isinstance(computed_types, str):
        m = re.search(r"serotype=([^,]+)", computed_types)
        if m:
            return m.group(1).strip(), "computed"
    if isinstance(serovar, str) and serovar.strip():
        return serovar.strip(), "submitted"
    return "unknown", "none"


def read_raw(raw_dir):
    frames = []
    for f in sorted(Path(raw_dir).glob("*.tsv")):
        d = pd.read_csv(f, sep="\t", dtype=str)
        d.columns = [c.lstrip("#") for c in d.columns]
        d["source_file"] = f.name
        frames.append(d)
    if not frames:
        raise SystemExit(f"No .tsv files found in {raw_dir}")
    if len({tuple(d.columns) for d in frames}) != 1:
        raise SystemExit("Raw files do not have identical columns. Re-download with the same columns ticked.")
    return pd.concat(frames, ignore_index=True)


def check_file_labels(raw):
    """Each file name (drug_R.tsv / drug_S.tsv) must agree with its AST column. Returns list of problems."""
    problems = []
    for fname, g in raw.groupby("source_file"):
        stem = fname.replace(".tsv", "")
        if "_" not in stem:
            continue
        drug, expected = stem.rsplit("_", 1)
        got = g["AST phenotypes"].apply(lambda s: parse_ast(s).get(drug))
        bad = int((got != expected).sum())
        if bad:
            problems.append(f"{fname}: {bad} of {len(g)} rows do not have {drug}={expected} in the AST column")
    return problems


# ----------------------------------------------------------------------------- main pipeline
def curate(raw):
    log = []          # (step, n_remaining)
    removed = []      # DataFrames with a 'reason' column

    n_rows = len(raw)
    log.append(("Rows read from all raw files", n_rows))

    # 1. one row per isolate (the same isolate appears in several downloads)
    merged = raw.drop_duplicates(subset="Isolate").copy()
    files_per_isolate = raw.groupby("Isolate")["source_file"].apply(lambda s: ";".join(sorted(set(s))))
    merged["found_in_files"] = merged["Isolate"].map(files_per_isolate)
    log.append(("Unique isolates (same isolate in several files counted once)", len(merged)))

    # numeric columns
    for c in ["N50", "Length", "Contigs"]:
        merged[c] = pd.to_numeric(merged[c], errors="coerce")

    # 2. duplicate BioSample: keep record with assembly and highest N50
    merged["_has_asm"] = merged["Assembly"].notna().astype(int)
    merged = merged.sort_values(["BioSample", "_has_asm", "N50"], ascending=[True, False, False])
    dup = merged.duplicated(subset="BioSample", keep="first")
    if dup.any():
        r = merged[dup][["Isolate", "BioSample", "Assembly"]].copy()
        r["reason"] = "duplicate_biosample (another record of the same BioSample was kept)"
        removed.append(r)
    merged = merged[~dup].copy()
    log.append(("After removing duplicate BioSamples", len(merged)))

    # 3. same assembly accession under two isolates (should not happen; checked, not assumed)
    dup_asm = merged["Assembly"].notna() & merged.duplicated(subset="Assembly", keep="first")
    if dup_asm.any():
        r = merged[dup_asm][["Isolate", "BioSample", "Assembly"]].copy()
        r["reason"] = "duplicate_assembly_accession"
        removed.append(r)
    merged = merged[~dup_asm].copy()
    log.append(("After removing duplicate assembly accessions", len(merged)))

    # 4. no assembly accession
    no_asm = merged["Assembly"].isna()
    if no_asm.any():
        r = merged[no_asm][["Isolate", "BioSample", "Assembly"]].copy()
        r["reason"] = "no_assembly_accession"
        removed.append(r)
    merged = merged[~no_asm].copy()
    log.append(("After removing isolates without an assembly accession", len(merged)))

    # 5. assembly quality sanity limits
    bad_qc = (
        merged["Length"].isna() | merged["Contigs"].isna() | merged["N50"].isna()
        | (merged["Length"] < MIN_LENGTH) | (merged["Length"] > MAX_LENGTH)
        | (merged["Contigs"] > MAX_CONTIGS) | (merged["N50"] < MIN_N50)
    )
    if bad_qc.any():
        r = merged[bad_qc][["Isolate", "BioSample", "Assembly"]].copy()
        r["reason"] = "assembly_quality_limits"
        removed.append(r)
    merged = merged[~bad_qc].copy()
    log.append(("After assembly quality sanity limits", len(merged)))

    # ---- clean columns
    ast = merged["AST phenotypes"].apply(parse_ast)
    out = pd.DataFrame(index=merged.index)
    out["isolate"] = merged["Isolate"]
    out["biosample"] = merged["BioSample"]
    out["assembly"] = merged["Assembly"]
    st = [serotype_from_row(c, s) for c, s in zip(merged["Computed types"], merged["Serovar"])]
    out["serotype"] = [x[0] for x in st]
    out["serotype_source"] = [x[1] for x in st]
    out["country"] = merged["Location"].fillna("unknown").str.split(":").str[0].str.strip()
    year = pd.to_numeric(merged["Collection date"].str[:4], errors="coerce")
    out["collection_year"] = year.astype("Int64")
    out["bioproject"] = merged["BioProject"]
    out["isolation_type"] = merged["Isolation type"]
    out["host_group"] = np.where(merged["Isolation type"] == "clinical", "clinical", "non_clinical")
    out["isolation_source_raw"] = merged["Isolation source"]
    out["source_category"] = [
        source_category(t, i, it)
        for t, i, it in zip(merged["Isolation source"], merged["IFSAC category"], merged["Isolation type"])
    ]
    out["snp_cluster"] = merged["SNP cluster"]
    out["cluster_group"] = merged["SNP cluster"].where(merged["SNP cluster"].notna(), "SINGLETON_" + merged["Isolate"])
    out["n50"] = merged["N50"]
    out["genome_length"] = merged["Length"]
    out["contigs"] = merged["Contigs"]
    out["assembler"] = merged["Method"]
    out["amrfinder_version"] = merged["AMRFinderPlus version"]
    out["amrfinder_catalog_version"] = merged["PD Ref Gene Catalog version"]
    out["in_primary_set"] = (merged["PD Ref Gene Catalog version"] == PRIMARY_CATALOG)
    for d in ALL_DRUGS:
        calls = ast.apply(lambda x, d=d: x.get(d))
        out[f"ast_{d}"] = calls
        out[f"y_{d}"] = calls.apply(to_label)
    out["amr_genes_core"] = merged["AMR genotypes core"]
    out["ast_phenotypes_raw"] = merged["AST phenotypes"]
    out["found_in_files"] = merged["found_in_files"]
    out = out.sort_values("isolate").reset_index(drop=True)

    removed_df = pd.concat(removed, ignore_index=True) if removed else pd.DataFrame(columns=["Isolate", "BioSample", "Assembly", "reason"])
    return out, removed_df, log


def source_mapping_table(raw_unique, curated):
    m = raw_unique[["Isolate", "Isolation source", "IFSAC category", "Isolation type"]].merge(
        curated[["isolate", "source_category"]], left_on="Isolate", right_on="isolate", how="inner"
    )
    m["isolation_source_clean"] = m["Isolation source"].fillna("MISSING").str.strip().str.lower()
    t = (m.groupby(["isolation_source_clean", "Isolation type", "IFSAC category", "source_category"], dropna=False)
           .size().reset_index(name="n_isolates").sort_values("n_isolates", ascending=False))
    return t


# ----------------------------------------------------------------------------- report
def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def label_table(df):
    """Resistant / susceptible / excluded counts per antibiotic and host group."""
    rows = []
    for d in ALL_DRUGS:
        for grp in ["all", "non_clinical", "clinical"]:
            g = df if grp == "all" else df[df["host_group"] == grp]
            call = g[f"ast_{d}"]
            r = int((call == "R").sum()); s = int((call == "S").sum())
            excl = len(g) - r - s
            ncl = g.loc[g[f"y_{d}"] == 1, "cluster_group"].nunique()
            rows.append([d + (" (exploratory)" if d in EXPLORATORY_DRUGS else ""), grp, r, s, excl,
                         f"{100 * r / (r + s):.1f}%" if (r + s) else "NA", ncl])
    return pd.DataFrame(rows, columns=["Antibiotic", "Group", "R (label 1)", "S (label 0)",
                                       "Excluded (I/NS/ND/missing)", "R share", "Clusters containing R"])


def write_report(path, log, removed, curated, mapping, problems, raw_files):
    L = []
    L.append("# Metadata curation report\n")
    L.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}  \nScript: `workflow/scripts/01_curate_metadata.py`  \nRaw files: {', '.join(raw_files)}\n")
    L.append("## 1. Records kept after each step\n")
    L.append(md_table(pd.DataFrame(log, columns=["Step", "Isolates remaining"])) + "\n")
    L.append("Removed isolates by reason:\n")
    if len(removed):
        L.append(md_table(removed["reason"].value_counts().rename_axis("reason").reset_index(name="n")) + "\n")
    else:
        L.append("None.\n")
    L.append("Quality limits applied: genome length %s-%s bp, at most %s contigs, N50 at least %s bp. They are sanity limits chosen from typical *Salmonella* genome size.\n" % (f"{MIN_LENGTH:,}", f"{MAX_LENGTH:,}", MAX_CONTIGS, f"{MIN_N50:,}"))
    L.append("## 2. File-versus-label check\n")
    L.append("Each raw file should contain only the class in its name according to its own AST column.\n")
    L.append(("Problems found:\n- " + "\n- ".join(problems) + "\n") if problems else "No problems found: every row agrees with its file name.\n")
    L.append("## 3. Labels per antibiotic (after curation, all isolates)\n")
    L.append(md_table(label_table(curated)) + "\n")
    primary = curated[curated["in_primary_set"]]
    L.append(f"## 3b. Labels per antibiotic, primary dataset only (catalog {PRIMARY_CATALOG}; {len(primary)} of {len(curated)} isolates)\n")
    L.append("The primary dataset uses isolates whose NCBI AMR results come from one software and one reference-gene-catalog version. "
             "The other isolates are kept for sensitivity analyses.\n")
    L.append(md_table(label_table(primary)) + "\n")
    L.append("## 4. Independent groups (SNP clusters)\n")
    vc = curated["cluster_group"].value_counts()
    singles = int(curated["cluster_group"].str.startswith("SINGLETON_").sum())
    L.append(f"- Groups: {len(vc)} (including {singles} isolates without a SNP cluster, each its own group)\n- Largest group: {vc.iloc[0]} isolates ({100 * vc.iloc[0] / len(curated):.1f}%); top 5 groups: {100 * vc.iloc[:5].sum() / len(curated):.1f}%\n")
    L.append("## 5. Source, place, time, study\n")
    L.append("Source categories:\n\n" + md_table(curated["source_category"].value_counts().rename_axis("source_category").reset_index(name="n")) + "\n")
    L.append("Countries (top 8):\n\n" + md_table(curated["country"].value_counts().head(8).rename_axis("country").reset_index(name="n")) + "\n")
    L.append("Collection years (grouped):\n\n" + md_table(curated["collection_year"].dropna().astype(int).apply(lambda y: "2002-2009" if y <= 2009 else ("2010-2014" if y <= 2014 else ("2015-2018" if y <= 2018 else ("2019-2021" if y <= 2021 else "2022+")))).value_counts().rename_axis("years").reset_index(name="n").sort_values("years")) + f"\n\nMissing collection year: {int(curated['collection_year'].isna().sum())}\n")
    L.append("BioProjects (top 5 of %d):\n\n" % curated["bioproject"].nunique() + md_table(curated["bioproject"].value_counts().head(5).rename_axis("bioproject").reset_index(name="n")) + "\n")
    L.append("Assemblers (top 5):\n\n" + md_table(curated["assembler"].fillna("missing").value_counts().head(5).rename_axis("assembler").reset_index(name="n")) + "\n")
    L.append("AMRFinderPlus versions used by NCBI:\n\n" + md_table(curated["amrfinder_version"].fillna("missing").value_counts().rename_axis("version").reset_index(name="n")) + "\n")
    L.append("Reference gene catalog (database) versions used by NCBI:\n\n" + md_table(curated["amrfinder_catalog_version"].fillna("missing").value_counts().rename_axis("catalog_version").reset_index(name="n")) + "\n")
    L.append("## 6. Source text that could not be grouped confidently (top 15)\n")
    unclear = mapping[mapping["source_category"].isin(["other_or_unclear", "poultry_unspecified", "unknown"])].head(15)
    L.append(md_table(unclear[["isolation_source_clean", "source_category", "n_isolates"]]) + "\n" if len(unclear) else "None.\n")
    L.append("## 7. Known limitations\n")
    L.append("- Antibiotic labels are submitter-supplied and not vetted by NCBI beyond basic checks; the test method is not in the table.\n- Isolates are heavily clustered, so isolate counts overstate independent evidence; use `cluster_group` for splitting.\n- Most isolates are from the USA.\n- The 'clinical' group is small (human isolates only).\n- `source_category` comes from text rules (`data/processed/source_mapping.tsv`) and should be reviewed.\n")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(L), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--processed-dir", default="data/processed")
    ap.add_argument("--interim-dir", default="data/interim")
    ap.add_argument("--report", default="reports/curation_report.md")
    a = ap.parse_args()

    raw = read_raw(a.raw_dir)
    problems = check_file_labels(raw)
    curated, removed, log = curate(raw)
    uniq = raw.drop_duplicates(subset="Isolate")
    mapping = source_mapping_table(uniq, curated)

    Path(a.processed_dir).mkdir(parents=True, exist_ok=True)
    Path(a.interim_dir).mkdir(parents=True, exist_ok=True)
    curated.to_csv(Path(a.processed_dir) / "metadata_curated.tsv", sep="\t", index=False)
    mapping.to_csv(Path(a.processed_dir) / "source_mapping.tsv", sep="\t", index=False)
    removed.to_csv(Path(a.interim_dir) / "removed_isolates.tsv", sep="\t", index=False)
    write_report(a.report, log, removed, curated, mapping, problems, sorted(raw["source_file"].unique()))
    print(f"Curated isolates: {len(curated)}  |  removed: {len(removed)}  |  file/label problems: {len(problems)}")
    print(f"Report: {a.report}")


if __name__ == "__main__":
    main()
