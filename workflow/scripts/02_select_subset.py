#!/usr/bin/env python3
"""
02_select_subset.py

Pick a small, fair subset of isolates for re-running AMRFinderPlus ourselves
(to check NCBI's pre-computed AMR results and to test the pipeline end to end).

Reads : data/processed/metadata_curated.tsv   (made by 01_curate_metadata.py)
Writes: data/processed/subset_300.tsv                  isolates chosen, with labels and metadata
        data/processed/subset_300_accessions.txt       assembly accessions, one per line (Linux line endings)
        data/processed/subset_pilot_accessions.txt     first 20 of them, to test the workflow quickly

How the choice is made (all rules are written in the printed summary):
  * fixed random seed, so everyone gets the same subset
  * 100 human clinical + 200 non-clinical isolates
  * at most MAX_PER_GROUP isolates per SNP-cluster group within each host group (so up to 2 x MAX_PER_GROUP overall),
    so one big family cannot fill the subset
  * the pilot has 6 clinical + 14 non-clinical isolates

Run from the repository root:
    python workflow/scripts/02_select_subset.py
"""
import argparse
from pathlib import Path

import pandas as pd

SEED = 42
N_CLINICAL = 100
N_NON_CLINICAL = 200
MAX_PER_GROUP = 3
PILOT_CLINICAL = 6
PILOT_NON_CLINICAL = 14


def draw(df, n, seed, max_per_group):
    """Random draw of up to n rows with at most max_per_group rows per cluster_group."""
    shuffled = df.sample(frac=1.0, random_state=seed)
    capped = shuffled.groupby("cluster_group", sort=False).head(max_per_group)
    return capped.head(n)


def select_subset(meta, seed=SEED, n_clin=N_CLINICAL, n_non=N_NON_CLINICAL, max_per_group=MAX_PER_GROUP,
                  pilot_clin=PILOT_CLINICAL, pilot_non=PILOT_NON_CLINICAL):
    meta = meta[meta["assembly"].notna()]
    clin = draw(meta[meta["host_group"] == "clinical"], n_clin, seed, max_per_group)
    non = draw(meta[meta["host_group"] == "non_clinical"], n_non, seed, max_per_group)
    pilot = pd.concat([clin.head(pilot_clin), non.head(pilot_non)])
    rest = pd.concat([clin.iloc[pilot_clin:], non.iloc[pilot_non:]]).sort_values("isolate")
    # pilot first (sorted), then the rest (sorted)
    ordered = pd.concat([pilot.sort_values("isolate"), rest]).reset_index(drop=True)
    ordered["in_pilot"] = [True] * len(pilot) + [False] * len(rest)
    return ordered


def write_lines(path, items):
    # newline="\n" keeps Linux line endings even when this script runs on Windows
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(items) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default="data/processed/metadata_curated.tsv")
    ap.add_argument("--out-dir", default="data/processed")
    a = ap.parse_args(argv)

    meta = pd.read_csv(a.metadata, sep="\t", dtype=str)
    sub = select_subset(meta)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sub.to_csv(out / "subset_300.tsv", sep="\t", index=False, lineterminator="\n")
    write_lines(out / "subset_300_accessions.txt", sub["assembly"].tolist())
    write_lines(out / "subset_pilot_accessions.txt", sub.loc[sub["in_pilot"], "assembly"].tolist())

    print(f"Subset: {len(sub)} isolates (pilot {int(sub['in_pilot'].sum())})  seed={SEED}  max per cluster group={MAX_PER_GROUP}")
    print("Host groups:", sub["host_group"].value_counts().to_dict())
    print("Largest cluster group in subset:", int(sub["cluster_group"].value_counts().iloc[0]), "isolates (cap is applied separately in each host group)")
    for d in ["tetracycline", "ampicillin", "ceftriaxone", "ciprofloxacin"]:
        c = sub[f"ast_{d}"].value_counts(dropna=False).to_dict()
        print(f"  {d}: {c}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
