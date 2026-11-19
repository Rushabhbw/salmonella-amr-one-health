"""Tests for workflow/scripts/01_curate_metadata.py (run from the repository root: pytest)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "01_curate_metadata.py"
spec = importlib.util.spec_from_file_location("curate_metadata", SCRIPT)
cm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cm)


def row(isolate, biosample, assembly="auto", n50=100000, length=4_900_000, contigs=50,
        ast="ampicillin=R,tetracycline=S,ceftriaxone=I,ciprofloxacin=NS", snp="PDS1.1",
        itype="environmental/other", source="chicken breast", ifsac=None, src_file="ampicillin_R.tsv"):
    if assembly == "auto":
        assembly = f"GCA_{isolate}.1"  # every fake isolate gets its own assembly accession
    return {
        "Isolate": isolate, "BioSample": biosample, "Assembly": assembly, "N50": str(n50), "Length": str(length),
        "Contigs": str(contigs), "AST phenotypes": ast, "Computed types": "antigen_formula=4:i:-,serotype=Typhimurium",
        "Serovar": None, "Location": "USA:NY", "Collection date": "2015-03-02", "BioProject": "PRJNA1",
        "Isolation type": itype, "Isolation source": source, "IFSAC category": ifsac, "SNP cluster": snp,
        "Method": "SKESA v. 2.2", "AMRFinderPlus version": "4.2.5", "PD Ref Gene Catalog version": "2025-12-03.1", "AMR genotypes core": "blaCMY-2=COMPLETE",
        "source_file": src_file,
    }


def test_parse_ast_and_labels():
    d = cm.parse_ast("ampicillin=R,tetracycline=S,ceftriaxone=I,trimethoprim-sulfamethoxazole=NS")
    assert d["ampicillin"] == "R" and d["trimethoprim-sulfamethoxazole"] == "NS"
    assert cm.to_label("R") == 1.0 and cm.to_label("S") == 0.0
    assert np.isnan(cm.to_label("I")) and np.isnan(cm.to_label("NS")) and np.isnan(cm.to_label(None))
    assert cm.parse_ast(float("nan")) == {}


def test_source_category():
    f = cm.source_category
    assert f("stool", None, "clinical") == "human_clinical"
    assert f("ground turkey", None, "environmental/other") == "turkey"
    assert f("Chicken Breasts", None, "environmental/other") == "chicken"
    assert f("animal-swine-sow", None, "environmental/other") == "pork_swine"
    assert f("ground beef", None, "environmental/other") == "beef_cattle"
    assert f("wings", None, "environmental/other") == "poultry_unspecified"
    assert f("wings", "chicken", "environmental/other") == "chicken"
    assert f(None, None, "environmental/other") == "unknown"
    assert f("shrimp (sh)", None, "environmental/other") == "other_or_unclear"


def test_serotype_fallback():
    assert cm.serotype_from_row("antigen_formula=4:i:-,serotype=Typhimurium", "x") == ("Typhimurium", "computed")
    assert cm.serotype_from_row(None, "Kentucky") == ("Kentucky", "submitted")
    assert cm.serotype_from_row(None, None) == ("unknown", "none")


def test_curate_removes_duplicates_and_missing_assembly_and_keeps_singletons():
    raw = pd.DataFrame([
        row("PDT1", "SAMN1"),
        row("PDT1", "SAMN1", src_file="ceftriaxone_R.tsv"),           # same isolate in a second file
        row("PDT2", "SAMN2", assembly=None),                             # no assembly -> removed
        row("PDT3", "SAMN3"), row("PDT4", "SAMN3", n50=500000),         # same BioSample -> keep higher N50
        row("PDT5", "SAMN5", snp=None),                                  # no SNP cluster -> own group
        row("PDT6", "SAMN6", length=9_000_000),                          # absurd genome length -> removed
    ])
    out, removed, log = cm.curate(raw)
    assert sorted(out["isolate"]) == ["PDT1", "PDT4", "PDT5"]
    assert set(removed["reason"].str.split(" ").str[0]) == {"no_assembly_accession", "duplicate_biosample", "assembly_quality_limits"}
    assert out.loc[out.isolate == "PDT5", "cluster_group"].iloc[0] == "SINGLETON_PDT5"
    r = out.iloc[0]
    assert r["y_ampicillin"] == 1.0 and r["y_tetracycline"] == 0.0
    assert np.isnan(r["y_ceftriaxone"]) and np.isnan(r["y_ciprofloxacin"])
    assert out["host_group"].eq("non_clinical").all()


def test_file_label_check_detects_mismatch():
    raw = pd.DataFrame([row("PDT1", "SAMN1", ast="ampicillin=S", src_file="ampicillin_R.tsv")])
    assert len(cm.check_file_labels(raw)) == 1


def test_duplicate_assembly_accession_removed():
    raw = pd.DataFrame([row("PDT1", "SAMN1", assembly="GCA_X.1"), row("PDT2", "SAMN2", assembly="GCA_X.1")])
    out, removed, _ = cm.curate(raw)
    assert len(out) == 1
    assert removed["reason"].iloc[0] == "duplicate_assembly_accession"


def test_primary_set_flag():
    raw = pd.DataFrame([
        row("PDT1", "SAMN1"),
        {**row("PDT2", "SAMN2"), "PD Ref Gene Catalog version": "2026-08-07.1"},
    ])
    out, _, _ = cm.curate(raw)
    flags = dict(zip(out["isolate"], out["in_primary_set"]))
    assert flags == {"PDT1": True, "PDT2": False}
