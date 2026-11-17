"""Tests for workflow/scripts/03_compare_amrfinder.py (run from the repository root: pytest)."""
import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cmp_amr", ROOT / "workflow" / "scripts" / "03_compare_amrfinder.py")
ca = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ca)

REAL = ROOT / "tests" / "data" / "GCA_006497535.1.amrfinder.tsv"   # real output from AMRFinderPlus 4.2.7
NCBI = "aph(3'')-Ib=COMPLETE,aph(6)-Id=COMPLETE,blaTEM-1=COMPLETE,sul2=COMPLETE,tet(B)=COMPLETE"


def test_parse_ncbi_core():
    d = ca.parse_ncbi_core("blaTEM-1=COMPLETE,gyrA_S83F=POINT,aac(6')-Iaa=PARTIAL_CONTIG_END")
    assert d == {"blaTEM-1": "COMPLETE", "gyrA_S83F": "POINT", "aac(6')-Iaa": "PARTIAL_CONTIG_END"}
    assert ca.parse_ncbi_core(float("nan")) == {} and ca.parse_ncbi_core(None) == {}


def test_real_amrfinder_output_matches_ncbi_string():
    local = ca.read_local(REAL)
    assert set(local) == {"tet(B)", "aph(6)-Id", "aph(3'')-Ib", "sul2", "blaTEM-1"}   # 10 table rows, 5 distinct elements
    c = ca.compare_sets(local, ca.parse_ncbi_core(NCBI))
    assert c["identical"] and c["jaccard"] == 1.0 and c["n_shared"] == 5


def test_differences_are_detected():
    local = ca.read_local(REAL)
    ncbi = ca.parse_ncbi_core("aph(3'')-Ib=COMPLETE,aph(6)-Id=COMPLETE,blaTEM-1=COMPLETE,tet(B)=COMPLETE,blaCMY-2=COMPLETE")
    c = ca.compare_sets(local, ncbi)
    assert not c["identical"]
    assert c["only_local"] == ["sul2"] and c["only_ncbi"] == ["blaCMY-2"]
    assert round(c["jaccard"], 3) == round(4 / 6, 3)


def test_empty_result_file_means_no_elements(tmp_path):
    p = tmp_path / "x.amrfinder.tsv"
    p.write_text("Protein id\tContig id\tElement symbol\tScope\tType\tMethod\n", encoding="utf-8")
    assert ca.read_local(p) == {}
    assert ca.compare_sets({}, {})["identical"]


def test_end_to_end(tmp_path):
    amr = tmp_path / "amr"; amr.mkdir()
    (amr / "GCA_006497535.1.amrfinder.tsv").write_bytes(REAL.read_bytes())
    (amr / "versions.txt").write_text("Database directory: /home/someone/secret_folder/2026-10-29.1\nDatabase version: 2026-10-29.1\n", encoding="utf-8")
    meta = pd.DataFrame([{"assembly": "GCA_006497535.1", "isolate": "PDT000161115.2", "host_group": "clinical",
                          "serotype": "I 4", "amrfinder_version": "4.2.5", "amrfinder_catalog_version": "2025-12-03.1",
                          "amr_genes_core": NCBI}])
    mp = tmp_path / "meta.tsv"; meta.to_csv(mp, sep="\t", index=False)
    rep, tab = tmp_path / "rep.md", tmp_path / "tab.tsv"
    code = ca.main(["--amr-dir", str(amr), "--metadata", str(mp), "--report", str(rep), "--table", str(tab)])
    assert code == 0
    text = rep.read_text(encoding="utf-8")
    assert "Identical sets of resistance elements: 1 (100.0%)" in text and "Database version: 2026-10-29.1" in text
    assert "secret_folder" not in text
    assert pd.read_csv(tab, sep="\t")["identical"].all()
