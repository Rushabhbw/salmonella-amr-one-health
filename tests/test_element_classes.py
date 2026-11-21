"""Tests for workflow/scripts/05_map_elements_to_classes.py (run from the repository root: pytest)."""
import importlib.util
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "05_map_elements_to_classes.py"
spec = importlib.util.spec_from_file_location("map_classes", SCRIPT)
mc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mc)

HEADER = "#allele\tgene_family\twhitelisted_taxa\tproduct_name\tscope\ttype\tsubtype\tclass\tsubclass\tdb_version\n"
ROWS = [
    "blaTEM-1\tblaTEM\t\tTEM-1 beta-lactamase\tcore\tAMR\tAMR\tBETA-LACTAM\tBETA-LACTAM\t2026-10-29.1",
    "blaCMY-2\tblaCMY\t\tCMY-2 beta-lactamase\tcore\tAMR\tAMR\tBETA-LACTAM\tCEPHALOSPORIN\t2026-10-29.1",
    "tet(A)\ttet(A)\t\ttetracycline efflux\tcore\tAMR\tAMR\tTETRACYCLINE\tTETRACYCLINE\t2026-10-29.1",
    "gyrA_D87Y\tgyrA\tSalmonella\tquinolone resistant GyrA\tcore\tAMR\tPOINT\tQUINOLONE\tQUINOLONE\t2026-10-29.1",
    "\tblaCTX-M\t\tCTX-M family beta-lactamase\tcore\tAMR\tAMR\tBETA-LACTAM\tCEPHALOSPORIN\t2026-10-29.1",
    "sul1\tsul1\t\tsulfonamide-resistant DHPS\tcore\tAMR\tAMR\tSULFONAMIDE\tSULFONAMIDE\t2026-10-29.1",
    "weird\tweird\t\tx\tcore\tAMR\tAMR\tBETA-LACTAM\tCEPHALOSPORIN\t2026-10-29.1",
    "weird\tweird\t\tx\tplus\tSTRESS\tBIOCIDE\tBIOCIDE\tBIOCIDE\t2026-10-29.1",
]


def make_catalog(tmp_path):
    p = tmp_path / "cat.txt"
    p.write_text(HEADER + "\n".join(ROWS) + "\n", encoding="utf-8")
    return p


def test_drug_rules():
    plain = lambda pairs: {k: v for k, v in mc.drugs_for(pairs).items() if "multidrug" not in k}
    assert plain([("BETA-LACTAM", "BETA-LACTAM")]) == {"tetracycline": 0, "ampicillin": 1, "ceftriaxone": 0, "ciprofloxacin": 0}
    assert mc.drugs_for([("BETA-LACTAM", "CEPHALOSPORIN")])["ceftriaxone"] == 1
    assert mc.drugs_for([("BETA-LACTAM", "CARBAPENEM")])["ceftriaxone"] == 0
    assert mc.drugs_for([("TETRACYCLINE", "TETRACYCLINE")])["tetracycline"] == 1
    assert mc.drugs_for([("QUINOLONE", "QUINOLONE")])["ciprofloxacin"] == 1
    assert sum(mc.drugs_for([("SULFONAMIDE", "SULFONAMIDE")]).values()) == 0


def test_mapping_allele_family_and_unmatched(tmp_path):
    cat = mc.read_catalog(make_catalog(tmp_path))
    d = pd.DataFrame({"element": ["blaTEM-1", "blaCMY-2", "tet(A)", "gyrA_D87Y", "blaCTX-M", "sul1", "unknownX", "weird"]})
    m = mc.map_elements(d, cat).set_index("element")
    assert m.loc["blaTEM-1", ["ampicillin", "ceftriaxone", "tetracycline"]].tolist() == [1, 0, 0]
    assert m.loc["blaCMY-2", ["ampicillin", "ceftriaxone"]].tolist() == [1, 1]
    assert m.loc["tet(A)", "tetracycline"] == 1 and m.loc["gyrA_D87Y", "ciprofloxacin"] == 1
    assert m.loc["blaCTX-M", "matched_on"] == "gene_family" and m.loc["blaCTX-M", "ceftriaxone"] == 1
    assert m.loc["sul1", ["tetracycline", "ampicillin", "ceftriaxone", "ciprofloxacin"]].sum() == 0
    assert m.loc["unknownX", "matched_on"] == "none"
    assert m.loc["weird", "scope"] == "core"      # the core row is preferred over the plus row


def test_missing_column_stops_with_clear_message(tmp_path):
    p = tmp_path / "bad.txt"; p.write_text("#allele\tfoo\nx\ty\n", encoding="utf-8")
    try:
        mc.read_catalog(p)
        assert False, "should have stopped"
    except SystemExit as e:
        assert "missing expected columns" in str(e)


def test_end_to_end(tmp_path):
    cat = make_catalog(tmp_path)
    dic = tmp_path / "dic.tsv"
    pd.DataFrame({"element": ["blaTEM-1", "tet(A)", "unknownX"], "kind": ["gene"] * 3,
                  "n_isolates_primary": [5, 4, 1], "n_groups_primary": [3, 2, 1]}).to_csv(dic, sep="\t", index=False)
    out, rep = tmp_path / "ec.tsv", tmp_path / "rep.md"
    code = mc.main(["--catalog", str(cat), "--dictionary", str(dic), "--uncertain", str(tmp_path / "none.tsv"),
                    "--out", str(out), "--report", str(rep)])
    assert code == 0 and len(pd.read_csv(out, sep="\t")) == 3
    assert "Not found in the catalog: 1" in rep.read_text(encoding="utf-8")
    assert mc.main(["--catalog", str(tmp_path / "nope.txt"), "--dictionary", str(dic), "--uncertain", str(tmp_path / "none.tsv"),
                    "--out", str(out), "--report", str(rep)]) == 1


def test_family_symbol_with_mixed_alleles_is_not_cephalosporin_evidence():
    # like blaTEM: one entry is TEM-1-like (no cephalosporin activity), another is a cephalosporin-active allele
    pairs = [("BETA-LACTAM", "BETA-LACTAM"), ("BETA-LACTAM", "CEPHALOSPORIN")]
    r = mc.drugs_for(pairs)
    assert r["ampicillin"] == 1 and r["ceftriaxone"] == 0
    assert mc.drugs_for([("BETA-LACTAM", "CEPHALOSPORIN"), ("BETA-LACTAM", "CEPHALOSPORIN/TANIBORBACTAM")])["ceftriaxone"] == 1
    assert mc.drugs_for([])["tetracycline"] == 0                       # no entries: never evidence


def test_multidrug_class_is_sensitivity_only():
    r = mc.drugs_for([("MULTIDRUG", "AMPICILLIN/CHLORAMPHENICOL/FLUOROQUINOLONE/RIFAMPIN/TETRACYCLINE")])
    assert r["tetracycline"] == 0 and r["tetracycline_with_multidrug"] == 1
    assert r["ampicillin"] == 0 and r["ampicillin_with_multidrug"] == 1


def test_weak_only_elements_are_added(tmp_path):
    cat = mc.read_catalog(make_catalog(tmp_path))
    d = pd.DataFrame({"element": ["tet(A)"], "kind": ["gene"]})
    m = mc.map_elements(d, cat, extra_elements=["tet(A)", "blaTEM-1"]).set_index("element")
    assert list(m.index) == ["tet(A)", "blaTEM-1"] and m.loc["blaTEM-1", "kind"] == "weak_only" and m.loc["blaTEM-1", "ampicillin"] == 1
