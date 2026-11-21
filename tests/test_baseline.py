"""Tests for workflow/scripts/06_baseline_rules.py (run from the repository root: pytest)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "06_baseline_rules.py"
spec = importlib.util.spec_from_file_location("baseline", SCRIPT)
bl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bl)


def test_metrics_from_counts_known_values():
    m = bl.metrics_from_counts([8, 2, 85, 5])          # TP, FN, TN, FP
    assert round(float(m["sensitivity"]), 3) == 0.8 and round(float(m["specificity"]), 4) == round(85 / 90, 4)
    assert round(float(m["vme"]), 3) == 0.2 and round(float(m["me"]), 4) == round(5 / 90, 4)
    expected_mcc = (8 * 85 - 5 * 2) / np.sqrt((8 + 5) * (8 + 2) * (85 + 5) * (85 + 2))
    assert abs(float(m["mcc"]) - expected_mcc) < 1e-9
    z = bl.metrics_from_counts([0, 0, 10, 0])
    assert np.isnan(z["sensitivity"]) and z["specificity"] == 1.0


def make_data():
    # 10 isolates, 5 groups. tet(A) is evidence for tetracycline; blaTEM-1 for ampicillin.
    meta = pd.DataFrame({
        "isolate": [f"I{i}" for i in range(10)],
        "in_primary_set": [True] * 9 + [False],
        "host_group": ["non_clinical"] * 7 + ["clinical"] * 3,
        "cluster_group": ["G0", "G0", "G1", "G1", "G2", "G3", "G3", "G4", "G4", "G4"],
        "y_tetracycline": [1, 1, 1, 0, 0, 0, 1, 0, 1, 1],
        "y_ampicillin": [1, 0, 0, 0, 0, 0, 0, 1, 1, 0],
        "y_ceftriaxone": [np.nan] * 10,
        "y_ciprofloxacin": [0] * 10,
    })
    feats = pd.DataFrame({
        "isolate": meta["isolate"], "assembly": "x", "in_primary_set": meta["in_primary_set"],
        "host_group": meta["host_group"], "cluster_group": meta["cluster_group"],
        "tet(A)": [1, 1, 0, 1, 0, 0, 0, 0, 1, 1], "blaTEM-1": [1, 0, 0, 0, 0, 0, 0, 1, 1, 0],
    })
    mapping = pd.DataFrame({"element": ["tet(A)", "blaTEM-1"], "tetracycline": [1, 0], "ampicillin": [0, 1],
                            "ceftriaxone": [0, 0], "ciprofloxacin": [0, 0]})
    return meta, feats, mapping


def test_predictions_and_counts():
    meta, feats, mapping = make_data()
    preds, sections = bl.build(meta, feats, mapping)
    assert preds["pred_tetracycline"].tolist() == [1, 1, 0, 1, 0, 0, 0, 0, 1, 1]
    sec = {s["drug"]: s for s in sections}
    tab = sec["tetracycline"]["table"]
    row = tab[tab["Group"] == "all (primary dataset)"].iloc[0]
    # primary isolates I0..I8: y = 1,1,1,0,0,0,1,0,1 ; pred = 1,1,0,1,0,0,0,0,1
    assert row["TP/FN/TN/FP"] == "3/2/3/1"
    assert row["R / S isolates"] == "5 / 4"
    sec_amp = sec["ampicillin"]["table"]
    assert sec_amp[sec_amp["Group"] == "all (primary dataset)"].iloc[0]["TP/FN/TN/FP"] == "3/0/6/0"   # blaTEM-1 matches ampicillin labels exactly
    # ceftriaxone has no labels and no mapped elements: must not crash
    assert sec["ceftriaxone"]["elems"] == []


def test_bootstrap_is_reproducible_and_brackets_point_estimate():
    counts = np.array([[5, 1, 20, 2]] * 12 + [[1, 3, 8, 1]] * 6, dtype=float)
    a = bl.cluster_bootstrap(counts, n_boot=500, seed=1)
    b = bl.cluster_bootstrap(counts, n_boot=500, seed=1)
    assert a == b
    point = bl.metrics_from_counts(counts.sum(axis=0))
    lo, hi = a["sensitivity"]
    assert lo <= float(point["sensitivity"]) <= hi
    single = bl.cluster_bootstrap(np.array([[3, 1, 5, 1]], dtype=float))
    assert all(np.isnan(x) for pair in single.values() for x in pair)    # one group: no interval


def test_end_to_end(tmp_path):
    meta, feats, mapping = make_data()
    meta.to_csv(tmp_path / "m.tsv", sep="\t", index=False); feats.to_csv(tmp_path / "f.tsv", sep="\t", index=False)
    mapping.to_csv(tmp_path / "map.tsv", sep="\t", index=False)
    code = bl.main(["--metadata", str(tmp_path / "m.tsv"), "--features", str(tmp_path / "f.tsv"), "--mapping", str(tmp_path / "map.tsv"),
                    "--out", str(tmp_path / "p.tsv"), "--report", str(tmp_path / "r.md")])
    assert code == 0 and "## tetracycline" in (tmp_path / "r.md").read_text(encoding="utf-8")
    assert bl.main(["--mapping", str(tmp_path / "missing.tsv")]) == 1


def test_sensitivity_rows_partial_calls_and_bioproject_table():
    meta, feats, mapping = make_data()
    meta["bioproject"] = ["PRJ1"] * 4 + ["PRJ2"] * 6
    mapping["tetracycline_with_multidrug"] = [1, 0]; mapping["ampicillin_with_multidrug"] = [0, 1]
    mapping["ceftriaxone_with_multidrug"] = [0, 0]; mapping["ciprofloxacin_with_multidrug"] = [0, 0]
    # isolate I2 is resistant (y=1) but has no strong tet(A); give it a partial tet(A) call
    unc = pd.DataFrame({"isolate": meta["isolate"], "assembly": "x", "in_primary_set": meta["in_primary_set"],
                        "host_group": meta["host_group"], "cluster_group": meta["cluster_group"], "tet(A)": [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]})
    preds, sections = bl.build(meta, feats, mapping, unc)
    sec = {s["drug"]: s for s in sections}["tetracycline"]
    groups = sec["table"]["Group"].tolist()
    assert "primary, also counting partial calls" in groups and "primary, also counting MULTIDRUG-class elements" in groups
    assert sec["fn_weak"] == 1                                     # the partial call explains one missed resistant isolate
    part = sec["table"][sec["table"]["Group"] == "primary, also counting partial calls"].iloc[0]
    assert part["TP/FN/TN/FP"] == "4/1/3/1"                        # I2 moves from FN to TP
    assert sec["bioproject"] is not None and set(sec["bioproject"].columns) >= {"bioproject", "errors", "error_rate"}
