"""Tests for workflow/scripts/02_select_subset.py (run from the repository root: pytest)."""
import importlib.util
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "02_select_subset.py"
spec = importlib.util.spec_from_file_location("select_subset", SCRIPT)
ss = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ss)


def fake_meta(n_clin=150, n_non=600):
    rows = []
    for i in range(n_clin):
        rows.append({"isolate": f"C{i}", "assembly": f"GCA_C{i}.1", "host_group": "clinical", "cluster_group": f"PDS_C{i // 2}"})
    for i in range(n_non):
        rows.append({"isolate": f"N{i}", "assembly": f"GCA_N{i}.1", "host_group": "non_clinical", "cluster_group": f"PDS_N{i // 20}"})
    rows.append({"isolate": "X0", "assembly": None, "host_group": "clinical", "cluster_group": "PDS_X"})
    return pd.DataFrame(rows)


def test_sizes_cap_and_pilot():
    sub = ss.select_subset(fake_meta())
    assert (sub["host_group"] == "clinical").sum() == 100
    # 600 non-clinical isolates in 30 groups, cap of 3 per group -> only 90 available, so fewer than 200 are returned
    assert (sub["host_group"] == "non_clinical").sum() == 90
    assert sub["cluster_group"].value_counts().max() <= ss.MAX_PER_GROUP
    pilot = sub[sub["in_pilot"]]
    assert len(pilot) == 20 and (pilot["host_group"] == "clinical").sum() == 6
    assert sub["assembly"].notna().all() and sub["isolate"].is_unique


def test_reproducible_with_same_seed():
    a = ss.select_subset(fake_meta())["isolate"].tolist()
    b = ss.select_subset(fake_meta())["isolate"].tolist()
    assert a == b


def test_accession_file_has_linux_line_endings(tmp_path):
    p = tmp_path / "acc.txt"
    ss.write_lines(p, ["GCA_1.1", "GCA_2.1"])
    assert p.read_bytes() == b"GCA_1.1\nGCA_2.1\n"
