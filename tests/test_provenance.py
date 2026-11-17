"""Tests for workflow/scripts/00_log_data_provenance.py (run from the repository root: pytest)."""
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "00_log_data_provenance.py"
spec = importlib.util.spec_from_file_location("log_prov", SCRIPT)
lp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lp)


def make_tsv(path, n, cols=("#Organism group", "Isolate")):
    lines = ["\t".join(cols)] + ["\t".join(["Salmonella enterica", f"PDT{i}"]) for i in range(n)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_manifest(path, rows):
    lines = ["file\tquery\texpected_rows\texpected_date\tdownload_date\tnotes"]
    for f, exp in rows:
        lines.append(f"{f}\tquery for {f}\t{exp}\t2026-11-16\t\t")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(tmp_path, files, manifest_rows, **kw):
    raw = tmp_path / "raw"; raw.mkdir()
    for name, n in files.items():
        make_tsv(raw / name, n)
    man = tmp_path / "manifest.tsv"; make_manifest(man, manifest_rows)
    out = tmp_path / "log.md"
    code = lp.main(["--raw-dir", str(raw), "--manifest", str(man), "--out", str(out)])
    return code, out.read_text(encoding="utf-8")


def test_exact_counts_ok(tmp_path):
    code, text = run(tmp_path, {"a_R.tsv": 5, "a_S.tsv": 7}, [("a_R.tsv", 5), ("a_S.tsv", 7)])
    assert code == 0 and "OK (exact)" in text and "Problems found\n\nNone." in text
    assert "Organism group" in text and "#Organism" not in text


def test_small_drift_accepted(tmp_path):
    code, text = run(tmp_path, {"a_R.tsv": 1000}, [("a_R.tsv", 1005)])
    assert code == 0 and "within 1.0%" in text


def test_big_drift_flagged(tmp_path):
    code, text = run(tmp_path, {"a_R.tsv": 100}, [("a_R.tsv", 200)])
    assert code == 1 and "CHECK" in text


def test_missing_and_unlisted_files_flagged(tmp_path):
    code, text = run(tmp_path, {"extra.tsv": 3}, [("missing.tsv", 3)])
    assert code == 1
    assert "missing from" in text and "not in the manifest" in text
