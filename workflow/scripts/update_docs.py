#!/usr/bin/env python3
"""
update_docs.py

One command that refreshes every auto-generated document:
  1. docs/data_provenance.md     (from data/raw and config/download_manifest.tsv)
  2. reports/curation_report.md  (plus the curated tables in data/processed)

Run from the repository root:
    python workflow/scripts/update_docs.py

Then look at `git status` / `git diff`, read what changed, and commit.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = ["00_log_data_provenance.py", "01_curate_metadata.py"]


def main():
    worst = 0
    for step in STEPS:
        print(f"\n=== {step} ===")
        code = subprocess.call([sys.executable, str(HERE / step)])
        worst = max(worst, code)
        if code != 0:
            print(f"{step} reported problems (exit code {code}). Read the messages above before committing.")
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
