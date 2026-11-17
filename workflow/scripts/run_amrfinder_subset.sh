#!/usr/bin/env bash
# run_amrfinder_subset.sh - download assemblies for the chosen subset and run AMRFinderPlus on each.
#
# Run from the repository root, inside WSL, with the environment active:
#     conda activate salmonella-amr-tools
#     bash workflow/scripts/run_amrfinder_subset.sh pilot     # 20 isolates, a quick test
#     bash workflow/scripts/run_amrfinder_subset.sh full      # 300 isolates
#
# Big files (the genomes) go to $AMR_DATA_DIR (default ~/amr_data), OUTSIDE OneDrive.
# Small result tables go to data/interim/amrfinder_local/ (ignored by Git).
# Safe to stop and run again: finished genomes and results are skipped.
#
# Optional settings:  JOBS=2 (genomes at once)  THREADS=2 (threads per genome)  AMR_DATA_DIR=/some/path
set -euo pipefail

MODE="${1:-pilot}"
case "$MODE" in
  pilot) LIST="data/processed/subset_pilot_accessions.txt" ;;
  full)  LIST="data/processed/subset_300_accessions.txt" ;;
  *) echo "Usage: bash $0 [pilot|full]"; exit 2 ;;
esac
[ -f "$LIST" ] || { echo "Missing $LIST. In Windows PowerShell run: python workflow\\scripts\\02_select_subset.py"; exit 1; }

DATA_DIR="${AMR_DATA_DIR:-$HOME/amr_data}"
ASM_DIR="$DATA_DIR/assemblies"
OUT_DIR="data/interim/amrfinder_local"
JOBS="${JOBS:-2}"
THREADS="${THREADS:-2}"
mkdir -p "$ASM_DIR" "$OUT_DIR/logs"

for tool in datasets amrfinder unzip; do
  command -v "$tool" >/dev/null 2>&1 || { echo "Tool '$tool' not found. Did you run: conda activate salmonella-amr-tools ?"; exit 1; }
done

# record versions used for this run
{ date -u +"%Y-%m-%dT%H:%M:%SZ"; echo "amrfinder --version:"; amrfinder --version 2>&1 || true
  echo "amrfinder -l:"; amrfinder -l 2>&1 | head -n 6 || true
  echo "datasets --version:"; datasets --version 2>&1 || true; } > "$OUT_DIR/versions.txt"

# ---------------------------------------------------------------- 1. download
MISSING="$DATA_DIR/missing_${MODE}.txt"
while read -r acc; do
  [ -s "$ASM_DIR/$acc.fna" ] || echo "$acc"
done < "$LIST" > "$MISSING"

fetch_unpack() {   # $1 = zip file; moves every *_genomic.fna to $ASM_DIR/<accession>.fna
  local zip="$1" tmp="$1.dir"
  rm -rf "$tmp"; mkdir -p "$tmp"
  unzip -q -o "$zip" -d "$tmp"
  find "$tmp" -name "*_genomic.fna" | while read -r f; do
    mv "$f" "$ASM_DIR/$(basename "$(dirname "$f")").fna"
  done
  rm -rf "$tmp" "$zip"
}

if [ -s "$MISSING" ]; then
  echo ">>> Downloading $(wc -l < "$MISSING") assemblies from NCBI (one package) ..."
  PKG="$DATA_DIR/pkg_${MODE}.zip"
  if datasets download genome accession --inputfile "$MISSING" --include genome --filename "$PKG"; then
    fetch_unpack "$PKG"
  else
    echo "Batch download reported a problem; trying the missing ones one by one."
  fi
  # one-by-one retry for anything still missing
  while read -r acc; do
    [ -s "$ASM_DIR/$acc.fna" ] && continue
    echo "    retrying $acc"
    ONE="$DATA_DIR/one_${acc}.zip"
    if datasets download genome accession "$acc" --include genome --filename "$ONE"; then fetch_unpack "$ONE"; fi
  done < "$MISSING"
fi

FAILS="$OUT_DIR/download_failures_${MODE}.txt"
while read -r acc; do [ -s "$ASM_DIR/$acc.fna" ] || echo "$acc"; done < "$LIST" > "$FAILS"
echo ">>> Assemblies present: $(( $(wc -l < "$LIST") - $(wc -l < "$FAILS") )) of $(wc -l < "$LIST")  (failures listed in $FAILS)"

# ---------------------------------------------------------------- 2. AMRFinderPlus
run_one() {
  local acc="$1"
  local out="$OUT_DIR/${acc}.amrfinder.tsv"
  [ -s "$out" ] && return 0                       # already done
  [ -s "$ASM_DIR/$acc.fna" ] || return 0          # no genome to analyse
  if amrfinder -n "$ASM_DIR/$acc.fna" -O Salmonella --threads "$THREADS" -o "$out.tmp" 2> "$OUT_DIR/logs/${acc}.log"; then
    mv "$out.tmp" "$out"
  else
    rm -f "$out.tmp"; echo "$acc" >> "$OUT_DIR/amrfinder_failures_${MODE}.txt"
  fi
}
export -f run_one
export ASM_DIR OUT_DIR THREADS MODE

echo ">>> Running AMRFinderPlus: $JOBS genomes at a time, $THREADS threads each ..."
xargs -a "$LIST" -P "$JOBS" -I{} bash -c 'run_one {}'

DONE=$(find "$OUT_DIR" -maxdepth 1 -name "*.amrfinder.tsv" | wc -l)
echo ">>> Finished. Result tables so far: $DONE  (in $OUT_DIR)"
[ -f "$OUT_DIR/amrfinder_failures_${MODE}.txt" ] && echo "    Some genomes failed: see $OUT_DIR/amrfinder_failures_${MODE}.txt and logs/"
exit 0
