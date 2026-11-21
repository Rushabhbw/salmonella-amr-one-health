#!/usr/bin/env bash
# fetch_catalog.sh - copy the AMRFinderPlus Reference Gene Catalog (and release notes) into data/external/.
#
# Run from the repository root, inside WSL / Linux, with the environment active:
#     conda activate salmonella-amr-tools
#     bash workflow/scripts/fetch_catalog.sh
#
# The catalog gives each resistance element its antibiotic class. It is used by 05_map_elements_to_classes.py.
# First choice: the copy inside your installed AMRFinderPlus database. Fallback: download the latest copy from NCBI.
set -euo pipefail

command -v amrfinder >/dev/null 2>&1 || { echo "amrfinder not found. Did you run: conda activate salmonella-amr-tools ?"; exit 1; }
DEST="data/external"
mkdir -p "$DEST"

DB="$(amrfinder -l 2>&1 | sed -n 's/^Database directory: *//p' | head -n 1)"
if [ -n "$DB" ] && [ -f "$DB/ReferenceGeneCatalog.txt" ]; then
  cp "$DB/ReferenceGeneCatalog.txt" "$DEST/"
  for f in changes.txt version.txt; do [ -f "$DB/$f" ] && cp "$DB/$f" "$DEST/"; done
  echo "Copied from your local database: $DB"
else
  echo "Catalog not found in the local database folder; downloading the latest copy from NCBI ..."
  curl -fL -o "$DEST/ReferenceGeneCatalog.txt" \
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt"
fi

echo "Database version: $(cat "$DEST/version.txt" 2>/dev/null || echo 'not recorded (downloaded copy)')"
echo "Rows (including header): $(wc -l < "$DEST/ReferenceGeneCatalog.txt")"
echo "Columns in the header:"
head -n 1 "$DEST/ReferenceGeneCatalog.txt" | tr '\t' '\n' | sed 's/^/   /' | head -n 25
