# Data provenance log

Purpose: record exactly where every raw data file came from, so anyone can repeat the download.
NCBI Pathogen Detection is updated about daily, so counts can drift slightly. Always record the download date.

Source: NCBI Pathogen Detection, Isolates Browser (https://www.ncbi.nlm.nih.gov/pathogens/isolates/)
Download method: Download button, then tab-delimited (TSV) table.
Raw files are stored in `data/raw/` and are NOT committed to Git (see `.gitignore`).

| File name | Exact search query | Download date and time | Rows in file (excluding header) | Rows expected from browser (16 Nov 2026) | Notes |
|---|---|---|---|---|---|
| ceftriaxone_R.tsv | taxgroup_name:"Salmonella enterica" AND AST_phenotypes:ceftriaxone=R | | | 1,235 | |
| ceftriaxone_S.tsv | taxgroup_name:"Salmonella enterica" AND AST_phenotypes:ceftriaxone=S | | | 8,388 | |
| ampicillin_R.tsv | taxgroup_name:"Salmonella enterica" AND AST_phenotypes:ampicillin=R | | | 2,528 | |
| ampicillin_S.tsv | taxgroup_name:"Salmonella enterica" AND AST_phenotypes:ampicillin=S | | | 7,350 | |
| tetracycline_R.tsv | taxgroup_name:"Salmonella enterica" AND AST_phenotypes:tetracycline=R | | | 4,989 | |
| tetracycline_S.tsv | taxgroup_name:"Salmonella enterica" AND AST_phenotypes:tetracycline=S | | | 4,827 | |

All 51 available columns are ticked in "Choose columns" before every download (so all files have identical columns):

Organism group; Strain; Isolate identifiers; Serovar; Isolate; AMR genotypes core; IFSAC category; Source type; Host; K-mer group; Scientific name; AST phenotypes; Virulence genotypes; Stress genotypes; AMRFinderPlus version; PD Ref Gene Catalog version; AMRFinderPlus analysis type; BioProject; Collected by; Collection date; Lat/Lon; N50; Length; Contigs; Library layout; Platform; Run; SRA Center; SRA release date; Species TaxID; Outbreak; Host disease; Level; Method; PFGE primary enzyme pattern; PFGE secondary enzyme pattern; TaxID; WGS prefix; WGS accession; Create date; Location; Isolation source; Isolation type; Food origin; SNP cluster; Min-same; Min-diff; BioSample; Assembly; AMR genotypes; Computed types

First download (17 columns, ceftriaxone_R.tsv) was replaced by the 51-column version (ceftriaxone_R.tsv, saved as ceftriaxone_R_1.tsv at upload).

Known limitations: AST values are submitter-supplied and not vetted by NCBI beyond basic quality control.
