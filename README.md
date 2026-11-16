# Salmonella AMR: resistome profiling and machine-learning prediction (One Health)

**Status: work in progress (project set-up stage). No results yet.**

## Goal
Compare a rule-based baseline (AMRFinderPlus genes and mutations) with machine-learning models for predicting antimicrobial resistance in *Salmonella enterica*, using public NCBI Pathogen Detection data, and test how well models trained on food/animal/environmental isolates transfer to human clinical isolates.

## Research questions
See [`docs/00_project_charter.md`](docs/00_project_charter.md) for the full charter and [`docs/feasibility_notes.md`](docs/feasibility_notes.md) for data feasibility.

## Planned workflow
1. Data curation (metadata, labels, host categories)
2. Assembly QC
3. AMR detection (AMRFinderPlus)
4. Rule-based baseline
5. Machine learning with cluster-aware validation
6. Statistics and figures
7. Report

## Repository layout
| Folder | Purpose |
|---|---|
| `config/` | Workflow settings |
| `data/` | Raw, interim and processed data (large files are not committed) |
| `workflow/` | Snakemake rules and scripts |
| `notebooks/` | Exploratory analysis |
| `results/` | Figures and tables |
| `reports/` | Final report |
| `tests/` | Automated tests |
| `envs/` | Environment files |
| `docs/` | Charter, feasibility notes, documentation |

## Reproducibility
Environment file: `environment.yml` (versions will be pinned after the first successful install and test).

## Data sources
NCBI Pathogen Detection (Isolates Browser). Phenotype data are submitter-supplied and not vetted beyond basic quality control by NCBI.

## Author
Rushabh Wakade

## License
MIT (see `LICENSE`).
