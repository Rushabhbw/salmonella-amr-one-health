# Project Charter

**Title:** Resistome profiling and machine-learning prediction of antimicrobial resistance in livestock/food-associated and human-associated *Salmonella enterica*: a reproducible One Health workflow

**Author:** Rushabh Wakade  |  **Version:** 0.2 (scope locked after feasibility check)  |  **Date:** 16 Nov 2026

> Status: scope locked on 16 Nov 2026 using NCBI Isolates Browser counts recorded in `feasibility_notes.md`. Assembly availability is still to be verified in the data-curation step. No results are claimed in this document.

## 1. Background and problem statement
Antimicrobial resistance (AMR) spreads between animals, food and people, so it is studied under the One Health idea. Surveillance programmes (for example NARMS in the USA) test bacterial isolates in the laboratory and also sequence their genomes, which creates public datasets that pair genomes with resistance results. Rule-based tools such as AMRFinderPlus detect known resistance genes and mutations, but a detected gene does not always mean a resistant bacterium, and rules can miss unknown mechanisms. Machine learning can learn patterns from many labelled genomes, but results are easily inflated when closely related isolates appear in both training and test data. This project builds a transparent, reproducible workflow that compares a rule-based baseline with machine-learning models for predicting resistance, and compares resistance profiles between livestock/food-associated and human-associated isolates.

## 2. Research questions
- **RQ1.** How accurately do AMRFinderPlus-detected genes and mutations predict laboratory resistance phenotypes for the chosen antibiotics (rule-based baseline)?
- **RQ2.** Do machine-learning models using gene/mutation presence-absence features outperform the baseline when evaluated with a relatedness-aware split (not only a random split)?
- **RQ3.** How well does a model trained on non-clinical isolates (food, animal, environment) transfer to human clinical isolates, which are kept as a separate, cluster-independent test set?
- **RQ4 (secondary, descriptive).** Do resistome profiles differ between clinical and non-clinical isolates, and how much of any difference is explained by serotype, country, year and study?

## 3. Hypotheses (to be tested, not assumed)
- **H1.** The baseline performs well for antibiotics with well-known resistance genes and worse where mechanisms are complex or rare.
- **H2.** ML matches or improves on the baseline for some antibiotics, and performance under a random split is higher than under a relatedness-aware split (showing leakage inflation).
- **H3.** Performance on clinical isolates is lower than on held-out non-clinical isolates for at least some antibiotics, because of differences in serotype and source mix.
- **H4.** Apparent resistome differences between clinical and non-clinical isolates shrink after accounting for serotype and study.

## 4. Objectives and deliverables
1. A curated, documented dataset (metadata table with provenance and a frozen query date).
2. A Snakemake workflow: QC -> AMR detection -> feature table -> analysis.
3. Baseline vs ML evaluation with proper validation, class-imbalance handling and error analysis.
4. Statistical comparison of host categories.
5. Publication-quality figures, a project report, README, and a resume-ready description.

## 5. Data sources
- NCBI Pathogen Detection Isolates Browser (isolate metadata, AST phenotypes, AMRFinderPlus genotypes, assemblies). Organism: *Salmonella enterica* (non-typhoidal). (counts in `feasibility_notes.md`)
- Optional cross-check: BV-BRC curated AST collection.
- Note: AST values in NCBI are submitter-supplied and not vetted by NCBI staff beyond basic checks (see Limitations).

## 6. Methods overview (short)
- **Data curation:** filter to isolates with assembly + AST + usable host/source; harmonise labels (S/I/R -> resistant vs not resistant, rule written down); define host categories from the source field with a documented mapping.
- **Preprocessing/QC:** assembly quality checks (size, contig count, N50), removal of failed assemblies.
- **AMR detection:** AMRFinderPlus on assemblies (full set if feasible, otherwise precomputed NCBI results for the full set plus a local re-run on a subset to demonstrate and cross-check the pipeline).
- **Baseline:** gene/mutation -> phenotype rules; report accuracy, sensitivity, specificity, major errors and very major errors.
- **ML:** logistic regression, random forest, gradient boosting on presence/absence features; nested or grouped cross-validation; metrics suited to imbalance (PR-AUC, MCC, balanced accuracy); feature importance reported with caution.
- **Statistics:** host-category comparison with appropriate tests, multiple-testing correction (BH-FDR), effect sizes.

## 7. Success metrics
- Entire workflow reproducible from a clean clone with one documented command.
- Frozen, versioned metadata and environment files.
- Baseline and ML compared on the same held-out data, with a relatedness-aware split.
- All numbers in the report traceable to code and data.

## 7b. Scope decisions (v0.2)
- Organism: *Salmonella enterica* (NCBI taxgroup).
- Core antibiotics: tetracycline, ampicillin, ceftriaxone.
- Exploratory only: ciprofloxacin (resistant isolates are dominated by one SNP cluster, so evaluation is at cluster level and clearly caveated).
- Dropped: trimethoprim-sulfamethoxazole (too few independent resistant clusters; one cluster holds most resistant isolates).
- Host categories: "clinical" vs "environmental/other" as defined by NCBI isolation type. The clinical group is small (about 2-5% of AST isolates), so it is used mainly as an independent test set.
- Unit of independent evidence is the SNP cluster; splits are grouped by cluster.

## 8. Scope and non-goals
Non-goals: new laboratory work; predicting exact MIC values; multiple organisms; clinical/diagnostic claims; metagenomic resistome (possible future module).

## 9. Risks and mitigations
| Risk | Mitigation |
|---|---|
| Too few isolates in a class | Choose antibiotics only after counting; drop weak ones |
| Label noise / mixed AST methods | Record method info; filter where possible; state as limitation |
| Related isolates in train and test (leakage) | Group by serotype / SNP cluster; compare with random split |
| Inconsistent host/source text | Documented mapping table, manual review of unclear terms, exclude unclear |
| Sampling bias (country, year, outbreaks) | Report distribution; include as covariates or stratify |
| Compute/storage limits | Start with a subset; use precomputed AMR results for the full set |
| One huge SNP cluster dominates resistant isolates | Cluster-grouped splits; cap isolates per cluster; report cluster-level counts |
| Few clinical isolates (about 51-231 per class) | Use as test set only; bootstrap confidence intervals by cluster; state wide uncertainty |
| Some isolates have no assembly accession | Require an assembly for the local pipeline; report how many isolates are lost |
| Host group confounded with serotype, country and study | Include as covariates; stratified comparisons; state as limitation |

## 10. Timeline (about 5 weeks)
Week 1 charter, feasibility, repo/environment. Week 2 data curation and QC. Week 3 AMR detection and baseline. Week 4 ML and statistics. Week 5 figures, report, README, release.
