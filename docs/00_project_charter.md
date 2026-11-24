# Project Charter

**Title:** Resistome profiling and machine-learning prediction of antimicrobial resistance in livestock/food-associated and human-associated *Salmonella enterica*: a reproducible One Health workflow

**Author:** Rushabh Wakade  |  **Version:** 0.6 (machine-learning design fixed before any model results)  |  **Date:** 16 Nov 2026

> Status: scope locked on 16 Nov 2026 and updated on 17 Nov 2026 with the curated dataset (9,783 isolates) and the 300-isolate local AMRFinderPlus check. Counts are in `reports/curation_report.md` and `docs/feasibility_notes.md`. No model results are claimed in this document.

## 1. Background and problem statement
Antimicrobial resistance (AMR) spreads between animals, food and people, so it is studied under the One Health idea. Surveillance programmes (for example NARMS in the USA) test bacterial isolates in the laboratory and also sequence their genomes, which creates public datasets that pair genomes with resistance results. Rule-based tools such as AMRFinderPlus detect known resistance genes and mutations, but a detected gene does not always mean a resistant bacterium, and rules can miss unknown mechanisms. Machine learning can learn patterns from many labelled genomes, but results are easily inflated when closely related isolates appear in both training and test data. This project builds a transparent, reproducible workflow that compares a rule-based baseline with machine-learning models for predicting resistance, and compares resistance profiles between livestock/food-associated and human-associated isolates.

## 2. Research questions
- **RQ1.** How accurately do AMRFinderPlus-detected genes and mutations predict laboratory resistance phenotypes for the chosen antibiotics (rule-based baseline)?
- **RQ2.** Do machine-learning models using gene/mutation presence-absence features outperform the baseline when evaluated with a relatedness-aware split (not only a random split)?
- **RQ3.** How well does a model trained on non-clinical isolates (food, animal, environment) transfer to human clinical isolates, which are kept as a separate, cluster-independent test set? (In the primary dataset 222 of the 270 clinical isolates come from one BioProject and 204 were collected in 2014, so this is largely transfer to one human collection. This is stated as a limitation.)
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
- **AMR detection:** NCBI's stored AMRFinderPlus core results for all curated isolates (used as features), plus a local re-run with a single software and database version on a 300-isolate subset to check agreement (done: 295 of 300 identical).
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
- Ciprofloxacin: descriptive, rule-based baseline only, no machine learning. After curation there are only 12 resistant human isolates (7 groups) and, in the primary dataset, resistant isolates come from only 29 SNP-cluster groups (23 non-clinical), which is too few for cluster-grouped validation or a human transfer test. It is analysed on all isolates and reported with that caveat.
- Curated dataset: 9,783 isolates (288 human clinical, 9,495 non-clinical; about 85% of the non-clinical group is poultry; about 97% are from the USA).
- Human clinical resistant isolates available for RQ3 (primary dataset): tetracycline 145, ampicillin 129, ceftriaxone 45 (from 25 groups). Confidence intervals for ceftriaxone will be wide and are reported as such.
- Labels come from the `AST phenotypes` column: S = 0, R = 1; I, NS, ND and missing values are excluded for that drug.
- Evaluation splits (fixed once in `07_make_splits.py`, seed 42, primary dataset only): (1) grouped 5-fold cross-validation by SNP-cluster group (no group in both train and test; new-lineage performance); (2) clean human transfer test: train on the 5,158 non-clinical isolates that share no group with a clinical isolate, test on the 270 clinical isolates; the naive version is reported only to show its leak (162 of 270 clinical isolates, 60%, have a relative in the training data); (3) time split: train up to 2016, test later (87.6% of test isolates share a group with training, so this is a forward-prediction check and not a new-lineage test).
- Rule-based baseline: an isolate is predicted resistant if it carries an element that NCBI's Reference Gene Catalog links to the antibiotic (tetracycline: TETRACYCLINE class; ampicillin: BETA-LACTAM class; ceftriaxone: BETA-LACTAM class with CEPHALOSPORIN subclass; ciprofloxacin: QUINOLONE). An element counts only if every catalog entry it matches satisfies the rule (a family-level symbol such as blaTEM is not cephalosporin evidence), and MULTIDRUG-class elements are not drug-specific evidence; both choices were refined after a first run (see `docs/feasibility_notes.md`, section I) and their alternatives are reported as sensitivity rows. Reported with sensitivity, specificity, MCC, very major and major error rates and cluster-bootstrap 95% intervals. Machine-learning models must beat the baseline on the same splits to be reported as an improvement.
- Ciprofloxacin cannot be evaluated with the time split (17 resistant isolates in training, 584 in the test period).
- Machine learning (`08_train_evaluate_models.py`; tetracycline, ampicillin and ceftriaxone only): logistic regression, random forest and histogram gradient boosting on element presence/absence, a gradient-boosting variant that also sees partial-call indicators, and a shortcut test (random forest without the elements the catalog links to the antibiotic). Hyperparameters are fixed in advance, nothing is tuned on test data, threshold 0.5. Models are judged on the grouped 5-fold, clean human-transfer and time splits; a deliberately leaky random 5-fold is shown only to measure how much unfair evaluation inflates results. Each model is compared with the rule-based baseline on the same test isolates with a paired cluster-bootstrap interval for the difference in MCC; a model is reported as an improvement only if that interval excludes zero. "No clear difference" and "worse than the rules" are reported as results.
- Features: presence/absence of AMRFinderPlus core elements from NCBI's stored results. Primary dataset = the 9,242 of 9,783 isolates (270 of 288 clinical) whose stored results come from one software and catalog version (4.2.5, 2025-12-03.1); the remaining 541 isolates are used only in a sensitivity analysis. A local re-run (software 4.2.7, database 2026-10-29.1) on 300 isolates gave identical element sets for 295 (98.3%; 95% CI for the difference rate 0.7-3.8%).
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
