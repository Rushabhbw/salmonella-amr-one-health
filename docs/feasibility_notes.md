# Feasibility Notes

**Prepared:** 16 Nov 2026 (literature-based part). The "Your counts" part must be filled from the live NCBI website by the project owner. Counts change over time, so record the query date.

## A. Evidence from the literature (checked via web search, not from live NCBI)
| Finding | Source |
|---|---|
| NARMS-based dataset of ~4,500 assembled non-typhoidal *Salmonella* genomes with MIC data for 15 antibiotics was used for ML | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10819212/ |
| NARMS isolates (2002-2016) from raw retail meat/poultry and livestock at slaughter used for XGBoost MIC prediction for 15 antibiotics | https://pmc.ncbi.nlm.nih.gov/articles/PMC6355527 |
| 640 *Salmonella* (43 serotypes; 104 human, 536 retail meat) with broth microdilution AST and WGS | https://pmc.ncbi.nlm.nih.gov/articles/PMC4997858 |
| BV-BRC curated AST collection: 67,817 genomes, 128 antibiotics | https://pmc.ncbi.nlm.nih.gov/articles/PMC8575023 |
| NCBI states AST values are submitter-supplied and not vetted beyond basic QC | https://www.ncbi.nlm.nih.gov/pathogens/ast |
| Isolates Browser supports `AST_phenotypes:` and `AMR_genotypes:` searches, table and assembly download | https://ncbi.nlm.nih.gov/core/assets/pathogens/files/HowTo/search_ast_data.pdf |

Conclusion so far: public data pairing *Salmonella* genomes with laboratory AST, from both food/animal and human sources, exists and has been used for ML. This supports feasibility but the exact numbers for our filters are NOT yet known.

## B. Counts from the NCBI Isolates Browser (read from screenshots by the project owner, 16 Nov 2026)
Search pattern: `taxgroup_name:"Salmonella enterica" AND AST_phenotypes:<antibiotic>=R` (or `=S`). Isolate counts are the "Displaying 1 - 20 of N" value in the Matched Isolates panel. Cluster counts are the "Matched Clusters" count. Intermediate (I) isolates are not included. Assembly availability not yet verified for each isolate.

| Antibiotic | R isolates | S isolates | R share | R clusters | S clusters | Biggest single cluster's share of R isolates |
|---|---|---|---|---|---|---|
| Ampicillin | 2,528 | 7,350 | 25.6% | 264 | 807 | 343 / 2,528 = 13.6% |
| Ceftriaxone | 1,235 | 8,388 | 12.8% | 159 | 795 | 309 / 1,235 = 25.0% |
| Ciprofloxacin | 758 | 8,828 | 7.9% | 74 | 846 | 538 / 758 = 71.0% |
| Tetracycline | 4,989 | 4,827 | 50.8% | 419 | 651 | 620 / 4,989 = 12.4% |
| Trimethoprim-sulfamethoxazole | 333 | 9,314 | 3.5% | 54 | 849 | 226 / 333 = 67.9% |

The biggest R cluster is the same one every time: PDS000265741.125 (17,704 isolates in total, almost all non-clinical).

Still to collect: clinical vs environmental/other split, countries/years, AST method information, assembly availability.

## B2. Host split (Isolation Type filter, read from screenshots, 16 Nov 2026)
"Clinical" = human clinical isolates. "Environmental/other" = food, environment and animal sources. In every search the two counts add up exactly to the totals in Section B (no unlabelled isolates).

| Antibiotic and class | Environmental/other | Clinical | Clinical share | Clinical clusters |
|---|---|---|---|---|
| Tetracycline R | 4,835 | 154 | 3.1% | 75 |
| Tetracycline S | 4,709 | 118 | 2.4% | 73 |
| Ampicillin R | 2,393 | 135 | 5.3% | 58 |
| Ampicillin S | 7,187 | 163 | 2.2% | 100 |
| Ceftriaxone R | 1,184 | 51 | 4.1% | 27 |
| Ceftriaxone S | 8,157 | 231 | 2.8% | 115 |

Resistant share within each group (sampling-dependent, not a prevalence estimate): tetracycline clinical 56.6% vs environmental/other 50.7%; ampicillin 45.3% vs 25.0%; ceftriaxone 18.1% vs 12.7%.

Observation to check: several top clinical rows show an empty Assembly column, so some isolates may have no assembly accession. To be quantified when the metadata table is downloaded.

## C. Findings and revised decision rule
- The unit of independent evidence is the SNP cluster, not the isolate. Isolates in one cluster are near-identical, so thousands of isolates can hide a small number of independent examples.
- Original rule (smaller class has a few hundred isolates) was too weak. Revised rule: keep an antibiotic as a core ML target only if the resistant class spans many independent clusters (at least about 100) and no single cluster holds more than about 30% of resistant isolates.
- Against this rule: tetracycline and ampicillin pass; ceftriaxone is borderline (159 clusters, 25% in one cluster) and kept as core because of clinical importance; ciprofloxacin (74 clusters, 71% in one) and trimethoprim-sulfamethoxazole (54 clusters, 68% in one) fail. Proposed: drop trimethoprim-sulfamethoxazole; treat ciprofloxacin as exploratory only, with cluster-level evaluation and an explicit limitation.
- Host balance: clinical isolates are only about 2-5% of AST isolates (about 270-300 per drug for the core antibiotics). A host-vs-host comparison with equal groups is not possible, so RQ3 was changed to a transfer test (train non-clinical, test clinical) and the comparison became a secondary descriptive question.
- Clinical isolates come from many countries, serotypes (e.g. Typhi, Typhimurium, Luciana, Saintpaul) and sources, so host group is confounded with serotype, place and study.
- NCBI AST values are submitter-supplied; label noise is a stated limitation.

## D. What the full-column test download showed (ceftriaxone R, 51 columns, 16 Nov 2026)
Checked with code on `ceftriaxone_R_1.tsv`: 1,235 rows (matches the browser count), same isolates as the 17-column file, and every row has ceftriaxone=R in its AST column.

- **Labels:** the `AST phenotypes` column lists every tested antibiotic for each isolate as `drug=S/R/I/NS/ND`. Labels for all drugs can therefore be built from this column. Planned rule: S = 0, R = 1; I, NS and ND are excluded for that drug and counted in the report. Most isolates have 14-15 drugs tested. The column has no MIC values or test method, so label noise cannot be examined further.
- **Time:** `Collection date` has values for 1,224 of 1,235 isolates (2002-2024). `Create date` is the NCBI record date, not the sampling date.
- **Study effect:** 22 BioProjects; the largest holds 598 of 1,235 isolates (48%). Candidate validation schemes: grouped by SNP cluster, by time, by BioProject.
- **Assembly QC columns:** `N50`, `Length`, `Contigs`. In this file: length 4.63-5.45 Mb, contigs 1-320, N50 31 kb to 5.1 Mb. Assemblers are mixed (SKESA 1,052; CLC 98; long-read assemblers about 64). Proposed filters will be set in the curation step.
- **AMR calls:** `AMR genotypes core` lists resistance genes only; `AMR genotypes` also lists stress/efflux genes (e.g. mdsA, mdsB). Call types: COMPLETE, POINT (mutations; 335 isolates), PARTIAL_END_OF_CONTIG, PARTIAL, MISTRANSLATION, HMM. Rule for PARTIAL calls to be decided and documented.
- **Reproducibility:** NCBI's stored AMR results come from two AMRFinderPlus software versions and five reference-gene-catalog (database) versions. After curation (9,783 isolates): software 4.2.5 for 9,242 and 4.2.7 for 541; catalog 2025-12-03.1 for 9,242, 2026-08-07.1 for 471, 2026-05-15.1 for 66, 2026-01-21.1 for 3, 2026-03-24.1 for 1. A local re-run with a single version on a subset is therefore justified.
- **Host and source fields:** `Host` is empty for 97% and `Source type` for 56%; `IFSAC category` for 66%. The free-text `Isolation source` (26 missing) will need a documented mapping to categories.
- **Duplicates:** one BioSample (strain CVM N45955, Dublin) appears twice under different isolate IDs; duplicates will be removed.
- **Geography:** 1,214 of 1,235 isolates are labelled USA (98%). Conclusions apply to US surveillance data.
- **Assemblies:** 14 of 1,235 have no assembly accession (5 of 51 clinical; 9 of 1,184 non-clinical).
- **Rough sanity check (not a result):** about 96.7% of these ceftriaxone-resistant isolates carry at least one gene from a hand-picked AmpC/ESBL list. Only the resistant class was examined, so false-alarm behaviour is unknown.

## E. Local AMRFinderPlus pilot (20 isolates, 17 Nov 2026)
- Tools in WSL (Ubuntu): AMRFinderPlus software 4.2.7, database 2026-10-29.1, NCBI datasets 18.41.1.
- Subset: 300 isolates chosen with seed 42 (100 human clinical, 200 non-clinical, at most 3 per SNP-cluster group within each host group); pilot = first 20 (6 clinical, 14 non-clinical).
- Download: 20 of 20 assemblies (29.5 MB package), 0 failures. 20 result tables produced with `amrfinder -n <assembly> -O Salmonella` (core elements only, no `--plus`).
- First result inspected (isolate PDT000161115.2, human clinical): our run reports 5 distinct elements, tet(B), aph(6)-Id, aph(3'')-Ib, sul2 and blaTEM-1, each on two contigs (10 table rows). NCBI's stored string lists the same 5 elements once each. Element sets are identical. Copy number is therefore not comparable, and features will be presence/absence.
- Our database (2026-10-29.1) is newer than the catalog used for 9,242 of the 9,783 curated isolates (2025-12-03.1), so small differences are expected and are not necessarily errors.
- Concordance over all 20 pilot isolates: 20 of 20 have identical sets of resistance elements (mean Jaccard 1.000). Their NCBI catalog versions were 2025-12-03.1 (17 isolates) and 2026-05-15.1 (3); none used our newer database. With 0 differences in 20 isolates, the true difference rate could still be as high as about 15% (rule of three: 3/20), so the 300-isolate run is needed before claiming the stored results are reproducible.

## F. Local AMRFinderPlus run on 300 isolates (17 Nov 2026)
- Our run: AMRFinderPlus 4.2.7, database 2026-10-29.1, `amrfinder -n <assembly> -O Salmonella` (core elements only). 300 of 300 genomes processed.
- Agreement with NCBI's stored results: 295 of 300 isolates have identical sets of resistance elements (98.3%); mean Jaccard 0.995. The 95% confidence interval for the difference rate is 0.7-3.8%. The subset is stratified (100 clinical, 200 non-clinical, at most 3 per cluster group within each host group), so it is not a random sample of all 9,783 isolates.
- The five differing isolates (all had catalog 2025-12-03.1 in NCBI's result):
  - dfrA51 only in NCBI's result (COMPLETE): PDT000161187.2 and PDT000161200.2 (clinical, serotype Potsdam, same BioProject, USA, 2014)
  - gar only in NCBI's result (HMM call): PDT000121333.2 (non-clinical)
  - mrx(A) only in our result (EXACTX): PDT000598976.1 (clinical, Stanley)
  - ramR_Q19Ter only in our result (POINTX): PDT000598974.1 (clinical, Senftenberg)
- None of the five differences involves a beta-lactamase or tetracycline-resistance gene, so no effect is expected on the three core antibiotics in this sample.
- 4 of the 5 differing isolates are clinical (4 of 100 clinical vs 1 of 200 non-clinical; Fisher exact p = 0.044). This is a post-hoc observation from small numbers and is treated as a hypothesis only.
- Explanation, partly verified (21 Nov 2026): dfrA51 does not exist in the current Reference Gene Catalog (database 2026-10-29.1), while 63 primary-dataset isolates carry it in NCBI's 2025-12-03.1 results. This explains the two dfrA51 differences. The other three (gar, mrx(A), ramR_Q19Ter) are not yet explained; the latest `changes.txt` contains no mention of dfrA51, mrx(A) or ramR, so it does not help.
- Decision: features come from NCBI's stored results. The primary dataset is the 9,242 isolates with catalog 2025-12-03.1 (software 4.2.5); the other 541 (software 4.2.7, four newer catalogs) are used only in a sensitivity analysis. Primary-dataset clinical resistant isolates: tetracycline 145, ampicillin 129, ceftriaxone 45, ciprofloxacin 11.
- Cost of restricting: ciprofloxacin-resistant isolates fall from 78 to 29 SNP-cluster groups, so ciprofloxacin is analysed descriptively only (no machine learning). Core antibiotics lose little (resistant groups: tetracycline 528 to 501, ampicillin 331 to 313, ceftriaxone 205 to 198).
- Re-running all 9,783 genomes locally was rejected: the 300-isolate run took roughly an hour, so the full set would take about 30 hours and about 15 GB of downloads, for a small expected gain.

## G. Feature table (17-19 Nov 2026)
Built from NCBI's stored AMRFinderPlus core results with `04_build_features.py` (report: `reports/feature_report.md`).
- Call types in the stored results: COMPLETE 23,412; POINT 896; PARTIAL_END_OF_CONTIG 655; PARTIAL 356; MISTRANSLATION 29; HMM 6. Features use COMPLETE and POINT only; weak calls are kept in a separate table.
- 149 elements have at least one COMPLETE/POINT call; 9 more occur only with weak calls and are not features. In the primary dataset, 127 elements occur in at least 1 SNP-cluster group, 77 in at least 3, 38 in at least 10. Median 2 elements per isolate (maximum 21).
- 2,936 isolates (30.0%) have no element. Among them almost all are susceptible: resistant share 1.2% for tetracycline (36 of 2,910 with a label), 0.5% for ampicillin, 0.3% for ceftriaxone and 0.6% for ciprofloxacin. So known elements explain nearly all resistance; the few resistant isolates without any element are unexplained (unknown mechanism, label noise or a database gap).
- The point mutation gyrA_D87Y occurs in 722 primary-dataset isolates but only 19 independent groups, which again shows how few independent examples back ciprofloxacin resistance.
- No feature is filtered at this stage; rare-feature filtering will be done inside cross-validation.

## H. Evaluation design findings (20 Nov 2026, `reports/splits_report.md`)
- 64 SNP-cluster groups contain both clinical and non-clinical isolates. 162 of the 270 clinical isolates (60%) sit in such groups, so a naive "train on non-clinical, test on clinical" split would be heavily leaky. The clean version removes 3,814 of 8,972 non-clinical isolates from training (43%), leaving 5,158 for training. Resistant/susceptible counts in the clean transfer split: tetracycline train 2,706/2,427, test 145/113; ampicillin 1,176/3,982 and 129/141; ceftriaxone 639/4,434 and 45/224.
- The clinical test set is concentrated: 222 of 270 isolates come from BioProject PRJNA287785 and 204 were collected in 2014 (232 clinical isolates have a collection year). The transfer result therefore mostly reflects transfer to one 2014 human collection.
- Time split (train up to 2016: 5,489 isolates; test 2017 onward: 3,715): 3,253 of the 3,715 test isolates (87.6%) share a group with training. Ciprofloxacin has 17 resistant isolates in training and 584 in the test period, so it is not evaluated this way.
- Grouped 5-fold: 1,848-1,849 isolates per fold, 0 test isolates with a relative in training. Fewest resistant isolates in any fold: tetracycline 535, ampicillin 250, ceftriaxone 126.

## I. Element mapping and first baseline run (21 Nov 2026)
Mapping of the 149 elements with the Reference Gene Catalog (database 2026-10-29.1, downloaded from NCBI): 58 matched on allele, 81 on gene family, 10 not found, 10 with more than one class/subclass entry. Not found: dfrA51 (63 primary isolates, trimethoprim, not one of our antibiotics), eight ramR frameshift or nonsense mutations and aac(3)-C322 (0 primary isolates each).

First baseline run (rules version 1, any matching catalog entry counted; primary dataset; all four antibiotics): tetracycline sensitivity 0.986 and specificity 0.992 (9 elements); ampicillin 0.970 and 0.995 (24); ceftriaxone 0.961 and 0.990 (13); ciprofloxacin 0.993 and 0.996 (22, wide intervals because of few independent resistant groups).

Problems found in that run, and what was changed (version 2 of the rules):
- Ceftriaxone: 66 of its 77 major errors came from the family-level symbol blaTEM. A family symbol stands for many alleles, and TEM-1 does not break cephalosporins, but one catalog entry for the family does list CEPHALOSPORIN. Rule changed: an element counts only if every catalog entry it matches satisfies the rule. Side effect: the family-level blaCTX-M symbol (10 isolates, 1 group) no longer counts for ceftriaxone.
- Tetracycline and ciprofloxacin: ramR_T18P (class MULTIDRUG; 7 isolates) was present in 5 susceptible isolates for tetracycline and 6 for ciprofloxacin. oqxA and oqxB were also counted. Rule changed: MULTIDRUG-class elements are not drug-specific evidence; they are kept in a sensitivity row.
- Both changes were made after seeing the first report, so they are analyst choices informed by the same data. The effect is small, the original and alternative rules are both reported, and the evidence rules still come from NCBI's catalog and not from the labels.
- Ciprofloxacin: qnr genes (qnrB19 13, qnrS1 6) account for most major errors. Plasmid-mediated qnr genes lower susceptibility only a little, so an isolate called susceptible can carry one. Ciprofloxacin stays descriptive.
- Unexplained resistance (first run): tetracycline 68 resistant isolates called susceptible (34 with no element, 34 with other elements); ampicillin 71; ceftriaxone 46; ciprofloxacin 4.
New diagnostics added: sensitivity rows counting partial calls and MULTIDRUG-class elements, a count of missed resistant isolates that have a partial call for a linked element, and a table of the BioProjects with the most errors.

