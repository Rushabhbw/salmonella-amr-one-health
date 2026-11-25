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

## J. Baseline version 2 and where its errors come from (21 Nov 2026, `reports/baseline_report.md`)
Primary dataset, rules version 2 (family-level symbols need every catalog entry to qualify; MULTIDRUG-class elements excluded). Sensitivity / specificity / MCC: tetracycline 0.986 / 0.993 / 0.978; ampicillin 0.970 / 0.995 / 0.971; ceftriaxone 0.960 / 0.999 / 0.972; ciprofloxacin 0.993 (interval 0.875-1.000) / 0.997 / 0.976. Ceftriaxone specificity rose from 0.990 to 0.999 after the blaTEM fix (major errors 77 to 10). The confusion counts were identical on the author's Windows computer and on a Linux test machine.
- Errors (false alarms / missed resistance): tetracycline 31 / 68; ampicillin 32 / 71; ceftriaxone 10 / 47; ciprofloxacin 24 / 4.
- Partial calls matter: a partial call for a linked element is present in 13% of the missed tetracycline-resistant isolates, 24% of ampicillin and 30% of ceftriaxone, against 5% of all labelled isolates, and in 0% of the false alarms. Counting partial calls raises sensitivity to 0.987 (tetracycline), 0.977 (ampicillin) and 0.972 (ceftriaxone) with specificity almost unchanged. Part of the missed resistance is therefore assembly fragmentation.
- Ampicillin: 57 of the 71 missed resistant isolates carry other resistance elements but no beta-lactamase: tet(A) 24, aph(3'')-Ib 24, aph(6)-Id 24, tet(B) 15, sul2 14. This is a multi-resistance pattern without a detected beta-lactamase. A model could learn such co-occurrence, which would be an association and not a mechanism; the grouped and transfer splits can test whether it generalizes.
- Errors are not spread evenly. The largest study (PRJNA292661) holds 54% of primary isolates but 34-39% of the errors. The older studies PRJNA183847/-48/-50/-51 (isolates from 2002-2010) have error rates of 1.6-3.6% for tetracycline and 1.5-2.2% for ampicillin, against 1.1% overall. PRJNA292666 (2014-2016) has an ampicillin error rate of 3.2% (18 missed, 9 false alarms) and a partial-call rate of 9.3% against 4.8% in the largest study. The cause (older tests, label noise or sequencing batch) is not known.
- Ciprofloxacin: isolates with a target-site mutation (gyrA, gyrB, parC, parE) are resistant in 579 of 584 labelled cases (99.1%); isolates with neither a mutation nor a qnr gene are resistant in 4 of 8,352; isolates with only a qnr gene are resistant in 20 of 39 labelled cases, with another 29 called intermediate and excluded. qnr genes caused 19 of the 24 false alarms. BioProject PRJNA1018843 (10 Nigerian isolates, 2018) has 10 of 10 errors: all are labelled susceptible and carry only qnrB19 or qnrS1.
- Consequence for machine learning: the baseline leaves little room (sensitivity 0.96-0.99, specificity at least 0.993 for the core antibiotics). Any gain would have to come from the missed resistance (partial calls, co-occurring genes). The planned comparison therefore uses the same fixed splits, a random split as a contrast, and paired cluster-bootstrap intervals for the difference from the baseline.

## K. Machine-learning design decisions (23 Nov 2026, written before any result was discussed)
- Question: can models learn something the rules miss, and does it transfer? Candidate sources of gain named in advance: partial calls (assembly fragmentation) and co-occurring genes (for example the ampicillin-resistant isolates without a detected beta-lactamase).
- Antibiotics: tetracycline, ampicillin, ceftriaxone. Ciprofloxacin stays descriptive (too few independent resistant groups).
- Fixed in advance: model types and hyperparameters, threshold 0.5, no class re-weighting, no tuning on test data. A hyperparameter search would need nested grouped cross-validation and is left as future work unless a model clearly underfits.
- Added a leaky `random5` split (92.4% of its test isolates have a close relative in training) so that the size of the inflation from unfair evaluation is measured and reported, never used to judge a model.
- Evidence standard: a model counts as better only if the paired cluster-bootstrap interval for the MCC difference (model minus rules) excludes zero.
- Shortcut test: a random forest without the catalog-linked elements shows what linkage alone can predict; the gap between random5 and the grouped splits shows whether that shortcut survives new lineages.
- Runtime on a single CPU core: about 5 minutes for 5 model variants x 3 antibiotics x 4 splits.

## L. Machine-learning results (24 Nov 2026; `reports/ml_report.md`)
Same isolates, splits and labels as the baseline; verdicts use the paired cluster-bootstrap interval for the MCC difference (model minus rules).
- Headline: no model beat the rules on any split. Fair splits (fold5, transfer_disjoint, time), 4 pre-specified models x 3 antibiotics x 3 splits = 36 comparisons: 27 worse than the rules, 9 no clear difference (all tetracycline), 0 better. The shortcut-only forest was worse in all 9. Differences in MCC range from 0.000 to -0.34 (the larger deficits belong to gradient boosting and the shortcut-only forest).
- Tetracycline: logistic regression reproduced the rules exactly on the transfer and time splits (MCC 0.976 and 0.987). The rules are essentially the ceiling for this antibiotic.
- Ampicillin and ceftriaxone: the main models were worse by roughly 0.04-0.17 MCC. Gradient boosting was least stable on the time split (ampicillin MCC 0.749, ceftriaxone 0.645, against 0.982 for the rules).
- Human transfer, ceftriaxone: logistic regression and random forest missed 9 of 45 resistant clinical isolates (sensitivity 0.800), the rules missed 1. Eight of the nine carried a catalog-linked gene that training had never or hardly seen: blaCTX-M-15 on 4 missed isolates (0 training isolates carry it), blaCMY-4 (0), blaSHV-12 (3), blaCTX-M-1 (1). Only 16 of the 5,158 training isolates carry any blaCTX-M gene, while 11 of the 45 resistant clinical isolates do. The catalog knows these genes; the models only know what poultry isolates showed them.
- An unfair evaluation hides this. With the leaky random split the four main models come within 0.00-0.01 MCC of the rules, but on the grouped folds the deficit is up to 0.15 MCC larger: inflation (random5 minus fold5, MCC) was 0.000-0.007 for tetracycline, 0.016-0.090 for ampicillin and 0.016-0.147 for ceftriaxone (gradient boosting: 0.963 leaky against 0.816 grouped).
- Shortcut test (random forest without the catalog-linked elements): MCC 0.30-0.66 on the grouped folds, 0.36-0.70 on the time and transfer splits. Resistance can partly be guessed from co-occurring genes (co-resistance). Most of that signal survives new lineages (grouped MCC is 91%, 76% and 62% of the leaky value for tetracycline, ampicillin and ceftriaxone), so the shortcut does not simply collapse on new lineages. It did not help the models beat the rules.
- What the models rely on: weight on catalog-linked elements was 55% (logistic regression) and 63% (random forest) for tetracycline, 63% and 78% for ampicillin, 42% and 74% for ceftriaxone. Top non-linked elements include sul2, aph(3'')-Ib, aph(6)-Id, tet(A) (for ampicillin and ceftriaxone), floR and fosA7.2.
- Partial-call indicators did not help the models (gradient boosting with and without them differ by 0.00-0.01 MCC), although counting partial calls helped the rules a little (section J).
- Exploratory hybrids (rules feature plus all elements; added after the first report): logistic regression copied the rules almost everywhere (difference 0.000, with small deficits for ceftriaxone on the transfer and time splits); gradient boosting was 16 of 18 times no different and twice worse, never better. Given the rules, the gene table adds nothing measurable.
- Caveats: hyperparameters were fixed and not tuned; about 60 comparisons were made, and several worse-than-rules intervals only just exclude 0 (for example ampicillin gradient boosting on fold5), but the direction is the same in every split and model; the transfer test has 258-270 isolates, mostly one 2014 study, so "no clear difference" there means wide intervals and not equivalence; labels are submitter-supplied.
- Reproducibility: after the first run, only the fold5 rows differed between the author's Windows computer and a Linux test machine, because fold assignment depended on the pandas version. It now uses a seeded numpy shuffle with a stable sort. Transfer, time and random5 results had been identical on both machines.

