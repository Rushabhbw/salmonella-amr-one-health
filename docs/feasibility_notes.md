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
