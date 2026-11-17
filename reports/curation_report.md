# Metadata curation report

Generated: 2026-11-17 09:31  
Script: `workflow/scripts/01_curate_metadata.py`  
Raw files: ampicillin_R.tsv, ampicillin_S.tsv, ceftriaxone_R.tsv, ceftriaxone_S.tsv, tetracycline_R.tsv, tetracycline_S.tsv

## 1. Records kept after each step

| Step | Isolates remaining |
|---|---|
| Rows read from all raw files | 29317 |
| Unique isolates (same isolate in several files counted once) | 9882 |
| After removing duplicate BioSamples | 9881 |
| After removing duplicate assembly accessions | 9881 |
| After removing isolates without an assembly accession | 9783 |
| After assembly quality sanity limits | 9783 |

Removed isolates by reason:

| reason | n |
|---|---|
| no_assembly_accession | 98 |
| duplicate_biosample (another record of the same BioSample was kept) | 1 |

Quality limits applied: genome length 4,300,000-5,700,000 bp, at most 500 contigs, N50 at least 20,000 bp. They are sanity limits chosen from typical *Salmonella* genome size.

## 2. File-versus-label check

Each raw file should contain only the class in its name according to its own AST column.

No problems found: every row agrees with its file name.

## 3. Labels per antibiotic (after curation)

| Antibiotic | Group | R (label 1) | S (label 0) | Excluded (I/NS/ND/missing) | R share | Clusters containing R |
|---|---|---|---|---|---|---|
| tetracycline | all | 4944 | 4775 | 64 | 50.9% | 528 |
| tetracycline | non_clinical | 4797 | 4662 | 36 | 50.7% | 488 |
| tetracycline | clinical | 147 | 113 | 28 | 56.5% | 72 |
| ampicillin | all | 2494 | 7289 | 0 | 25.5% | 331 |
| ampicillin | non_clinical | 2363 | 7132 | 0 | 24.9% | 296 |
| ampicillin | clinical | 131 | 157 | 0 | 45.5% | 60 |
| ceftriaxone | all | 1220 | 8322 | 241 | 12.8% | 205 |
| ceftriaxone | non_clinical | 1174 | 8095 | 226 | 12.7% | 192 |
| ceftriaxone | clinical | 46 | 227 | 15 | 16.8% | 25 |
| ciprofloxacin (exploratory) | all | 735 | 8771 | 277 | 7.7% | 78 |
| ciprofloxacin (exploratory) | non_clinical | 723 | 8555 | 217 | 7.8% | 72 |
| ciprofloxacin (exploratory) | clinical | 12 | 216 | 60 | 5.3% | 7 |

## 4. Independent groups (SNP clusters)

- Groups: 1203 (including 298 isolates without a SNP cluster, each its own group)
- Largest group: 671 isolates (6.9%); top 5 groups: 23.9%

## 5. Source, place, time, study

Source categories:

| source_category | n |
|---|---|
| chicken | 4667 |
| turkey | 3260 |
| pork_swine | 836 |
| beef_cattle | 494 |
| human_clinical | 288 |
| poultry_unspecified | 177 |
| other_or_unclear | 47 |
| animal_or_farm_sample | 8 |
| unknown | 6 |

Countries (top 8):

| country | n |
|---|---|
| USA | 9508 |
| Mexico | 136 |
| Germany | 73 |
| unknown | 38 |
| Nigeria | 13 |
| South Korea | 8 |
| Pakistan | 3 |
| Oman | 3 |

Collection years (grouped):

| years | n |
|---|---|
| 2002-2009 | 2528 |
| 2010-2014 | 2179 |
| 2015-2018 | 2472 |
| 2019-2021 | 2558 |
| 2022+ | 8 |

Missing collection year: 38

BioProjects (top 5 of 31):

| bioproject | n |
|---|---|
| PRJNA292661 | 5016 |
| PRJNA292666 | 837 |
| PRJNA183851 | 781 |
| PRJNA183850 | 592 |
| PRJNA242614 | 526 |

Assemblers (top 5):

| assembler | n |
|---|---|
| SKESA v. 2.2 | 9086 |
| CLC Genomics Workbench v. 6.0.2 | 528 |
| SMRTLink Microbial Assembly v. SMRTLink9.0 | 71 |
| HGAP v. 4 | 44 |
| HGAP v. 3 | 11 |

AMRFinderPlus versions used by NCBI:

| version | n |
|---|---|
| 4.2.5 | 9242 |
| 4.2.7 | 541 |

## 6. Source text that could not be grouped confidently (top 15)

| isolation_source_clean | source_category | n_isolates |
|---|---|---|
| poultry | poultry_unspecified | 136 |
| whole - cut in lab | other_or_unclear | 9 |
| wings | poultry_unspecified | 9 |
| breast | poultry_unspecified | 9 |
| thighs | poultry_unspecified | 8 |
| mixed parts | other_or_unclear | 6 |
| missing | unknown | 6 |
| breasts | poultry_unspecified | 3 |
| legs | poultry_unspecified | 3 |
| ground | other_or_unclear | 3 |
| dead-in-shell embryo | other_or_unclear | 3 |
| gb-ground | other_or_unclear | 2 |
| cb-thighs | poultry_unspecified | 2 |
| shrimp (sh) | other_or_unclear | 2 |
| gt-ground | other_or_unclear | 2 |

## 7. Known limitations

- Antibiotic labels are submitter-supplied and not vetted by NCBI beyond basic checks; the test method is not in the table.
- Isolates are heavily clustered, so isolate counts overstate independent evidence; use `cluster_group` for splitting.
- Most isolates are from the USA.
- The 'clinical' group is small (human isolates only).
- `source_category` comes from text rules (`data/processed/source_mapping.tsv`) and should be reviewed.
