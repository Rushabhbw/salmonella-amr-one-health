# Rule-based baseline report

Generated: 2026-11-23 09:54  
Script: `workflow/scripts/06_baseline_rules.py`

## Method

- Prediction: RESISTANT if the isolate carries at least one element that NCBI links to the antibiotic (`config/element_classes.tsv`), otherwise SUSCEPTIBLE. No training.
- Labels: S = 0, R = 1; intermediate and other calls are excluded for that antibiotic.
- Very major error rate = resistant isolates predicted susceptible (FN / all R). Major error rate = susceptible isolates predicted resistant (FP / all S).
- Confidence intervals: cluster bootstrap, 2000 resamples of whole SNP-cluster groups, seed 42. 'Groups' is the number of independent groups behind each row.
- The last two rows of each table are sensitivity analyses on the primary dataset: they change the evidence rule, not the data.
- Ciprofloxacin is descriptive only (see the charter): few independent resistant groups.

## tetracycline

Elements used as evidence (6): tet(A), tet(B), tet(C), tet(G), tet(M), tet(D)

| Group | R / S isolates | Groups | TP/FN/TN/FP | Sensitivity | Specificity | MCC | Very major error (R called S) | Major error (S called R) |
|---|---|---|---|---|---|---|---|---|
| all (primary dataset) | 4718 / 4476 | 1113 | 4650/68/4445/31 | 0.986 (0.977-0.990) | 0.993 (0.989-0.996) | 0.978 (0.973-0.982) | 1.4% (1.0-2.3%) | 0.7% (0.4-1.1%) |
| non-clinical (primary) | 4573 / 4363 | 1039 | 4506/67/4334/29 | 0.985 (0.977-0.990) | 0.993 (0.989-0.996) | 0.979 (0.973-0.982) | 1.5% (1.0-2.3%) | 0.7% (0.4-1.1%) |
| clinical (primary) | 145 / 113 | 137 | 144/1/111/2 | 0.993 (0.976-1.000) | 0.982 (0.952-1.000) | 0.976 (0.942-1.000) | 0.7% (0.0-2.4%) | 1.8% (0.0-4.8%) |
| all isolates incl. non-primary (sensitivity check) | 4944 / 4775 | 1187 | 4865/79/4739/36 | 0.984 (0.975-0.989) | 0.992 (0.988-0.996) | 0.976 (0.971-0.980) | 1.6% (1.1-2.5%) | 0.8% (0.4-1.2%) |
| primary, also counting MULTIDRUG-class elements | 4718 / 4476 | 1113 | 4650/68/4440/36 | 0.986 (0.977-0.990) | 0.992 (0.987-0.995) | 0.977 (0.971-0.981) | 1.4% (1.0-2.3%) | 0.8% (0.5-1.3%) |
| primary, also counting partial calls | 4718 / 4476 | 1113 | 4658/60/4438/38 | 0.987 (0.978-0.992) | 0.992 (0.987-0.995) | 0.979 (0.973-0.983) | 1.3% (0.8-2.2%) | 0.8% (0.5-1.3%) |

Primary dataset errors: 31 major (S called R), 68 very major (R called S).

Elements most often present in susceptible isolates wrongly called resistant: tet(A) (20), tet(B) (9), tet(C) (1), tet(M) (1)

Resistant isolates called susceptible: 34 carry no element at all; 34 carry elements, but none linked to this antibiotic. 8 of the 68 have a partial call for an element linked to this antibiotic (the gene may be broken across contigs).

BioProjects with the most errors (overall error rate 1.1%). A rate far above that points to a study-specific problem such as label noise or a different test method:

| bioproject | labelled | errors | error_rate | share_of_all_errors |
|---|---|---|---|---|
| PRJNA292661 | 4987 | 38 | 0.8% | 38.4% |
| PRJNA183851 | 771 | 17 | 2.2% | 17.2% |
| PRJNA183847 | 394 | 14 | 3.6% | 14.1% |
| PRJNA292666 | 835 | 7 | 0.8% | 7.1% |
| PRJNA183848 | 436 | 7 | 1.6% | 7.1% |

## ampicillin

Elements used as evidence (24): blaTEM-1, blaCMY-2, blaCTX-M-65, blaHER-3, blaTEM, blaCARB-2, blaCMY, blaCTX-M, blaCTX-M-15, blaCTX-M-1, blaHER, blaSHV-12, blaCTX-M-14, blaTEM-135, blaCTX-M-55, blaLAP-2, blaOXA-1, blaCMY-4, blaKPC-2, blaLAP, blaNDM-1, blaOXA-9, blaPSE, blaSHV-2A

| Group | R / S isolates | Groups | TP/FN/TN/FP | Sensitivity | Specificity | MCC | Very major error (R called S) | Major error (S called R) |
|---|---|---|---|---|---|---|---|---|
| all (primary dataset) | 2363 / 6879 | 1121 | 2292/71/6847/32 | 0.970 (0.953-0.979) | 0.995 (0.993-0.997) | 0.971 (0.959-0.977) | 3.0% (2.1-4.7%) | 0.5% (0.3-0.7%) |
| non-clinical (primary) | 2234 / 6738 | 1042 | 2164/70/6706/32 | 0.969 (0.951-0.979) | 0.995 (0.993-0.997) | 0.969 (0.957-0.977) | 3.1% (2.1-4.9%) | 0.5% (0.3-0.7%) |
| clinical (primary) | 129 / 141 | 143 | 128/1/141/0 | 0.992 (0.969-1.000) | 1.000 (1.000-1.000) | 0.993 (0.973-1.000) | 0.8% (0.0-3.1%) | 0.0% (0.0-0.0%) |
| all isolates incl. non-primary (sensitivity check) | 2494 / 7289 | 1203 | 2414/80/7253/36 | 0.968 (0.952-0.977) | 0.995 (0.993-0.997) | 0.969 (0.958-0.975) | 3.2% (2.3-4.8%) | 0.5% (0.3-0.7%) |
| primary, also counting MULTIDRUG-class elements | 2363 / 6879 | 1121 | 2292/71/6843/36 | 0.970 (0.953-0.979) | 0.995 (0.992-0.997) | 0.969 (0.957-0.977) | 3.0% (2.1-4.7%) | 0.5% (0.3-0.8%) |
| primary, also counting partial calls | 2363 / 6879 | 1121 | 2308/55/6835/44 | 0.977 (0.964-0.984) | 0.994 (0.991-0.996) | 0.972 (0.962-0.977) | 2.3% (1.6-3.6%) | 0.6% (0.4-0.9%) |

Primary dataset errors: 32 major (S called R), 71 very major (R called S).

Elements most often present in susceptible isolates wrongly called resistant: blaTEM-1 (16), blaCMY-2 (6), blaHER-3 (4), blaCARB-2 (3), blaTEM (3), blaCMY (1)

Resistant isolates called susceptible: 14 carry no element at all; 57 carry elements, but none linked to this antibiotic. 16 of the 71 have a partial call for an element linked to this antibiotic (the gene may be broken across contigs).

BioProjects with the most errors (overall error rate 1.1%). A rate far above that points to a study-specific problem such as label noise or a different test method:

| bioproject | labelled | errors | error_rate | share_of_all_errors |
|---|---|---|---|---|
| PRJNA292661 | 5000 | 35 | 0.7% | 34.0% |
| PRJNA292666 | 837 | 27 | 3.2% | 26.2% |
| PRJNA183850 | 582 | 13 | 2.2% | 12.6% |
| PRJNA183851 | 780 | 12 | 1.5% | 11.7% |
| PRJNA219491 | 199 | 5 | 2.5% | 4.9% |

## ceftriaxone

Elements used as evidence (11): blaCMY-2, blaCTX-M-65, blaCMY, blaCTX-M-15, blaCTX-M-1, blaSHV-12, blaCTX-M-14, blaCTX-M-55, blaOXA-1, blaCMY-4, blaSHV-2A

| Group | R / S isolates | Groups | TP/FN/TN/FP | Sensitivity | Specificity | MCC | Very major error (R called S) | Major error (S called R) |
|---|---|---|---|---|---|---|---|---|
| all (primary dataset) | 1177 / 7974 | 1094 | 1130/47/7964/10 | 0.960 (0.925-0.977) | 0.999 (0.998-1.000) | 0.972 (0.954-0.982) | 4.0% (2.3-7.5%) | 0.1% (0.0-0.2%) |
| non-clinical (primary) | 1132 / 7750 | 1015 | 1086/46/7740/10 | 0.959 (0.917-0.976) | 0.999 (0.998-1.000) | 0.971 (0.949-0.981) | 4.1% (2.4-8.3%) | 0.1% (0.0-0.2%) |
| clinical (primary) | 45 / 224 | 142 | 44/1/224/0 | 0.978 (0.914-1.000) | 1.000 (1.000-1.000) | 0.987 (0.949-1.000) | 2.2% (0.0-8.6%) | 0.0% (0.0-0.0%) |
| all isolates incl. non-primary (sensitivity check) | 1220 / 8322 | 1133 | 1171/49/8312/10 | 0.960 (0.928-0.977) | 0.999 (0.998-1.000) | 0.972 (0.955-0.982) | 4.0% (2.3-7.2%) | 0.1% (0.0-0.2%) |
| primary, also counting MULTIDRUG-class elements | 1177 / 7974 | 1094 | 1130/47/7964/10 | 0.960 (0.925-0.977) | 0.999 (0.998-1.000) | 0.972 (0.954-0.982) | 4.0% (2.3-7.5%) | 0.1% (0.0-0.2%) |
| primary, also counting partial calls | 1177 / 7974 | 1094 | 1144/33/7963/11 | 0.972 (0.947-0.984) | 0.999 (0.998-0.999) | 0.978 (0.964-0.986) | 2.8% (1.6-5.3%) | 0.1% (0.1-0.2%) |

Primary dataset errors: 10 major (S called R), 47 very major (R called S).

Elements most often present in susceptible isolates wrongly called resistant: blaCMY-2 (8), blaCMY (1), blaCTX-M-65 (1)

Resistant isolates called susceptible: 8 carry no element at all; 39 carry elements, but none linked to this antibiotic. 14 of the 47 have a partial call for an element linked to this antibiotic (the gene may be broken across contigs).

BioProjects with the most errors (overall error rate 0.6%). A rate far above that points to a study-specific problem such as label noise or a different test method:

| bioproject | labelled | errors | error_rate | share_of_all_errors |
|---|---|---|---|---|
| PRJNA292661 | 4994 | 22 | 0.4% | 38.6% |
| PRJNA292666 | 837 | 7 | 0.8% | 12.3% |
| PRJNA183851 | 780 | 7 | 0.9% | 12.3% |
| PRJNA183850 | 582 | 6 | 1.0% | 10.5% |
| PRJNA183847 | 395 | 6 | 1.5% | 10.5% |

## ciprofloxacin

Elements used as evidence (21): gyrA_D87Y, qnrB19, gyrA_S83F, gyrA_D87N, qnrS1, gyrA_D87G, gyrA_S83Y, qnrA1, parC_S80I, qnrS2, aac(6')-Ib-cr5, gyrA_D82N, qnrB2, qnrS13, gyrB_S464Y, oqxA, oqxB, parE_A499T, qnrB1, qnrB77, qnrS4

| Group | R / S isolates | Groups | TP/FN/TN/FP | Sensitivity | Specificity | MCC | Very major error (R called S) | Major error (S called R) |
|---|---|---|---|---|---|---|---|---|
| all (primary dataset) | 603 / 8372 | 1086 | 599/4/8348/24 | 0.993 (0.875-1.000) | 0.997 (0.995-0.999) | 0.976 (0.721-0.992) | 0.7% (0.0-12.5%) | 0.3% (0.1-0.5%) |
| non-clinical (primary) | 592 / 8168 | 1033 | 588/4/8151/17 | 0.993 (0.867-1.000) | 0.998 (0.996-0.999) | 0.981 (0.749-0.995) | 0.7% (0.0-13.3%) | 0.2% (0.1-0.4%) |
| clinical (primary) | 11 / 204 | 112 | 11/0/197/7 | 1.000 (1.000-1.000) | 0.966 (0.933-0.988) | 0.768 (0.508-0.914) | 0.0% (0.0-0.0%) | 3.4% (1.2-6.7%) |
| all isolates incl. non-primary (sensitivity check) | 735 / 8771 | 1170 | 708/27/8744/27 | 0.963 (0.777-0.988) | 0.997 (0.994-0.999) | 0.960 (0.776-0.984) | 3.7% (1.2-22.3%) | 0.3% (0.1-0.6%) |
| primary, also counting MULTIDRUG-class elements | 603 / 8372 | 1086 | 599/4/8342/30 | 0.993 (0.875-1.000) | 0.996 (0.994-0.998) | 0.971 (0.689-0.990) | 0.7% (0.0-12.5%) | 0.4% (0.2-0.6%) |
| primary, also counting partial calls | 603 / 8372 | 1086 | 599/4/8348/24 | 0.993 (0.875-1.000) | 0.997 (0.995-0.999) | 0.976 (0.721-0.992) | 0.7% (0.0-12.5%) | 0.3% (0.1-0.5%) |

Primary dataset errors: 24 major (S called R), 4 very major (R called S).

Elements most often present in susceptible isolates wrongly called resistant: qnrB19 (13), qnrS1 (6), gyrA_D87Y (2), gyrA_D82N (1), gyrA_S83F (1), gyrA_D87G (1)

Resistant isolates called susceptible: 3 carry no element at all; 1 carry elements, but none linked to this antibiotic. 0 of the 4 have a partial call for an element linked to this antibiotic (the gene may be broken across contigs).

BioProjects with the most errors (overall error rate 0.3%). A rate far above that points to a study-specific problem such as label noise or a different test method:

| bioproject | labelled | errors | error_rate | share_of_all_errors |
|---|---|---|---|---|
| PRJNA1018843 | 10 | 10 | 100.0% | 35.7% |
| PRJNA292661 | 4851 | 5 | 0.1% | 17.9% |
| PRJNA292666 | 817 | 5 | 0.6% | 17.9% |
| PRJNA287785 | 179 | 4 | 2.2% | 14.3% |
| PRJNA324583 | 2 | 2 | 100.0% | 7.1% |

## How to read this

- A strong baseline leaves little room for machine learning. A model must beat these numbers on the same isolates and splits to be worth reporting.
- Wide intervals mean few independent groups, not necessarily a bad result.
- NCBI cautions that carrying a resistance gene does not always mean the isolate is resistant, and some resistance has no known gene.
- The evidence rules come from NCBI's catalog, not from these labels. Two rule choices were refined after reading a first version of this report (family-level symbols; MULTIDRUG-class elements); both alternatives are shown in the sensitivity rows.
