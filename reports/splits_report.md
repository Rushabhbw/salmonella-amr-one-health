# Evaluation splits report

Generated: 2026-11-24 10:13  
Script: `workflow/scripts/07_make_splits.py`  
Seed 42; folds 5; time cutoff 2016. Only the primary dataset is split.

## Sizes

| Split | Train | Test | Excluded |
|---|---|---|---|
| transfer_naive | 8972 | 270 | 0 |
| transfer_disjoint | 5158 | 270 | 3814 |
| time | 5489 | 3715 | 0 |
| fold5 and random5 (each fold is the test set once) | see folds table | 1848-1849 per fold | 0 |

## Leakage check: test isolates whose SNP-cluster group also appears in training

| Split | Test isolates | With a relative in train | Share |
|---|---|---|---|
| transfer_naive | 270 | 162 | 60.0% |
| transfer_disjoint | 270 | 0 | 0.0% |
| time | 3715 | 3253 | 87.6% |
| fold5 (all folds) | 9242 | 0 | 0.0% |
| random5 (all folds) | 9242 | 8539 | 92.4% |

`fold5` and `transfer_disjoint` must show 0. `transfer_naive` and `random5` show how leaky the simple splits would be. `time` is expected to share lineages across years; that is a realistic forward-prediction setting.

## Resistant / susceptible isolates in train and test

| Antibiotic | Split | Train R / S | Test R / S |
|---|---|---|---|
| tetracycline | transfer_disjoint | 2706 / 2427 | 145 / 113 |
| tetracycline | time | 2837 / 2618 | 1863 / 1842 |
| ampicillin | transfer_disjoint | 1176 / 3982 | 129 / 141 |
| ampicillin | time | 1628 / 3861 | 714 / 3001 |
| ceftriaxone | transfer_disjoint | 639 / 4434 | 45 / 224 |
| ceftriaxone | time | 765 / 4670 | 402 / 3276 |
| ciprofloxacin | transfer_disjoint | 32 / 5072 | 11 / 204 |
| ciprofloxacin | time | 17 / 5357 | 584 / 2993 |

## Folds (grouped 5-fold)

| Fold | Isolates | Groups | R tetracycline | R ampicillin | R ceftriaxone | R ciprofloxacin |
|---|---|---|---|---|---|---|
| 0 | 1849 | 224 | 1143 | 532 | 386 | 558 |
| 1 | 1849 | 225 | 1117 | 608 | 273 | 20 |
| 2 | 1848 | 224 | 570 | 598 | 271 | 14 |
| 3 | 1848 | 224 | 687 | 382 | 120 | 6 |
| 4 | 1848 | 224 | 1201 | 243 | 127 | 5 |

Every fold has at least 10 resistant isolates for each core antibiotic.

## Why the clean transfer test matters

- Groups containing both clinical and non-clinical isolates: 64
- Clinical isolates in such shared groups: 162 of 270
- Non-clinical isolates dropped from training in the clean version: 3814

## Concentration of the clinical test set

- Clinical isolates by collection year (top 3): 2014: 204, 2019: 19, 2016: 2 (of 232 with a year)

- Clinical isolates by BioProject (top 3): PRJNA287785: 222, PRJNA292661: 18, PRJNA562719: 9 (of 270)

- If one year or study dominates, the transfer test mostly measures transfer to that collection, not to human isolates in general.
