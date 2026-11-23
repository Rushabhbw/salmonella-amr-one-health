# Feature table report

Generated: 2026-11-23 09:53  
Script: `workflow/scripts/04_build_features.py`

## Rules

- Feature = a resistance gene or point mutation from NCBI's stored AMRFinderPlus core results.
- 1 = present with a COMPLETE or POINT call. Weak calls (PARTIAL, PARTIAL_END_OF_CONTIG, MISTRANSLATION, HMM) are kept in a separate table and are not used as features.
- Copy number is ignored. No feature is removed here; rare-feature filtering must be done inside cross-validation.

## Size

- Isolates: 9783
- Elements with at least one COMPLETE/POINT call: 149
- Elements seen only with weak calls (never COMPLETE/POINT; not features): 9
- Isolates with no element at all: 2936 (30.0%)
- Elements per isolate (median / max): 2 / 21

## Elements by how widely they occur (primary dataset)

| Elements present in | Number of elements |
|---|---|
| >= 1 independent groups | 127 |
| >= 3 independent groups | 77 |
| >= 10 independent groups | 38 |
| >= 20 independent groups | 25 |
| >= 10 isolates | 54 |
| >= 30 isolates | 40 |
| >= 100 isolates | 28 |

## Most common elements (primary dataset)

| element | kind | n_isolates_primary | n_groups_primary | n_clinical | n_non_clinical |
|---|---|---|---|---|---|
| tet(A) | gene | 2740 | 299 | 87 | 2788 |
| aph(3'')-Ib | gene | 2538 | 239 | 95 | 2588 |
| aph(6)-Id | gene | 2464 | 228 | 94 | 2501 |
| tet(B) | gene | 1817 | 183 | 52 | 1849 |
| sul2 | gene | 1472 | 196 | 88 | 1421 |
| sul1 | gene | 1363 | 195 | 57 | 1389 |
| aadA1 | gene | 1251 | 148 | 24 | 1323 |
| blaTEM-1 | gene | 1131 | 141 | 77 | 1114 |
| blaCMY-2 | gene | 806 | 168 | 31 | 810 |
| gyrA_D87Y | point_mutation | 722 | 19 | 13 | 745 |
| fosA7 | gene | 711 | 22 | 20 | 750 |
| aph(3')-Ia | gene | 623 | 87 | 11 | 650 |
| floR | gene | 549 | 68 | 50 | 555 |
| aac(3)-VIa | gene | 476 | 99 | 8 | 500 |
| aph(4)-Ia | gene | 472 | 32 | 5 | 509 |

## Sanity check: isolates with no element

If resistance is explained by known genes, isolates with no element should be mostly susceptible. Resistant isolates here have no known explanation in the data (unknown mechanism, label noise, or a gap in the database).

| Antibiotic | Isolates with | R | S | R share |
|---|---|---|---|---|
| tetracycline | no element | 36 | 2874 | 1.2% |
| tetracycline | at least one element | 4908 | 1901 | 72.1% |
| ampicillin | no element | 16 | 2920 | 0.5% |
| ampicillin | at least one element | 2478 | 4369 | 36.2% |
| ceftriaxone | no element | 9 | 2826 | 0.3% |
| ceftriaxone | at least one element | 1211 | 5496 | 18.1% |
| ciprofloxacin | no element | 18 | 2904 | 0.6% |
| ciprofloxacin | at least one element | 717 | 5867 | 10.9% |
