# AMRFinderPlus concordance: our run vs NCBI's stored result

Generated: 2026-11-17 16:01  
Script: `workflow/scripts/03_compare_amrfinder.py`

## Our run

```
2026-11-17T09:25:47Z
amrfinder --version:
4.2.7
amrfinder -l:
Running: amrfinder -l
Software version: 4.2.7
Database version: 2026-10-29.1

datasets --version:
datasets version: 18.41.1
```

## Overall

- Isolates compared: 300
- Identical sets of resistance elements: 295 (98.3%)
- Mean Jaccard similarity: 0.995 (1.000 means identical)
- Isolates with elements only in NCBI's result: 3
- Isolates with elements only in our result: 2

## By NCBI reference-gene-catalog version

| ncbi_catalog | same_db_as_ours | isolates | identical |
|---|---|---|---|
| 2025-12-03.1 | no | 284 | 279 |
| 2026-05-15.1 | no | 10 | 10 |
| 2026-08-07.1 | no | 6 | 6 |

## Differences

| isolate | host_group | ncbi_catalog | only_ncbi | only_local |
|---|---|---|---|---|
| PDT000161187.2 | clinical | 2025-12-03.1 | dfrA51 |  |
| PDT000161200.2 | clinical | 2025-12-03.1 | dfrA51 |  |
| PDT000121333.2 | non_clinical | 2025-12-03.1 | gar |  |
| PDT000598974.1 | clinical | 2025-12-03.1 |  | ramR_Q19Ter |
| PDT000598976.1 | clinical | 2025-12-03.1 |  | mrx(A) |

Differing elements (what the call or method was on the side that reported it):

| side | element | call_or_method | isolates |
|---|---|---|---|
| only in NCBI result | dfrA51 | COMPLETE | 2 |
| only in NCBI result | gar | HMM | 1 |
| only in our result | mrx(A) | EXACTX | 1 |
| only in our result | ramR_Q19Ter | POINTX | 1 |

## How to read this

- A difference can come from a different database version, a renamed element, or a PARTIAL call that one run reports and the other does not. Check the call or method column above before assuming an error.
- Only element names are compared, not copy number, position or call quality.
- This is a sample, so it supports reproducibility checking only. It does not measure prediction accuracy.
