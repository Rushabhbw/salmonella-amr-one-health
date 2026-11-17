# AMRFinderPlus concordance: our run vs NCBI's stored result

Generated: 2026-11-17 14:49  
Script: `workflow/scripts/03_compare_amrfinder.py`

## Our run

```
2026-11-17T05:23:42Z
amrfinder --version:
4.2.7
amrfinder -l:
Running: amrfinder -l
Software directory: /home/rushabh_ubuntu/miniforge3/envs/salmonella-amr-tools/bin/
Software version: 4.2.7
Database directory: /home/rushabh_ubuntu/miniforge3/envs/salmonella-amr-tools/share/amrfinderplus/data/2026-10-29.1
Database version: 2026-10-29.1

datasets --version:
datasets version: 18.41.1
```

## Overall

- Isolates compared: 20
- Identical sets of resistance elements: 20 (100.0%)
- Mean Jaccard similarity: 1.000 (1.000 means identical)
- Isolates with elements only in NCBI's result: 0
- Isolates with elements only in our result: 0

## By NCBI reference-gene-catalog version

| ncbi_catalog | same_db_as_ours | isolates | identical |
|---|---|---|---|
| 2025-12-03.1 | no | 17 | 17 |
| 2026-05-15.1 | no | 3 | 3 |

## Differences

None. Every compared isolate has the same set of elements.

## How to read this

- A difference can come from a different database version, a renamed element, or a PARTIAL call that one run reports and the other does not. Check the call or method column above before assuming an error.
- Only element names are compared, not copy number, position or call quality.
- This is a sample, so it supports reproducibility checking only. It does not measure prediction accuracy.
