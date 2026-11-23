# Element-to-class mapping report

Generated: 2026-11-23 09:53  
Script: `workflow/scripts/05_map_elements_to_classes.py`  
Catalog file: `ReferenceGeneCatalog.txt`

Catalog database version(s) seen for matched elements: 2026-10-29.1

## Rules

- tetracycline: class or subclass contains TETRACYCLINE
- ampicillin: class contains BETA-LACTAM (any beta-lactamase)
- ceftriaxone: class contains BETA-LACTAM and subclass contains CEPHALOSPORIN (NCBI's definition of CEPHALOSPORIN is the Lahey 2be group; CARBAPENEM subclass is not counted)
- ciprofloxacin: class or subclass contains QUINOLONE
- An element counts only if every catalog entry it matches satisfies the rule (a family-level symbol such as blaTEM is not evidence for ceftriaxone, because TEM-1 does not break cephalosporins).
- MULTIDRUG-class elements are not drug-specific evidence; they are listed separately and used only in a sensitivity analysis.
- NCBI cautions that carrying a resistance gene does not always mean the isolate is resistant.

## Matching

- Elements: 158
- Matched on allele: 58
- Matched on gene family: 87
- Not found in the catalog: 13
- With more than one class/subclass pair (ambiguous): 12

Not found (check these by hand):

| element | kind | n_isolates_primary |
|---|---|---|
| dfrA51 | gene | 63 |
| ramR_E146SfsTer42 | point_mutation | 0 |
| ramR_A171IfsTer12 | point_mutation | 0 |
| ramR_N36RfsTer11 | point_mutation | 0 |
| aac(3)-C322 | gene | 0 |
| ramR_D169GfsTer148 | point_mutation | 0 |
| ramR_E175IfsTer136 | point_mutation | 0 |
| ramR_L66PfsTer78 | point_mutation | 0 |
| ramR_L74RfsTer9 | point_mutation | 0 |
| ramR_Q23Ter | point_mutation | 0 |
| aph(2'')-I | weak_only | weak calls only |
| aph(2'')-Ia | weak_only | weak calls only |
| aph(3'') | weak_only | weak calls only |

## Class and subclass of the matched elements

| class | subclass | elements |
|---|---|---|
| QUINOLONE | QUINOLONE | 17 |
| AMINOGLYCOSIDE | STREPTOMYCIN | 14 |
| TRIMETHOPRIM | TRIMETHOPRIM | 11 |
| FOSFOMYCIN | FOSFOMYCIN | 11 |
| BETA-LACTAM | CEPHALOSPORIN | 10 |
| BETA-LACTAM | BETA-LACTAM | 10 |
| AMINOGLYCOSIDE | GENTAMICIN | 8 |
| PHENICOL | CHLORAMPHENICOL | 7 |
| TETRACYCLINE | TETRACYCLINE | 6 |
| RIFAMYCIN | RIFAMYCIN | 3 |
| SULFONAMIDE | SULFONAMIDE | 3 |
| MACROLIDE | ERYTHROMYCIN | 3 |
| NITROFURAN/PHENICOL/QUINOLONE/TETRACYCLINE|PHENICOL/QUINOLONE | NITROFURANTOIN/PHENICOL/QUINOLONE/TIGECYCLINE|PHENICOL/QUINOLONE | 2 |
| MACROLIDE | TYLOSIN | 2 |
| PHENICOL | CHLORAMPHENICOL/FLORFENICOL | 2 |
| LINCOSAMIDE | LINCOSAMIDE | 2 |
| BLEOMYCIN | BLEOMYCIN | 2 |
| BETA-LACTAM | CARBAPENEM | 2 |
| MACROLIDE | AZITHROMYCIN/ERYTHROMYCIN/SPIRAMYCIN/TELITHROMYCIN | 2 |
| AMINOGLYCOSIDE | KANAMYCIN | 2 |
| AMINOGLYCOSIDE | GENTAMICIN/KANAMYCIN/TOBRAMYCIN | 2 |
| QUINOLONE|QUINOLONE/TRICLOSAN | QUINOLONE|QUINOLONE/TRICLOSAN | 2 |
| MACROLIDE | AZITHROMYCIN|AZITHROMYCIN/CLARITHROMYCIN/ERYTHROMYCIN | 2 |
| AMINOGLYCOSIDE | APRAMYCIN/GENTAMICIN/TOBRAMYCIN | 1 |
| AMINOGLYCOSIDE | AMIKACIN/GENTAMICIN/KANAMYCIN/TOBRAMYCIN|AMIKACIN/KANAMYCIN/TOBRAMYCIN|AMINOGLYCOSIDE|GENTAMICIN|KANAMYCIN | 1 |

## Elements that count as evidence for each antibiotic

### tetracycline (6 elements)

| element | kind | class | subclass | n_isolates_primary | n_groups_primary |
|---|---|---|---|---|---|
| tet(A) | gene | TETRACYCLINE | TETRACYCLINE | 2740 | 299 |
| tet(B) | gene | TETRACYCLINE | TETRACYCLINE | 1817 | 183 |
| tet(C) | gene | TETRACYCLINE | TETRACYCLINE | 208 | 52 |
| tet(G) | gene | TETRACYCLINE | TETRACYCLINE | 50 | 14 |
| tet(M) | gene | TETRACYCLINE | TETRACYCLINE | 16 | 13 |
| tet(D) | gene | TETRACYCLINE | TETRACYCLINE | 7 | 7 |

Excluded as MULTIDRUG-class or not drug-specific (sensitivity analysis only): ramR_T18P

### ampicillin (27 elements)

| element | kind | class | subclass | n_isolates_primary | n_groups_primary |
|---|---|---|---|---|---|
| blaTEM-1 | gene | BETA-LACTAM | BETA-LACTAM | 1131 | 141 |
| blaCMY-2 | gene | BETA-LACTAM | CEPHALOSPORIN | 806 | 168 |
| blaCTX-M-65 | gene | BETA-LACTAM | CEPHALOSPORIN | 317 | 6 |
| blaHER-3 | gene | BETA-LACTAM | BETA-LACTAM | 90 | 20 |
| blaTEM | gene | BETA-LACTAM | AMOXICILLIN-CLAVULANIC ACID/PIPERACILLIN-TAZOBACTAM/TICARCILLIN-CLAVULANIC ACID|BETA-LACTAM|CEPHALOSPORIN | 74 | 17 |
| blaCARB-2 | gene | BETA-LACTAM | BETA-LACTAM | 59 | 19 |
| blaCMY | gene | BETA-LACTAM | AVIBACTAM/CEPHALOSPORIN|CEPHALOSPORIN|CEPHALOSPORIN/TANIBORBACTAM | 41 | 21 |
| blaCTX-M | gene | BETA-LACTAM | BETA-LACTAM|CEPHALOSPORIN | 10 | 1 |
| blaCTX-M-15 | gene | BETA-LACTAM | CEPHALOSPORIN | 6 | 3 |
| blaCTX-M-1 | gene | BETA-LACTAM | CEPHALOSPORIN | 6 | 5 |
| blaHER | gene | BETA-LACTAM | BETA-LACTAM | 5 | 4 |
| blaSHV-12 | gene | BETA-LACTAM | CEFIDEROCOL/CEPHALOSPORIN | 5 | 5 |
| blaTEM-135 | gene | BETA-LACTAM | BETA-LACTAM | 3 | 3 |
| blaCTX-M-14 | gene | BETA-LACTAM | CEPHALOSPORIN | 1 | 1 |
| blaLAP-2 | gene | BETA-LACTAM | BETA-LACTAM | 1 | 1 |
| blaOXA-1 | gene | BETA-LACTAM | CEPHALOSPORIN | 1 | 1 |
| blaPSE | gene | BETA-LACTAM | BETA-LACTAM | 1 | 1 |
| blaCMY-4 | gene | BETA-LACTAM | CEPHALOSPORIN | 1 | 1 |
| blaKPC-2 | gene | BETA-LACTAM | CARBAPENEM | 1 | 1 |
| blaNDM-1 | gene | BETA-LACTAM | CARBAPENEM | 1 | 1 |
| blaOXA-9 | gene | BETA-LACTAM | BETA-LACTAM | 1 | 1 |
| blaCTX-M-55 | gene | BETA-LACTAM | CEPHALOSPORIN | 0 | 0 |
| blaLAP | gene | BETA-LACTAM | BETA-LACTAM | 0 | 0 |
| blaSHV-2A | gene | BETA-LACTAM | CEPHALOSPORIN | 0 | 0 |
| bla | weak_only | BETA-LACTAM | BETA-LACTAM | weak calls only |  |
| blaDHA | weak_only | BETA-LACTAM | CEPHALOSPORIN | weak calls only |  |
| blaOXA | weak_only | BETA-LACTAM | BETA-LACTAM|CARBAPENEM|CEFIDEROCOL/CEPHALOSPORIN|CEPHALOSPORIN | weak calls only |  |

Excluded as MULTIDRUG-class or not drug-specific (sensitivity analysis only): ramR_T18P

### ceftriaxone (12 elements)

| element | kind | class | subclass | n_isolates_primary | n_groups_primary |
|---|---|---|---|---|---|
| blaCMY-2 | gene | BETA-LACTAM | CEPHALOSPORIN | 806 | 168 |
| blaCTX-M-65 | gene | BETA-LACTAM | CEPHALOSPORIN | 317 | 6 |
| blaCMY | gene | BETA-LACTAM | AVIBACTAM/CEPHALOSPORIN|CEPHALOSPORIN|CEPHALOSPORIN/TANIBORBACTAM | 41 | 21 |
| blaCTX-M-15 | gene | BETA-LACTAM | CEPHALOSPORIN | 6 | 3 |
| blaCTX-M-1 | gene | BETA-LACTAM | CEPHALOSPORIN | 6 | 5 |
| blaSHV-12 | gene | BETA-LACTAM | CEFIDEROCOL/CEPHALOSPORIN | 5 | 5 |
| blaCTX-M-14 | gene | BETA-LACTAM | CEPHALOSPORIN | 1 | 1 |
| blaOXA-1 | gene | BETA-LACTAM | CEPHALOSPORIN | 1 | 1 |
| blaCMY-4 | gene | BETA-LACTAM | CEPHALOSPORIN | 1 | 1 |
| blaCTX-M-55 | gene | BETA-LACTAM | CEPHALOSPORIN | 0 | 0 |
| blaSHV-2A | gene | BETA-LACTAM | CEPHALOSPORIN | 0 | 0 |
| blaDHA | weak_only | BETA-LACTAM | CEPHALOSPORIN | weak calls only |  |

### ciprofloxacin (22 elements)

| element | kind | class | subclass | n_isolates_primary | n_groups_primary |
|---|---|---|---|---|---|
| gyrA_D87Y | point_mutation | QUINOLONE | QUINOLONE | 722 | 19 |
| qnrB19 | gene | QUINOLONE | QUINOLONE | 51 | 28 |
| gyrA_S83F | point_mutation | QUINOLONE|QUINOLONE/TRICLOSAN | QUINOLONE|QUINOLONE/TRICLOSAN | 27 | 10 |
| gyrA_D87N | point_mutation | QUINOLONE | QUINOLONE | 20 | 6 |
| qnrS1 | gene | QUINOLONE | QUINOLONE | 14 | 6 |
| gyrA_D87G | point_mutation | QUINOLONE|QUINOLONE/TRICLOSAN | QUINOLONE|QUINOLONE/TRICLOSAN | 12 | 9 |
| gyrA_S83Y | point_mutation | QUINOLONE | QUINOLONE | 9 | 7 |
| parC_S80I | point_mutation | QUINOLONE | QUINOLONE | 4 | 3 |
| qnrA1 | gene | QUINOLONE | QUINOLONE | 3 | 3 |
| qnrS2 | gene | QUINOLONE | QUINOLONE | 3 | 2 |
| gyrA_D82N | point_mutation | QUINOLONE | QUINOLONE | 2 | 2 |
| qnrB2 | gene | QUINOLONE | QUINOLONE | 2 | 2 |
| aac(6')-Ib-cr5 | gene | AMINOGLYCOSIDE/QUINOLONE | AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN | 1 | 1 |
| gyrB_S464Y | point_mutation | QUINOLONE | QUINOLONE | 1 | 1 |
| oqxA | gene | NITROFURAN/PHENICOL/QUINOLONE/TETRACYCLINE|PHENICOL/QUINOLONE | NITROFURANTOIN/PHENICOL/QUINOLONE/TIGECYCLINE|PHENICOL/QUINOLONE | 1 | 1 |
| oqxB | gene | NITROFURAN/PHENICOL/QUINOLONE/TETRACYCLINE|PHENICOL/QUINOLONE | NITROFURANTOIN/PHENICOL/QUINOLONE/TIGECYCLINE|PHENICOL/QUINOLONE | 1 | 1 |
| qnrB1 | gene | QUINOLONE | QUINOLONE | 1 | 1 |
| qnrS4 | gene | QUINOLONE | QUINOLONE | 1 | 1 |
| qnrB77 | gene | QUINOLONE | QUINOLONE | 1 | 1 |
| qnrS13 | gene | QUINOLONE | QUINOLONE | 0 | 0 |
| parE_A499T | point_mutation | QUINOLONE | QUINOLONE | 0 | 0 |
| qnrB | weak_only | QUINOLONE | QUINOLONE | weak calls only |  |

Excluded as MULTIDRUG-class or not drug-specific (sensitivity analysis only): ramR_T18P
