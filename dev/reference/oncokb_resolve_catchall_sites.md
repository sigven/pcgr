# Primary site of OncoKB catch-all cancer types, relative to the query tumor

OncoKB's catch-all cancer types (e.g. 'All Solid Tumors', 'All Liquid
Tumors', 'All Tumors') have no tissue, and their evidence items are
extracted with BM_PRIMARY_SITE = 'Any' (pan-cancer). Catch-alls
restricted to liquid tumors are relabelled 'Liquid tumors' for a solid
query tumor, and catch-alls restricted to solid tumors 'Solid tumors'
for a liquid query tumor (Lymphoid/Myeloid), so that tier assignment
treats them as evidence from other tumor sites rather than as pan-cancer
evidence.

## Usage

``` r
oncokb_resolve_catchall_sites(biomarker_items = NULL, primary_site = "Any")
```

## Arguments

- biomarker_items:

  data.frame with biomarker evidence items

- primary_site:

  primary site of query tumor

## Value

data.frame with biomarker evidence items
