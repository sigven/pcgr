# Add OncoKB evidence items for other tumor types

Combines evidence items from a tumor type specific OncoKB query with
evidence items from a query without tumor type. Items from the tumor
type specific query take precedence: items from the query without tumor
type are only added if no item for the same variant, evidence type,
clinical significance, cancer type and therapy was returned by the tumor
type specific query. Each added item is labelled by its cancer
type/primary site (BM_CANCER_TYPE/BM_PRIMARY_SITE), which is what tier
assignment relies on.

## Usage

``` r
merge_oncokb_evidence_all_tumor_types(
  evidence_df = data.frame(),
  evidence_df_all = data.frame()
)
```

## Arguments

- evidence_df:

  evidence items from tumor type specific query

- evidence_df_all:

  evidence items from query without tumor type

## Value

data.frame with combined evidence items
