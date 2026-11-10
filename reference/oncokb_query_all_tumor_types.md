# Whether OncoKB should also be queried without a tumor type

When OncoKB is queried with a tumor type, it returns the treatments for
that tumor type, plus only a subset of the evidence for other tumor
types (re-leveled according to OncoKB's level propagation rules). A
second query without tumor type returns the evidence for all tumor
types, at their original levels. The second query is only needed when
the first query was restricted to a tumor type.

## Usage

``` r
oncokb_query_all_tumor_types(all_tumor_types = FALSE, oncotree_code = NULL)
```

## Arguments

- all_tumor_types:

  logical/integer, retrieve evidence for all tumor types

- oncotree_code:

  OncoTree code used for the tumor type specific query

## Value

logical
