# Deterministic ordering of gene-level records: cancer association ranks (TISSUE_ASSOC_RANK, then GLOBAL_ASSOC_RANK; highest first - as in the variant tables of the report), then optional extra keys, then SYMBOL

Deterministic ordering of gene-level records: cancer association ranks
(TISSUE_ASSOC_RANK, then GLOBAL_ASSOC_RANK; highest first - as in the
variant tables of the report), then optional extra keys, then SYMBOL

## Usage

``` r
rank_order_hl(d, first = NULL, then = NULL)
```

## Arguments

- d:

  data frame with column SYMBOL (rank columns optional)

- first:

  optional leading sort keys (numeric, larger = first)

- then:

  optional sort keys applied after the ranks (larger = first)

## Value

row order (integer)
