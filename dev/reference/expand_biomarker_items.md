# Function that expands biomarker evidence items with variant annotations

Function that expands biomarker evidence items with variant annotations

## Usage

``` r
expand_biomarker_items(
  callset = NULL,
  variant_origin = "somatic",
  target_genes = NULL
)
```

## Arguments

- callset:

  list object with 'variant' and 'biomarker_evidence' data frames

- variant_origin:

  'somatic' or 'germline'

- target_genes:

  data frame with target genes of interest
