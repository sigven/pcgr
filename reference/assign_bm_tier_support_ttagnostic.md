# For all biomarker evidence items with specified confidence level, assign whether these are tier-defining or providing additional support - tumor-type agnostic query

The function considers the strength of evidence (evidence levels) and
the match between biomarker site and primary site of query tumor. Mark
evidence considered tier-defining, while other biomarker evidence items
as additional support

## Usage

``` r
assign_bm_tier_support_ttagnostic(
  variants_tier_classified = NULL,
  vartype = NULL,
  etype_for_tiering = c("predictive"),
  biomarker_items = NULL
)
```

## Arguments

- variants_tier_classified:

  data frame with variant tier classifications based on support from
  biomarker evidence item

- vartype:

  variant type (e.g. 'snv_indel', 'cna', 'fusion')

- etype_for_tiering:

  evidence type(s) to consider for tier classification

- biomarker_items:

  data frame with biomarker evidence items

## Value

data frame with biomarker evidence items classified as tier-defining or
providing additional support
