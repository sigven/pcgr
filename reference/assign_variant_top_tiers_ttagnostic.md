# Assign top tiers of clinical significance (AMP/ASCO/CAP framework) to variants based on data from biomarker evidence items (tumor-type agnostic query)

The function considers the strength of evidence (evidence levels) and
the match between biomarker site and primary site of query tumor.

## Usage

``` r
assign_variant_top_tiers_ttagnostic(
  biomarker_items = NULL,
  biomarker_mapping_confidence = "high",
  var_df = NULL,
  etype_for_tiering = c("predictive")
)
```

## Arguments

- biomarker_items:

  data frame with biomarker evidence items

- biomarker_mapping_confidence:

  confidence level of variant-biomarker mapping resolution (e.g. 'high'
  or 'medium')

- var_df:

  data frame with variants (SNVs/InDels, CNAs, fusions)

- etype_for_tiering:

  evidence type(s) to consider for tier classification (e.g.
  'predictive', 'prognostic', 'diagnostic')

## Value

data frame with top tier classifications for variants based on biomarker
evidence items
