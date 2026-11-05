# Function that gathers data tables on biomarker variants (predictive, prognostic or diagnostic) for display in report

Shared implementation behind `prep_actble_display_tbl`,
`prep_progn_display_tbl` and `prep_diagn_display_tbl`.

## Usage

``` r
prep_bm_display_tbl(
  rep = NULL,
  tier = c(1, 2),
  etype_for_tiering = c("predictive"),
  clnsig = "therapeutic_sensitivity",
  tier_defining_eitems_only = TRUE,
  variant_category = "snv_indel",
  evidence_type = "predictive"
)
```

## Arguments

- rep:

  report object

- tier:

  tier level(s) to consider for display (e.g. 1, 2)

- etype_for_tiering:

  evidence type(s) to consider for tiering (e.g. 'predictive',
  'prognostic', 'diagnostic')

- clnsig:

  clinical significance to consider for tiering (e.g.
  'therapeutic_sensitivity', 'prognostic_poor', 'diagnostic_positive')

- tier_defining_eitems_only:

  consider only evidence items that were used for tiering (e.g. for tier
  1: only evidence items with A-level evidence, for tier 2: only
  evidence items with B-level evidence). If FALSE, all evidence items
  associated with each variant, not only the tier-defining evidence
  items, will be considered for display in the report.

- variant_category:

  cna, snv_indel, fusion, msi, or tmb (msi and tmb for predictive
  biomarkers only)

- evidence_type:

  type of biomarker table to prepare, one of 'predictive', 'prognostic'
  or 'diagnostic'

## Value

list with 'main' (one row per variant) and 'nested' (one row per
evidence item) data frames
