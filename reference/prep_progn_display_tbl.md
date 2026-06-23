# Function that gathers data tables on prognostic variants for display in report

Function that gathers data tables on prognostic variants for display in
report

## Usage

``` r
prep_progn_display_tbl(
  rep = NULL,
  tier = c(1, 2),
  etype_for_tiering = c("prognostic"),
  clnsig = "prognostic_better",
  tier_defining_eitems_only = TRUE,
  variant_category = "snv_indel"
)
```

## Arguments

- rep:

  report object

- tier:

  tier level(s) to consider for display (e.g. 1, 2)

- etype_for_tiering:

  evidence type(s) to consider for tiering (e.g. 'prognostic')

- clnsig:

  clinical significance to consider for tiering (e.g.
  'prognostic_better', 'prognostic_poor')

- tier_defining_eitems_only:

  consider only evidence items that were used for tiering (e.g. for tier
  1: only evidence items with A-level evidence, for tier 2: only
  evidence items with B-level evidence). If FALSE, all evidence items
  associated with each variant, not only the tier-defining evidence
  items, will be considered for display in the report.

- variant_category:

  cna, snv_indel, or fusion
