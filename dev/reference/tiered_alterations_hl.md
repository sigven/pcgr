# Somatic alterations with an actionability tier (SNVs/InDels, CNAs, RNA fusions) as display labels

Somatic alterations with an actionability tier (SNVs/InDels, CNAs, RNA
fusions) as display labels

## Usage

``` r
tiered_alterations_hl(
  report,
  tier = c(1, 2),
  clnsig = "therapeutic_sensitivity"
)
```

## Arguments

- report:

  PCGR report object

- tier:

  actionability tiers

- clnsig:

  type of biomarker evidence the tiers refer to -
  'therapeutic_sensitivity' (default, the tier of each variant record),
  or e.g. 'prognostic_poor', 'prognostic_better', 'diagnostic_positive'
  (tiers from the corresponding bm_evidence classification of each
  callset)

## Value

data frame with LABEL, SYMBOL, ALTERATION and ACTIONABILITY_TIER
