# Evidence items from sources that are not used for tier assignment

Evidence items from sources not used for tier assignment (e.g. CGI when
OncoKB is run) are kept as additional support for variants with a tier
(assigned from the other sources), and dropped for variants without a
tier.

## Usage

``` r
set_bm_support_sources_not_tiered(
  biomarker_items = NULL,
  sources_not_tiered = NULL
)
```

## Arguments

- biomarker_items:

  data frame with biomarker evidence items, with ACTIONABILITY_TIER and
  BM_ACTIONABILITY_SUPPORT assigned

- sources_not_tiered:

  biomarker sources (BM_SOURCE_DB) not used for tier assignment

## Value

data frame with biomarker evidence items
