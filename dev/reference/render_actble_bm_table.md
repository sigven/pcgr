# Build actionable (predictive) biomarker reactable

Build actionable (predictive) biomarker reactable

## Usage

``` r
render_actble_bm_table(
  rctbl_recs = NULL,
  variant_category = "snv_indel",
  clnsig = "therapeutic_sensitivity",
  color_palette = pcgrr::color_palette
)
```

## Arguments

- rctbl_recs:

  List with \$main and \$nested data frames, as returned by
  `prep_bm_display_tbl`. \$main must contain ACTIONABILITY_TIER with
  values 1 and 2.

- variant_category:

  One of "snv_indel", "cna", "fusion", "msi", "tmb" (msi and tmb for
  predictive biomarkers only)

- clnsig:

  One of "therapeutic_sensitivity" or "therapeutic_resistance",

- color_palette:

  color palette (default: pcgrr::color_palette)

## Value

A reactable object with the biomarker table
