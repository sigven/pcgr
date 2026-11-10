# Build biomarker reactable with category-aware styling Combines tier 1 and tier 2 (and for predictive biomarkers, tier 3) records in one table. Header uses tier 1 color; the summary cell (THERAPY_MATCH, PROGNOSTIC_OUTCOME or DIAGNOSTIC_EVIDENCE) background reflects the row's tier.

Shared implementation behind `render_actble_bm_table`,
`render_progn_bm_table` and `render_diagn_bm_table`.

## Usage

``` r
render_bm_table(
  rctbl_recs = NULL,
  variant_category = "snv_indel",
  clnsig = "therapeutic_sensitivity",
  color_palette = pcgrr::color_palette,
  evidence_type = "predictive"
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

  One of "therapeutic_sensitivity", "therapeutic_resistance",
  "prognostic_poor", "prognostic_better", "diagnostic_positive" or
  "diagnostic_negative"

- color_palette:

  color palette (default: pcgrr::color_palette)

- evidence_type:

  One of "predictive", "prognostic" or "diagnostic"

## Value

A reactable object with the biomarker table
