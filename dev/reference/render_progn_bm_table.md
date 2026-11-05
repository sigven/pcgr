# Build prognostic biomarker reactable

Build prognostic biomarker reactable

## Usage

``` r
render_progn_bm_table(
  rctbl_recs = NULL,
  variant_category = "snv_indel",
  clnsig = "prognostic_poor",
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

  One of "prognostic_poor" or "prognostic_better"

- color_palette:

  color palette (default: pcgrr::color_palette)

## Value

A reactable object with the biomarker table
