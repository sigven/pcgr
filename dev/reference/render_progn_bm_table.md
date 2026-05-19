# Build biomarker reactable with category-aware styling Combines tier 1 and tier 2 records in one table. Header uses tier 1 color; THERAPY_MATCH cell background reflects the row's tier (1 or 2).

Build biomarker reactable with category-aware styling Combines tier 1
and tier 2 records in one table. Header uses tier 1 color; THERAPY_MATCH
cell background reflects the row's tier (1 or 2).

## Usage

``` r
render_progn_bm_table(
  rctbl_recs = NULL,
  variant_category = "snv_indel",
  clnsig = "prognostic_poor"
)
```

## Arguments

- rctbl_recs:

  List with \$main and \$nested data frames. \$main must contain
  ACTIONABILITY_TIER with values 1 and 2.

- variant_category:

  One of "snv_indel", "cnv", "fusion"

- clnsig:

  One of "therapeutic_sensitivity" or "therapeutic_resistance",

## Value

A reactable object with the biomarker table
