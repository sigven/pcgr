# Biomarker evidence table (therapeutic sensitivity/resistance) for the PDF report - one row per matched alteration, with tier, matched therapies (tumor types), and evidence sources (cf. render_actble_bm_table() of the HTML report)

Biomarker evidence table (therapeutic sensitivity/resistance) for the
PDF report - one row per matched alteration, with tier, matched
therapies (tumor types), and evidence sources (cf.
render_actble_bm_table() of the HTML report)

## Usage

``` r
pdf_table_bm_evidence(
  rctbl_recs,
  clnsig = "therapeutic_sensitivity",
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- rctbl_recs:

  list with 'main' data frame, from prep_actble_display_tbl()

- clnsig:

  'therapeutic_sensitivity' or 'therapeutic_resistance'

- color:

  main report color (hex)

## Value

tinytable object, or NULL if there are no evidence records
