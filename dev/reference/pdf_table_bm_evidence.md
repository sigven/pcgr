# Biomarker evidence table (therapeutic sensitivity/resistance, prognostic or diagnostic) for the PDF report - one row per matched alteration, with tier, matched therapies / outcomes / diagnostic evidence (tumor types), and evidence sources (cf. render_actble_bm_table() of the HTML report)

Biomarker evidence table (therapeutic sensitivity/resistance, prognostic
or diagnostic) for the PDF report - one row per matched alteration, with
tier, matched therapies / outcomes / diagnostic evidence (tumor types),
and evidence sources (cf. render_actble_bm_table() of the HTML report)

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

  list with 'main' data frame, from prep_actble_display_tbl(),
  prep_progn_display_tbl() or prep_diagn_display_tbl()

- clnsig:

  evidence type - 'therapeutic_sensitivity', 'therapeutic_resistance',
  'prognostic_poor', 'prognostic_better', 'diagnostic_positive' or
  'diagnostic_negative' (tier colors)

- color:

  main report color (hex)

## Value

tinytable object, or NULL if there are no evidence records
