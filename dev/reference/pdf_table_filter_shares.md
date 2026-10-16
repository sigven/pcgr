# Table of variant classifications (tumor-only germline filtering, exonic status) with counts and shares (bars) for the PDF report

Table of variant classifications (tumor-only germline filtering, exonic
status) with counts and shares (bars) for the PDF report

## Usage

``` r
pdf_table_filter_shares(
  d,
  retained = retained_filter_classes(),
  settings = NULL,
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- d:

  data frame with CLASSIFICATION, N and SHARE (see filtering_summary())

- retained:

  classifications retained (in bold)

- settings:

  tumor-only settings (conf\$somatic_snv\$tumor_only)

- color:

  main report color (hex)

## Value

tinytable object
