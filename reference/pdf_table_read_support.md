# Table of the depth/allelic support filtering steps of SNVs/InDels for the PDF report (variant filtering section)

Table of the depth/allelic support filtering steps of SNVs/InDels for
the PDF report (variant filtering section)

## Usage

``` r
pdf_table_read_support(
  steps,
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- steps:

  data frame with STEP, THRESHOLD, APPLIED, NOTE, N_REMOVED and
  N_REMAINING (see filter_read_support())

- color:

  main report color (hex)

## Value

tinytable object
