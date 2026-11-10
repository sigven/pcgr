# RNA fusion table for the PDF report - fusion events ranked with rank_fusions(), with recurrence in the tumor site of the sample and overall recurrence (Mitelman DB), and actionability tier as colored pills

RNA fusion table for the PDF report - fusion events ranked with
rank_fusions(), with recurrence in the tumor site of the sample and
overall recurrence (Mitelman DB), and actionability tier as colored
pills

## Usage

``` r
pdf_table_fusions(
  fusions,
  site = "Any",
  max_rows = 25,
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- fusions:

  data frame with fusions (fusion callset - variant_display)

- site:

  primary tumor site of the sample

- max_rows:

  maximum number of fusions shown (after ranking)

- color:

  main report color (hex)

## Value

tinytable object
