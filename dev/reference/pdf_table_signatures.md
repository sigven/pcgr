# Table of mutational signature contributions (signature refitting) for the PDF report - signature (colored as in the contribution plot), contribution with bootstrap confidence interval, and aetiology

Table of mutational signature contributions (signature refitting) for
the PDF report - signature (colored as in the contribution plot),
contribution with bootstrap confidence interval, and aetiology

## Usage

``` r
pdf_table_signatures(
  per_signature,
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- per_signature:

  data frame with signature contributions (content
  mutational_signatures - result - contributions - per_signature)

- color:

  main report color (hex)

## Value

tinytable object
