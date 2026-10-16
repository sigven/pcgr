# Gene expression outliers in cancer genes for the PDF report - increased expression of proto-oncogenes, or reduced expression of tumor suppressor genes; percentile as colored pills (colors as in the HTML report), ranked by absolute Z-score

Gene expression outliers in cancer genes for the PDF report - increased
expression of proto-oncogenes, or reduced expression of tumor suppressor
genes; percentile as colored pills (colors as in the HTML report),
ranked by absolute Z-score

## Usage

``` r
pdf_table_expression_outliers(
  outliers,
  direction = "Increased expression",
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- outliers:

  data frame with expression outliers (content - expression

  - outliers)

- direction:

  'Increased expression' or 'Reduced expression'

- color:

  main report color (hex)

## Value

tinytable object, or NULL when there are no outliers
