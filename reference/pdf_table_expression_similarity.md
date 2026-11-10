# Most similar external samples (expression correlation) of a reference source for the PDF report - samples from the primary site of the query sample highlighted in the main report color

Most similar external samples (expression correlation) of a reference
source for the PDF report - samples from the primary site of the query
sample highlighted in the main report color

## Usage

``` r
pdf_table_expression_similarity(
  similarity,
  site = "Any",
  n = 10,
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- similarity:

  data frame with similarity results (content - expression -
  similarity_analysis - source)

- site:

  primary tumor site of the sample

- n:

  maximum number of samples shown

- color:

  main report color (hex)

## Value

tinytable object, or NULL when there are no results
