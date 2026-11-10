# Table of the reference signatures most similar to the mutational spectrum of the sample (cosine similarity) for the PDF report - signatures previously attributed to the tumor site in the main report color

Table of the reference signatures most similar to the mutational
spectrum of the sample (cosine similarity) for the PDF report -
signatures previously attributed to the tumor site in the main report
color

## Usage

``` r
pdf_table_signature_similarity(
  similarity,
  n = 10,
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- similarity:

  data frame with SIGNATURE_ID, SIMILARITY, AETIOLOGY_KEYWORD and
  SITE_SPECIFIC (content mutational_signatures - result -
  signature_similarity)

- n:

  maximum number of signatures shown

- color:

  main report color (hex)

## Value

tinytable object
