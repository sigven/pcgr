# Variant table (SNVs/InDels, or CNAs) for the PDF report - plain text columns, actionability tier and oncogenicity as colored pills

Variant table (SNVs/InDels, or CNAs) for the PDF report - plain text
columns, actionability tier and oncogenicity as colored pills

## Usage

``` r
pdf_table_variants(
  variants,
  color = pcgrr::color_palette$report_color$values[1],
  clnsig = "therapeutic_sensitivity",
  variant_class = TRUE,
  oncokb = FALSE
)
```

## Arguments

- variants:

  data frame with variants

- color:

  main report color (hex)

- clnsig:

  'therapeutic_sensitivity' or 'therapeutic_resistance' (tier colors)

- variant_class:

  show the variant class column

- oncokb:

  show the OncoKB oncogenicity column (ONCOGENICITY_OKB, as colored
  pills), when OncoKB annotations are present

## Value

tinytable object
