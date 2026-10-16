# Variant allelic fraction as a horizontal bar with the value (Typst markup), as in the HTML report tables

Variant allelic fraction as a horizontal bar with the value (Typst
markup), as in the HTML report tables

## Usage

``` r
pdf_vaf_bar(vaf, width = "1.3cm")
```

## Arguments

- vaf:

  allelic fractions (numeric or character, 0-1)

- width:

  bar width (Typst length)

## Value

character vector with Typst markup ("" for missing values)
