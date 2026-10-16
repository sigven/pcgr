# Colored pill (rounded box with bold text) in a PDF report table cell - raw Typst markup; the column must not be escaped by format_tt()

Colored pill (rounded box with bold text) in a PDF report table cell -
raw Typst markup; the column must not be escaped by format_tt()

## Usage

``` r
pdf_pill(label, fill, text_color = "#ffffff")
```

## Arguments

- label:

  pill text(s)

- fill:

  background color(s) (hex)

- text_color:

  text color(s) (hex)

## Value

character vector with Typst markup ("" for empty labels)
