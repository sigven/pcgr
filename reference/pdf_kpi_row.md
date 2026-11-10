# Row of KPI boxes (pcgr-kpi-box(), defined in typst/\_preamble.typ) in a PDF report section - as the bslib value boxes of the HTML report

Row of KPI boxes (pcgr-kpi-box(), defined in typst/\_preamble.typ) in a
PDF report section - as the bslib value boxes of the HTML report

## Usage

``` r
pdf_kpi_row(boxes)
```

## Arguments

- boxes:

  list of KPI boxes, each a list with 'title', 'value', and optionally
  'sub' and 'color' (hex; default: color_palette\$bg_dark)

## Value

raw Typst block (character), to be printed with cat() in a chunk with
output 'asis'
