# Legend line printed below a PDF report table (raw Typst block): the cancer hotspot symbol ('hotspot' - only when at least one of the variants is located in a hotspot), or the biomarker evidence source logos ('bm_sources')

Legend line printed below a PDF report table (raw Typst block): the
cancer hotspot symbol ('hotspot' - only when at least one of the
variants is located in a hotspot), or the biomarker evidence source
logos ('bm_sources')

## Usage

``` r
pdf_legend(
  type = c("hotspot", "bm_sources"),
  variants = NULL,
  sources = c("civic", "cgi", "oncokb"),
  muted = "pcgr-muted"
)
```

## Arguments

- type:

  'hotspot' or 'bm_sources'

- variants:

  data frame with the variants listed in the table (for 'hotspot')

- sources:

  evidence sources listed (for 'bm_sources'), any of 'civic', 'cgi' and
  'oncokb'

- muted:

  text color of the legend: name of a Typst color variable (e.g.
  'pcgr-muted', as defined in typst/\_preamble.typ) or hex color

## Value

character (raw Typst block; "" when there is nothing to explain)
