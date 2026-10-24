# Prognostic and diagnostic biomarker evidence of Tier 1/2 aberrations (SNVs/InDels, CNAs, RNA fusions) for the PDF report - one table per variant category and evidence type (poor/better outcome, diagnosis supported), same evidence selection as the prognostic and diagnostic tabs of the HTML report

Prognostic and diagnostic biomarker evidence of Tier 1/2 aberrations
(SNVs/InDels, CNAs, RNA fusions) for the PDF report - one table per
variant category and evidence type (poor/better outcome, diagnosis
supported), same evidence selection as the prognostic and diagnostic
tabs of the HTML report

## Usage

``` r
pdf_section_progn_diagn(
  report,
  tier = c(1, 2),
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- report:

  PCGR report object

- tier:

  tiers to include

- color:

  main report color (hex)

## Value

markdown with raw Typst tables (character), to be printed with cat() in
a chunk with output 'asis'
