# Biomarker evidence of Tier 1/2 aberrations (SNVs/InDels, CNAs, RNA fusions) for the actionability section of the PDF report - therapeutic sensitivity and resistance tables per variant category

Biomarker evidence of Tier 1/2 aberrations (SNVs/InDels, CNAs, RNA
fusions) for the actionability section of the PDF report - therapeutic
sensitivity and resistance tables per variant category

## Usage

``` r
pdf_section_bm_evidence(
  report,
  tier = c(1, 2),
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- report:

  PCGR report object

- tier:

  actionability tiers to include

- color:

  main report color (hex)

## Value

markdown with raw Typst tables (character), to be printed with cat() in
a chunk with output 'asis'
