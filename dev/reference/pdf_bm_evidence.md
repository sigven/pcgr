# Biomarker evidence tables (therapeutic sensitivity: Tier 1-3, predictive evidence; therapeutic resistance) of a variant category for the PDF report - same evidence selection as the HTML report

Biomarker evidence tables (therapeutic sensitivity: Tier 1-3, predictive
evidence; therapeutic resistance) of a variant category for the PDF
report - same evidence selection as the HTML report

## Usage

``` r
pdf_bm_evidence(
  report,
  variant_category = "snv_indel",
  tier = c(1, 2, 3),
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- report:

  PCGR report object

- variant_category:

  'snv_indel', 'cna', 'fusion', 'msi' or 'tmb'

- tier:

  actionability tiers to include (default: 1-3)

- color:

  main report color (hex)

## Value

list with 'sensitivity' and 'resistance' tinytable objects (NULL when
there is no evidence)
