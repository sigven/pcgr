# Key findings of a PCGR report (highlights) as markdown

The most critical findings come first: somatic alterations of strong or
potential clinical significance (AMP/ASCO/CAP Tier 1/2). These are
followed by supporting findings, each included only when informative:
prognostic and diagnostic biomarkers (Tier 1/2), oncogenic/likely
oncogenic SNV/InDel variants, composite genomic biomarkers (TMB, MSI),
copy number variants (amplified proto-oncogenes, deleted tumor
suppressors, genome doubling, HRD score), mutational signatures,
kataegis, germline findings (CPSR), and transcriptome findings (RNA
fusions, expression outliers).

## Usage

``` r
report_highlights(report, max_show = 5, links = FALSE)
```

## Arguments

- report:

  PCGR report object

- max_show:

  maximum number of items listed per finding

- links:

  finding titles as links to the corresponding report sections (PDF
  report - section labels as in the pcgr_quarto_report/pdf templates,
  e.g. 'pcgr-cna'); each finding is only included when its section is
  part of the report

## Value

list with 'lead' (markdown on actionable alterations), 'pills' (positive
composite biomarkers), 'pill_targets' (section label per pill),
'warnings' (markdown on actionable variants removed by variant
filtering, "" if none) and 'supporting' (markdown bullet list, "" if
none)
