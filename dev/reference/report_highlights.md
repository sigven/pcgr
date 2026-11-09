# Key findings of a PCGR report (highlights) as markdown

The most critical findings come first: somatic alterations of strong or
potential clinical significance (AMP/ASCO/CAP Tier 1/2). These are
followed by supporting findings, each included only when informative:
prognostic and diagnostic biomarkers (listed with the cancer types of
their evidence; untiered - see progn_diagn_tiers() for the evidence
considered; strong diagnostic hits matching a specified tumor site are
instead highlighted as pills, see diagnostic_pills_hl()),
oncogenic/likely oncogenic SNV/InDel variants, composite genomic
biomarkers (TMB, MSI, external HRD status), copy number variants
(amplified proto-oncogenes, deleted tumor suppressors, genome doubling,
HRD score), mutational signatures, kataegis, germline findings (CPSR),
and transcriptome findings (RNA fusions, expression outliers).

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

list with 'lead' (markdown on therapeutically actionable alterations
and, when a tumor site is specified, diagnostic evidence - strong
evidence only), 'pills' (positive composite biomarkers), 'pill_targets'
(section label per pill), 'pill_kinds' ('composite' per pill),
'warnings' (markdown on actionable variants removed by variant
filtering, "" if none) and 'supporting' (markdown bullet list, "" if
none)
