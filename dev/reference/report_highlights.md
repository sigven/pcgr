# Key findings of a PCGR report (highlights) as markdown

The most critical findings come first: somatic alterations of strong or
potential clinical significance (AMP/ASCO/CAP Tier 1/2). These are
followed by supporting findings, each included only when informative:
oncogenic/likely oncogenic SNV/InDel variants, composite genomic
biomarkers (TMB, MSI), copy number variants (amplified proto-oncogenes,
deleted tumor suppressors, genome doubling, HRD score), mutational
signatures, kataegis, germline findings (CPSR), and transcriptome
findings (RNA fusions, expression outliers).

## Usage

``` r
report_highlights(report, max_show = 5)
```

## Arguments

- report:

  PCGR report object

- max_show:

  maximum number of items listed per finding

## Value

list with 'lead' (markdown on actionable alterations), 'pills' (positive
composite biomarkers), 'warnings' (markdown on actionable variants
removed by variant filtering, "" if none) and 'supporting' (markdown
bullet list, "" if none)
