# Function that generates biomarker data for a complex biomarker (MSI-H or TMB-H) for the PCGR report

Evidence items are collected from CIViC/CGI (reference data) and from
OncoKB (web API, if enabled), and classified into tiers of clinical
significance (AMP/ASCO/CAP) for therapeutic sensitivity and resistance,
in the same manner as for other variant types (SNVs/InDels, CNAs,
fusions).

## Usage

``` r
generate_report_data_complex_biomarker(
  biomarker = "MSI-H",
  ref_data = NULL,
  settings = NULL
)
```

## Arguments

- biomarker:

  "MSI-H" or "TMB-H"

- ref_data:

  PCGR reference data object

- settings:

  PCGR run/configuration settings

## Value

list with the biomarker record ('variant', 'variant_display') and
biomarker evidence ('bm_evidence'), i.e. a 'callset' for the biomarker
