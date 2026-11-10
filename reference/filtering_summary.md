# Summary of SNV/InDel variant filtering, for the variant filtering sections of the HTML and PDF reports

Combines the depth/allelic support filtering (all assay modes, see
filter_read_support()) and - for tumor-only input - the classification
of variants as likely germline/somatic (SOMATIC_CLASSIFICATION) and the
exclusion of non-exonic variants.

## Usage

``` r
filtering_summary(report)
```

## Arguments

- report:

  PCGR report object

## Value

list with elements

- tumor_only - logical

- n_input - number of variants before filtering

- n_reported - number of variants retained (reported)

- read_support - data frame with the depth/allelic support filtering
  steps (STEP, THRESHOLD, APPLIED, NOTE, N_REMOVED, N_REMAINING), zero
  rows if no thresholds were set

- read_support_removed - actionable (Tier 1/2) variants and mutation
  hotspots removed by depth/allelic support filtering

- germline - (tumor-only) data frame with CLASSIFICATION, N, SHARE

- exonic - (tumor-only, non-exonic variants excluded) data frame with
  CLASSIFICATION, N, SHARE among somatic variants

- germline_removed - (tumor-only) actionable (Tier 1/2) variants
  filtered as likely germline

- tumor_only_settings - tumor-only filtering settings
