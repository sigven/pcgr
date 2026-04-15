# Function that makes a plot with TMB boxplots for reference cohorts, highlighting the TMB estimate for a given sample and the cohort/primary site of interest

Function that makes a plot with TMB boxplots for reference cohorts,
highlighting the TMB estimate for a given sample and the cohort/primary
site of interest

## Usage

``` r
plot_tmb_primary_site_tcga(
  tmb_reference,
  p_site = "Liver",
  tmb_estimates = NULL,
  tmb_display_type = "missense_only",
  tumor_only = FALSE
)
```

## Arguments

- tmb_reference:

  data frame with TMB estimates for reference samples (e.g. TCGA)

- p_site:

  primary tumor_site (sample)

- tmb_estimates:

  data frame with estimates of mutational burden (sample)

- tmb_display_type:

  TMB estimation used for display (e.g. missense_only)

- tumor_only:

  logical, if TRUE color using tumor-only color
