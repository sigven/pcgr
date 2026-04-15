# Function that generates mutational signatures data for PCGR report

Function that generates mutational signatures data for PCGR report

## Usage

``` r
generate_report_data_signatures(
  variant_set = NULL,
  vstats = NULL,
  ref_data = NULL,
  settings = NULL,
  n_bootstrap_iterations = 200,
  sig_contribution_cutoff = 0.01
)
```

## Arguments

- variant_set:

  Somatic callset (SNV)

- vstats:

  Variant statistics object

- ref_data:

  PCGR reference data object

- settings:

  PCGR run/configuration settings

- n_bootstrap_iterations:

  Number of bootstrap iterations for signature analysis

- sig_contribution_cutoff:

  Minimal signature contribution for reporting
