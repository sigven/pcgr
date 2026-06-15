# Function that reads and validates an annotated somatic SNV/InDel file from PCGR pre-reporting pipeline

Function that reads and validates an annotated somatic SNV/InDel file
from PCGR pre-reporting pipeline

## Usage

``` r
load_somatic_snv_indel(
  fname = NA,
  ref_data = NULL,
  settings = NULL,
  simulate_vaf_dp = FALSE
)
```

## Arguments

- fname:

  Path to file with pre-processed somatic SNV/InDel variants

- ref_data:

  PCGR reference data object

- settings:

  PCGR run/configuration settings

- simulate_vaf_dp:

  Internal/test use only. If TRUE and VAF_TUMOR is entirely missing,
  replace it with random values drawn from Uniform(0.01, 0.99). Never
  set this in production runs.
