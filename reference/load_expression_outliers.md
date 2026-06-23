# Load expression outlier results

Load expression outlier results

## Usage

``` r
load_expression_outliers(
  settings = NULL,
  ref_data = NULL,
  percentile_cutoff_high = 95,
  percentile_cutoff_low = 5,
  z_score_cutoff = 1.5
)
```

## Arguments

- settings:

  PCGR run/configuration settings

- ref_data:

  PCGR reference data object

- percentile_cutoff_high:

  numeric, percentile cutoff for high expression

- percentile_cutoff_low:

  numeric, percentile cutoff for low expression

- z_score_cutoff:

  numeric, z-score cutoff for expression outliers
