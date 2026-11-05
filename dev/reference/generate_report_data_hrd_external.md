# Function that generates report data for an external HRD result

Function that generates report data for an external HRD result

## Usage

``` r
generate_report_data_hrd_external(settings = NULL)
```

## Arguments

- settings:

  PCGR run/configuration settings (the normalised external HRD result is
  given by `settings$molecular_data$fname_hrd_tsv`)

## Value

list with elements 'eval', 'result' (one-row data frame) and
'threshold_user_defined' (logical: the HRD status was derived from the
score with a user-defined threshold, `--input_hrd_threshold`)
