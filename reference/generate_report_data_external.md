# Function that generates report data for the result of an external assay/tool (MSI or HRD)

Function that generates report data for the result of an external
assay/tool (MSI or HRD)

## Usage

``` r
generate_report_data_external(settings = NULL, type = c("msi", "hrd"))
```

## Arguments

- settings:

  PCGR run/configuration settings (the normalised external result is
  given by `settings$molecular_data$fname_<type>_tsv`)

- type:

  type of result, 'msi' or 'hrd'

## Value

list with elements 'eval', 'result' (one-row data frame) and
'threshold_user_defined' (logical: the threshold on the score used for
deriving the status was set by the user, `--input_<type>_threshold`)
