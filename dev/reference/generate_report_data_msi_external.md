# Function that generates report data for an external MSI result

Function that generates report data for an external MSI result

## Usage

``` r
generate_report_data_msi_external(settings = NULL)
```

## Arguments

- settings:

  PCGR run/configuration settings (the normalised external MSI result is
  given by `settings$molecular_data$fname_msi_tsv`)

## Value

list with elements 'eval', 'result' (one-row data frame) and
'threshold_user_defined' (logical: the threshold on the MSI score used
for deriving the MSI status was set by the user,
`--input_msi_threshold`)
