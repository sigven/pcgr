# Function that tells if the MSI status of the sample is MSI-high

Function that tells if the MSI status of the sample is MSI-high

## Usage

``` r
msi_status_is_high(msi_content = NULL, high_confidence_only = FALSE)
```

## Arguments

- msi_content:

  MSI report element (as returned by `generate_report_data_msi`)

- high_confidence_only:

  logical, if TRUE, return TRUE only if the MSI-H prediction is also of
  high confidence (P(MSI-H) \>= 0.9)

## Value

logical
