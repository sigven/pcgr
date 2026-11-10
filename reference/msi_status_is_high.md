# Function that tells if the MSI status of the sample is MSI-high

If a result from an external MSI caller is available (MSI-H or MSS),
this determines the MSI status (an external MSI-H call is regarded as
confident, as the confidence of the call is the responsibility of the
caller). An external MSS call means that the sample is not regarded as
MSI-high, also when the PCGR MSI classifier predicts MSI-H. If no
external result is available (or the external result is indeterminate),
the prediction of the PCGR MSI classifier is used.

## Usage

``` r
msi_status_is_high(msi_content = NULL, high_confidence_only = FALSE)
```

## Arguments

- msi_content:

  MSI report element (as returned by `generate_report_data_msi`)

- high_confidence_only:

  logical, if TRUE, return TRUE only if the MSI-H prediction is also of
  high confidence (P(MSI-H) \>= 0.9) - only relevant for the PCGR MSI
  classifier

## Value

logical
