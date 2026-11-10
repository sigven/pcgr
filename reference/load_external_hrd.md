# Function that loads the normalised result of an external HRD assay/tool

Reads the (single-row) tab-separated file with the HRD result, as
written by the PCGR Python workflow (columns SAMPLE_ID, HRD_STATUS,
HRD_STATUS_BASIS, HRD_CALLER, CALLER_VERSION, HRD_SCORE_TYPE,
HRD_SCORE_DESCRIPTION, HRD_SCORE, HRD_SCORE_THRESHOLD_USED, HRD_LOH,
HRD_LST, HRD_TAI).

## Usage

``` r
load_external_hrd(fname = NULL)
```

## Arguments

- fname:

  name of the file with the normalised external HRD result

## Value

data frame with one row
