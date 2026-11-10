# Function that loads the normalised result of an external MSI caller

Reads the (single-row) tab-separated file with the MSI result of an
external caller (e.g. MSIsensor-pro), as written by the PCGR Python
workflow (columns SAMPLE_ID, MSI_STATUS, MSI_STATUS_BASIS, MSI_CALLER,
CALLER_VERSION, MSI_SCORE, MSI_SCORE_TYPE, MSI_SCORE_THRESHOLD_USED,
N_LOCI_EVALUATED, N_LOCI_UNSTABLE).

## Usage

``` r
load_external_msi(fname = NULL)
```

## Arguments

- fname:

  name of the file with the normalised external MSI result

## Value

data frame with one row
