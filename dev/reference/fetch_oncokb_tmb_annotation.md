# Fetch OncoKB annotation for TMB-H (tumor mutational burden-high)

Fetch OncoKB annotation for TMB-H (tumor mutational burden-high)

## Usage

``` r
fetch_oncokb_tmb_annotation(
  oncotree_code = NULL,
  oncokb_token = NULL,
  base_api_url = NULL
)
```

## Arguments

- oncotree_code:

  OncoTree code/name (e.g., "COADREAD", "BLCA")

- oncokb_token:

  OncoKB API token

- base_api_url:

  Optional base URL for OncoKB API (default: oncokb_base_api_url)

## Value

List containing the complete JSON response from OncoKB API
