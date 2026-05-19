# Fetch OncoKB annotation for copy number alteration

Fetch OncoKB annotation for copy number alteration

## Usage

``` r
fetch_oncokb_cna_annotation(
  hugo_symbol = NULL,
  cna_type = NULL,
  oncotree_code = NULL,
  oncokb_token = NULL
)
```

## Arguments

- hugo_symbol:

  Gene symbol

- cna_type:

  Type of CNA: "Amplification" or "Deletion"

- oncotree_code:

  OncoTree code/name (e.g., "SKIN", "BREAST")

- oncokb_token:

  OncoKB API token

## Value

List containing the complete JSON response from OncoKB API
