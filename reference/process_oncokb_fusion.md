# Process OncoKB fusion output file and fetch complete annotations

Process OncoKB fusion output file and fetch complete annotations

## Usage

``` r
process_oncokb_fusion(
  fusion_file = NULL,
  var_calls = NULL,
  oncotree_code = NULL,
  oncokb_token = NULL,
  rate_limiting_delay = 1
)
```

## Arguments

- fusion_file:

  Path to OncoKB-annotated fusion file (TSV format)

- var_calls:

  Data frame with variant calls (used for mapping of variant alteration)

- oncotree_code:

  OncoTree code/name (e.g., "SKIN", "BREAST")

- oncokb_token:

  OncoKB API token

- rate_limiting_delay:

  Delay in seconds between API calls to respect rate limits (default: 1
  second)

## Value

List with fusion information and corresponding JSON annotations
