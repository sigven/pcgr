# Process OncoKB CNA output file and fetch complete annotations

Process OncoKB CNA output file and fetch complete annotations

## Usage

``` r
process_oncokb_cna(
  cna_file,
  oncotree_code = NULL,
  oncokb_token = NULL,
  var_calls = NULL,
  rate_limiting_delay = 1
)
```

## Arguments

- cna_file:

  Path to OncoKB-annotated CNA file (TSV format)

- oncotree_code:

  OncoTree code/name (e.g., "SKIN", "BREAST")

- oncokb_token:

  OncoKB API token

- var_calls:

  Data frame with variant calls (used for mapping of variant alteration)

- rate_limiting_delay:

  Delay in seconds between API calls to respect rate limits (default: 1
  seconds)

## Value

List with CNA information and corresponding JSON annotations
