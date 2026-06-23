# Process OncoKB MAF output files (both HGVSp and HGVSg) and fetch complete annotations

Process OncoKB MAF output files (both HGVSp and HGVSg) and fetch
complete annotations

## Usage

``` r
process_oncokb_maf(
  maf_file_hgvsp = NULL,
  maf_file_hgvsg = NULL,
  var_calls = NULL,
  oncotree_code = NULL,
  oncokb_token = NULL,
  oncokb_base_api_url = NULL,
  rate_limiting_delay = 1
)
```

## Arguments

- maf_file_hgvsp:

  Path to OncoKB-annotated MAF file using HGVSp (protein change)
  (optional)

- maf_file_hgvsg:

  Path to OncoKB-annotated MAF file using HGVSg (genomic change)
  (optional)

- var_calls:

  Data frame with variant calls (used for mapping of variant alteration)

- oncotree_code:

  OncoTree code used for OncoKB annotation

- oncokb_token:

  OncoKB API token

- oncokb_base_api_url:

  Optional base URL for OncoKB API (default: oncokb_base_api_url)

- rate_limiting_delay:

  Delay in seconds between API calls to respect rate limits (default: 1
  second)

## Value

Data frame with variant information and corresponding JSON annotations
