# Fetch OncoKB annotation for SNV/InDel via protein change

Fetch OncoKB annotation for SNV/InDel via protein change

## Usage

``` r
fetch_oncokb_hgvsp_annotation(
  hugo_symbol = "BRAF",
  protein_change = "V600E",
  oncotree_code = "SKIN",
  oncokb_token = NULL,
  base_api_url = NULL,
  reference_genome = "GRCh38"
)
```

## Arguments

- hugo_symbol:

  Gene symbol (e.g., "BRAF")

- protein_change:

  Protein change in short format (e.g., "V600E")

- oncotree_code:

  OncoTree code/name (e.g., "SKIN", "BREAST")

- oncokb_token:

  OncoKB API token

- base_api_url:

  Optional base URL for OncoKB API (default: oncokb_base_api_url)

- reference_genome:

  Genome build, either "GRCh37" or "GRCh38" (default: "GRCh38")

## Value

List containing the complete JSON response from OncoKB API
