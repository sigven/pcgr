# Fetch OncoKB annotation for gene fusion

Fetch OncoKB annotation for gene fusion

## Usage

``` r
fetch_oncokb_fusion_annotation(
  hugo_symbol_a = "EML4",
  hugo_symbol_b = "ALK",
  oncotree_code = "LUNG",
  oncokb_token = NULL,
  base_api_url = NULL,
  structural_variant_type = "FUSION",
  is_functional_fusion = TRUE
)
```

## Arguments

- hugo_symbol_a:

  First gene in fusion (5' partner)

- hugo_symbol_b:

  Second gene in fusion (3' partner)

- oncotree_code:

  OncoTree code/name (e.g., "SKIN", "BREAST")

- oncokb_token:

  OncoKB API token

- base_api_url:

  Optional base URL for OncoKB API (default: oncokb_base_api_url)

- structural_variant_type:

  Type of structural variant (default: "FUSION")

- is_functional_fusion:

  Whether fusion is functional (default: TRUE)

## Value

List containing the complete JSON response from OncoKB API
