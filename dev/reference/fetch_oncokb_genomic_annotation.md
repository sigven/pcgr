# Fetch OncoKB annotation for SNV/InDel via genomic change

Fetch OncoKB annotation for SNV/InDel via genomic change

## Usage

``` r
fetch_oncokb_genomic_annotation(
  hgvsg = "7:g.140753336A>T",
  oncotree_code = NULL,
  variant_origin = "somatic",
  oncokb_token = NULL,
  base_api_url = NULL,
  reference_genome = "GRCh38"
)
```

## Arguments

- hgvsg:

  Genomic change in HGVSg format (e.g., "7:g.140753336A\>T")

- oncotree_code:

  Tumor type name

- variant_origin:

  somatic/germline

- oncokb_token:

  OncoKB API token

- base_api_url:

  Optional base URL for OncoKB API (default: pcgrr::oncokb_base_api_url)

- reference_genome:

  Genome build, either "GRCh37" or "GRCh38" (default: "GRCh38")

## Value

List containing the complete JSON response from OncoKB API
