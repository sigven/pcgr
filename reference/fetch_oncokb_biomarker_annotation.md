# Fetch OncoKB annotation for a complex biomarker (MSI-H or TMB-H)

Microsatellite instability-high (MSI-H) and tumor mutational burden-high
(TMB-H) are annotated by OncoKB as "atypical alterations", through the
protein change endpoint. No gene is provided - OncoKB maps these
biomarkers to the pseudo-gene "Other Biomarkers". The response has the
same structure as for other alteration types (summaries, treatments with
levels of evidence). The evidence extraction helpers (e.g.
`extract_complete_annotation`) support these biomarkers with
`vartype = "msi"` or `vartype = "tmb"`.

## Usage

``` r
fetch_oncokb_biomarker_annotation(
  biomarker = "MSI-H",
  oncotree_code = NULL,
  oncokb_token = NULL,
  base_api_url = NULL
)
```

## Arguments

- biomarker:

  Complex biomarker, either "MSI-H" or "TMB-H"

- oncotree_code:

  OncoTree code/name (e.g., "COADREAD", "BLCA")

- oncokb_token:

  OncoKB API token

- base_api_url:

  Optional base URL for OncoKB API (default: oncokb_base_api_url)

## Value

List containing the complete JSON response from OncoKB API
