# Get OncoKB therapeutic evidence items for a complex biomarker (MSI-H or TMB-H)

Queries the OncoKB API for the biomarker and returns the evidence items
in the PCGR biomarker evidence data model (columns with prefix `BM_`),
as for other alteration types.

## Usage

``` r
fetch_oncokb_biomarker_eitems(
  biomarker = "MSI-H",
  oncotree_code = NULL,
  oncokb_token = NULL,
  base_api_url = NULL,
  all_tumor_types = FALSE
)
```

## Arguments

- biomarker:

  Complex biomarker, either "MSI-H" or "TMB-H"

- oncotree_code:

  OncoTree code of the tumor (NULL/NA: query not restricted to a tumor
  type)

- oncokb_token:

  OncoKB API token

- base_api_url:

  Optional base URL for OncoKB API (default: oncokb_base_api_url)

- all_tumor_types:

  logical, also add evidence for other tumor types (second query without
  tumor type, see `oncokb_query_all_tumor_types`)

## Value

Data frame with evidence items (empty if none, or if the query failed)
