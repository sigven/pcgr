# Function that maps biomarker identifiers from VCF variant annotation to full biomarker data tables for display in report

Function that maps biomarker identifiers from VCF variant annotation to
full biomarker data tables for display in report

## Usage

``` r
map_biomarker_data(
  varcalls = NULL,
  ref_data = NULL,
  settings = NULL,
  variant_origin = "Somatic",
  vartype = "snv_indel"
)
```

## Arguments

- varcalls:

  variant calls data frame with biomarker identifiers

- ref_data:

  reference data object containing biomarker data

- settings:

  PCGR settings object

- variant_origin:

  variant origin for filtering biomarker evidence items (e.g. "Somatic",
  "Germline", "Any")

- vartype:

  variant type for filtering biomarker evidence items (e.g. "snv_indel",
  "cna", "fusion")

## Value

data frame with biomarker data for report display
