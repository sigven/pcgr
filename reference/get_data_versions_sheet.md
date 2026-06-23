# Build the DATA_VERSIONS sheet for the Excel workbook

Returns a data frame with columns DATABASE / VERSION / DESCRIPTION / URL
/ LICENSE covering all reference datasets used in the PCGR/CPSR run,
mirroring the "Dataset versions" callout in the HTML Documentation
section.

## Usage

``` r
get_data_versions_sheet(
  report = NULL,
  wflow_pattern = "pcgr",
  tool_label = "PCGR"
)
```

## Arguments

- report:

  PCGR or CPSR report object

- wflow_pattern:

  regex pattern to filter the wflow column (default `"pcgr"`; use
  `"cpsr"` for CPSR reports)

- tool_label:

  label shown in the bundle row DATABASE column (default `"PCGR"`)

## Value

data.frame
