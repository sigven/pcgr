# Color encodings for report elements of PCGR/CPSR

Color encodings for report elements of PCGR/CPSR

## Usage

``` r
color_palette
```

## Format

A list object with different report elements that are color-coded in
PCGR/CPSR reports. Each list element have two vectors: 'levels' and
'values'. Currently, the following list elements are included:

- *pathogenicity* - Colors for five-level pathogenicity levels (CPSR)

- *clinical_evidence* - Colors for strength of evidence of
  cancer-variant associations (A-E)

- *tier* - Colors for tier levels for variant prioritization (PCGR)

- *report_color* - Colors for PCGR assay mode (tumor-control vs.
  tumor-only)

- *warning* - Color for warning (low confidence in PCGR analysis output)

- *success* - Color for success (no evident uncertainty in PCGR analysis
  output)
