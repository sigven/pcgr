# Function that reads CPSR-classified variants from a TSV file

Function that reads CPSR-classified variants from a TSV file

## Usage

``` r
load_cpsr_classified_variants(
  fname_cpsr_tsv = NA,
  fname_cpsr_yaml = NA,
  cols = NULL,
  ignore_vus = FALSE,
  ref_data = NULL
)
```

## Arguments

- fname_cpsr_tsv:

  Path to raw input file with CPSR-classified SNVs/InDels

- fname_cpsr_yaml:

  Path to YAML configuration file for CPSR analysis

- cols:

  column type definitions of raw input file

- ignore_vus:

  logical indicating if VUS should be ignored in report

- ref_data:

  PCGR reference data object
