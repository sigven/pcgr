# Function that reads and validates fully annotated CNA data (segments and genes) from PCGR pre-reporting pipeline

Function that reads and validates fully annotated CNA data (segments and
genes) from PCGR pre-reporting pipeline

## Usage

``` r
load_somatic_cna(
  fname_cna_segment = NULL,
  fname_cna_gene = NULL,
  ref_data = NULL,
  settings = NULL
)
```

## Arguments

- fname_cna_segment:

  Path to file with pre-processed CNA segments

- fname_cna_gene:

  Path to file with pre-processed CNA gene-level data

- ref_data:

  PCGR reference data object

- settings:

  PCGR run/configuration settings
