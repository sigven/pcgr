# Function that generates data for rainfall plot (mutation density along genome, considering SNVs only)

Function that generates data for rainfall plot (mutation density along
genome, considering SNVs only)

## Usage

``` r
generate_report_data_rainfall(
  variant_set,
  colors = NULL,
  autosomes = FALSE,
  build = NULL
)
```

## Arguments

- variant_set:

  data frame with SNVs/InDels (must contain "CHROM", "POS","REF","ALT")

- colors:

  character vector of six color codes (denoting the different mutation
  types)

- autosomes:

  logical indicating if plotting should only consider autosomes

- build:

  genome assembly (grch37/grch38)
