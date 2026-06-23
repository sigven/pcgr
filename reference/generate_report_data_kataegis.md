# Function that generates data frame with potential kataegis events

Function that generates data frame with potential kataegis events

## Usage

``` r
generate_report_data_kataegis(
  variant_set,
  sample_name = "SampleX",
  build = "grch37"
)
```

## Arguments

- variant_set:

  data frame with SNVs/InDels (must contain 'CHROM', 'POS','REF','ALT')

- sample_name:

  name of tumor sample

- build:

  genome assembly (grch37/grch38)
