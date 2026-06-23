# Function that generates dense and tiered annotated variant datasets

Function that generates dense and tiered annotated variant datasets

## Usage

``` r
generate_tier_tsv(variant_set, config, annotation_tags, sample_name = "test")
```

## Arguments

- variant_set:

  List with tiered variants

- config:

  PCGR configuration settings

- annotation_tags:

  List with display columns

- sample_name:

  Sample identifier

## Value

tsv_variants data frame with tier-annotated list of variants for
tab-separated output
