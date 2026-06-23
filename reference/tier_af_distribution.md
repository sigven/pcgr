# Function that plots a histogram of the the variant allelic support (tumor) - grouped by tiers

Function that plots a histogram of the the variant allelic support
(tumor) - grouped by tiers

## Usage

``` r
tier_af_distribution(tier_df, bin_size = 0.05)
```

## Arguments

- tier_df:

  data frame with somatic mutations

- bin_size:

  size of bins for allelic frequency

## Value

p geom_histogram plot from ggplot2
