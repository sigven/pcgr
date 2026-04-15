# Function that plots a histogram of the the variant allelic support (tumor)

Function that plots a histogram of the the variant allelic support
(tumor)

## Usage

``` r
af_distribution(var_df, bin_size = 0.05)
```

## Arguments

- var_df:

  data frame with somatic mutations

- bin_size:

  size of bins for allelic frequency

## Value

p geom_histogram plot from ggplot2
