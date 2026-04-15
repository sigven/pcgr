# Function that makes plots of mutational signature contributions in a given sample (both ggplot and plotly)

Function that makes plots of mutational signature contributions in a
given sample (both ggplot and plotly)

## Usage

``` r
plot_signature_contributions(
  signature_contributions = NULL,
  per_signature = TRUE
)
```

## Arguments

- signature_contributions:

  A list with two data frames: 'per_group' and 'per_signature'.

- per_signature:

  Logical. If TRUE, the plot will show the contribution per signature.
  If FALSE, the plot will show the contribution per group.
