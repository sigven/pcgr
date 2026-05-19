# Function that sets gnomAD_AF_ABOVE_TOLERATED to TRUE for variants if any gnomAD population frequency exceeds max_tolerated_af

Function that sets gnomAD_AF_ABOVE_TOLERATED to TRUE for variants if any
gnomAD population frequency exceeds max_tolerated_af

## Usage

``` r
assign_germline_popfreq_status(sample_calls, max_tolerated_af = 0.01)
```

## Arguments

- sample_calls:

  data frame with variants

- max_tolerated_af:

  max tolerated germline allele frequency

## Value

sample_calls
