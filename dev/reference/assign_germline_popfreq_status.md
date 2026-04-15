# Function that sets STATUS_POPFREQ_1KG_ABOVE_TOLERATED/ STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED to TRUE for variants if any population frequency exceeds max_tolerated_af

Function that sets STATUS_POPFREQ_1KG_ABOVE_TOLERATED/
STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED to TRUE for variants if any
population frequency exceeds max_tolerated_af

## Usage

``` r
assign_germline_popfreq_status(
  sample_calls,
  pop = "NFE",
  dbquery = c("gnomADe", "gnomADg"),
  max_tolerated_af = 0.01
)
```

## Arguments

- sample_calls:

  data frame with variants

- pop:

  population code (1000 Genomes/gnomAD)

- dbquery:

  character vector with germline db sources, e.g. c("gnomADe","gnomADg")

- max_tolerated_af:

  max tolerated germline allele frequency

## Value

sample_calls
