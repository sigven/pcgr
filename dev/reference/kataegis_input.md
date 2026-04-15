# Function that detects kataegis events from a data frame with genomic cooordinates of mutations

Function that detects kataegis events from a data frame with genomic
cooordinates of mutations

## Usage

``` r
kataegis_input(
  variant_set,
  chr = "chr",
  pos = "pos",
  ref = "ref",
  alt = "alt",
  build = NULL,
  context_size = 10
)
```

## Arguments

- variant_set:

  data frame with raw set of somatic mutations

- chr:

  column name in data that denotes chromosome

- pos:

  column name in data that denotes position

- ref:

  column name in data that denotes reference allele

- alt:

  column name in data that denotes alternate allele

- build:

  genome build (grch37 or hg38)

- context_size:

  size of neighbouring sequence context
