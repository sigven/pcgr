# Function that detects kataegis events from a data frame with genomic cooordinates of mutations

Function that detects kataegis events from a data frame with genomic
cooordinates of mutations

## Usage

``` r
kataegis_detect(
  data,
  sample_id = "sample_id",
  build = "grch37",
  min.mut = 6,
  max.dis = 1000,
  chr = "chr",
  pos = "pos",
  txdb = NULL
)
```

## Arguments

- data:

  data frame with somatic mutations, as produced by kataegis_input

- sample_id:

  sample identifier

- build:

  genomoe assembly build

- min.mut:

  minimum number of mutations in localized hypermutated region

- max.dis:

  maximum distance of kataegis event (basepairs)

- chr:

  column name in data that denotes chromosome

- pos:

  column name in data that denotes position

- txdb:

  transcript database (txdb)

## Value

kataegis_df data frame with potential kataegis events
