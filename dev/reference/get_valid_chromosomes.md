# Checks for valid chromosome names in data frame of variants

Checks for valid chromosome names in data frame of variants

## Usage

``` r
get_valid_chromosomes(vcf_data_df, chromosome_column = "CHROM", bsg = NULL)
```

## Arguments

- vcf_data_df:

  data frame

- chromosome_column:

  name of chromosome column

- bsg:

  BSGenome object

## Value

vcf_data_df valid data frame with valid mutations
