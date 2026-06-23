# Function that orders genomic aberrations according to order of chromosomes and chromosomal position

Function that orders genomic aberrations according to order of
chromosomes and chromosomal position

## Usage

``` r
order_variants(vcf_df, chrom_var = "CHROM", pos_var = "POS")
```

## Arguments

- vcf_df:

  data frame

- chrom_var:

  variable name of chromosome in data frame

- pos_var:

  variable name for chromosomal position

## Value

vcf_df data frame with ordered mutations
