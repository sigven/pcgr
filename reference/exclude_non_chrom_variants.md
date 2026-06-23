# Function that excludes genomic aberrations from non-nuclear chromosomes

Function that excludes genomic aberrations from non-nuclear chromosomes

## Usage

``` r
exclude_non_chrom_variants(vcf_df, chrom_var = "CHROM")
```

## Arguments

- vcf_df:

  data frame

- chrom_var:

  variable name of chromosome in data frame

## Value

vcf_df data frame with mutations from nuclear chromosomes only
