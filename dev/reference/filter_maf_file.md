# Function that takes a MAF file generated with vcf2maf and filters out variants that are presumably germline (tumor-only run)

Function that takes a MAF file generated with vcf2maf and filters out
variants that are presumably germline (tumor-only run)

## Usage

``` r
filter_maf_file(callset, settings)
```

## Arguments

- callset:

  Callset with pre-processed somatic SNV/InDel variants

- settings:

  PCGR run/configuration settings
