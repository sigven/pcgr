# A function that detects whether the sample name in variant data frame is unique (as present in column name VCF_SAMPLE_ID), throws an error if multiple sample names are present for the CPSR workflow

A function that detects whether the sample name in variant data frame is
unique (as present in column name VCF_SAMPLE_ID), throws an error if
multiple sample names are present for the CPSR workflow

## Usage

``` r
detect_vcf_sample_name(df, sample_name = NULL, cpsr = FALSE)
```

## Arguments

- df:

  VCF data frame

- sample_name:

  name of sample identifier

- cpsr:

  logical indicating CPSR workflow

## Value

df Vranges object
