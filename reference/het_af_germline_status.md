# Function that assigns a logical (STATUS_LIKELY_GERMLINE_HETEROZYGOUS) reflecting whether a variant is likely heterozygous (germline) - based on allelic fraction (VAF_TUMOR), presence in gnomAD and dbSNP, and no presence in TCGA and COSMIC

Function that assigns a logical (STATUS_LIKELY_GERMLINE_HETEROZYGOUS)
reflecting whether a variant is likely heterozygous (germline) - based
on allelic fraction (VAF_TUMOR), presence in gnomAD and dbSNP, and no
presence in TCGA and COSMIC

## Usage

``` r
het_af_germline_status(sample_calls)
```

## Arguments

- sample_calls:

  data frame with sample variant calls
