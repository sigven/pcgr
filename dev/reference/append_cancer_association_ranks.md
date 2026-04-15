# Function that appends cancer gene evidence links

Function that appends cancer gene evidence links

## Usage

``` r
append_cancer_association_ranks(
  var_data_df = NULL,
  ref_data = NULL,
  primary_site = "Any"
)
```

## Arguments

- var_data_df:

  Data frame of sample variants from VCF

- ref_data:

  PCGR/CPSR reference data object

- primary_site:

  Primary tumor site

## Value

vcf_data_df
