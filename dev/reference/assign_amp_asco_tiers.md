# Function that assigns tier classifications to somatic CNA segments and SNVs/InDels, based on the presence of biomarker evidence found in the variant set

Function that assigns tier classifications to somatic CNA segments and
SNVs/InDels, based on the presence of biomarker evidence found in the
variant set

## Usage

``` r
assign_amp_asco_tiers(
  vartype = "snv_indel",
  primary_site = "Any",
  variants_df = NULL,
  biomarker_items = NULL
)
```

## Arguments

- vartype:

  variant type ('snv_indel' or 'cna')

- primary_site:

  primary tumor site

- variants_df:

  data frame with variants (SNVs/InDels or CNAs)

- biomarker_items:

  data frame with biomarker evidence items
