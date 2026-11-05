# Function that joins targeted drugs (any label, and on-label for the primary site) to a list of variants, by gene (SYMBOL)

Function that joins targeted drugs (any label, and on-label for the
primary site) to a list of variants, by gene (SYMBOL)

## Usage

``` r
join_targeted_drugs(
  var_df,
  primary_site = "Breast",
  ref_data = NULL,
  drop_cols = NULL,
  empty_cols = NULL
)
```

## Arguments

- var_df:

  data frame with variants (columns VAR_ID and SYMBOL)

- primary_site:

  primary tumor sample site

- ref_data:

  PCGR/CPSR reference data object

- drop_cols:

  drug columns removed after the join

- empty_cols:

  drug columns set to NA when no variant has targeted drugs
