# Clean extracted evidence data from OncoKB

Clean extracted evidence data from OncoKB

## Usage

``` r
clean_oncokb_evidence(
  evidence_df = NULL,
  gene = NA,
  oncokb_root_gene = NA,
  alteration = NA,
  vartype = NA,
  profile_name = NA,
  oncotree_code = NA
)
```

## Arguments

- evidence_df:

  data.frame with OncoKB evidence

- gene:

  gene symbol

- oncokb_root_gene:

  root gene symbol from OncoKB annotation (for reference)

- alteration:

  name of alteration

- vartype:

  variant type (e.g., "snv_indel", "fusion", "cna")

- profile_name:

  molecular profile name (for display)

- oncotree_code:

  OncoTree code for tumor type

## Value

Cleaned data.frame with standardized fields for evidence level, clinical
significance, and molecular profile formatting
