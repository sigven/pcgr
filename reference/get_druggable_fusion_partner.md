# Identify druggable fusion partners

This function identifies druggable fusion partners from a given data
frame of fusion genes. It uses a reference data bundle to append
targeted drug annotations based on the specified primary site.

## Usage

``` r
get_druggable_fusion_partner(
  df = NULL,
  partner = "5P",
  ref_data = NULL,
  variant_display = FALSE,
  primary_site = "Lung"
)
```

## Arguments

- df:

  A data frame containing fusion gene information with columns "VAR_ID"
  and "FUSION_GENE_5P" or "FUSION_GENE_3P".

- partner:

  A character string indicating whether to analyze the 5' or 3 ' fusion
  partner. Default is "5P".

- ref_data:

  A list containing reference data, including drug annotations.

- variant_display:

  A logical value indicating whether to use variant-level drug
  annotations. Default is FALSE.

- primary_site:

  A character string specifying the primary site for drug annotation.
  Default is "Lung".

## Value

A data frame with druggable fusion partners and their associated
targeted drug annotations.
