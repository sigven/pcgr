# A function that generates a HTML link for selected identifiers (DBSNP, COSMIC, CLINVAR, ENTREZ)

A function that generates a HTML link for selected identifiers (DBSNP,
COSMIC, CLINVAR, ENTREZ)

## Usage

``` r
generate_annotation_link(
  var_df,
  vardb = "DBSNP",
  vartype = "snv_indel",
  group_by_var = c("VAR_ID", "ENTREZGENE"),
  url_prefix = "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",
  link_key_var = "DBSNP_RSID",
  link_display_var = "DBSNP_RSID"
)
```

## Arguments

- var_df:

  data frame

- vardb:

  type of database

- vartype:

  'snv_indel' or 'cna'

- group_by_var:

  variable used for grouping (VAR_ID)

- url_prefix:

  url prefix for link generation

- link_key_var:

  variable used in url for linking

- link_display_var:

  variable used in url for display

## Value

df_annotation_links
