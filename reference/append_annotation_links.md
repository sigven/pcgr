# Function that appends multiple HTML annotation links to variant identifiers e.g. COSMIC, CLINVAR, REFSEQ etc

Function that appends multiple HTML annotation links to variant
identifiers e.g. COSMIC, CLINVAR, REFSEQ etc

## Usage

``` r
append_annotation_links(var_data_df, vartype = "snv_indel", skip = NULL)
```

## Arguments

- var_data_df:

  data frame with variant entries

- vartype:

  'snv_indel' or 'cna', or 'exp'

- skip:

  elements to be ignored during annotation
