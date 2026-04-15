# Function that loads specific set of clinical variant evidence items (CIViC + CGI) based on given parameters (mutation type, variant origin, tumor type etc)

Function that loads specific set of clinical variant evidence items
(CIViC + CGI) based on given parameters (mutation type, variant origin,
tumor type etc)

## Usage

``` r
load_eitems(
  eitems_raw = NULL,
  ontology = NULL,
  alteration_types = c("MUT"),
  origin = "Somatic",
  tumor_type_specificity = NULL,
  tumor_type = NULL
)
```

## Arguments

- eitems_raw:

  complete set of clinical variant evidence items

- ontology:

  phenotype ontology data frame

- alteration_types:

  types of alteration ('MUT', 'CNA', 'MUT_LOF')

- origin:

  variant origin ('Somatic','Germline')

- tumor_type_specificity:

  tumor type specificity ('any', 'specific')

- tumor_type:

  primary tumor site

## Value

eitems variant evidence items
