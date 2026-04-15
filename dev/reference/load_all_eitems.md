# Function that loads all evidence items from CIViC and CGI, and combines them in a unified data.frame

Function that loads all evidence items from CIViC and CGI, and combines
them in a unified data.frame

## Usage

``` r
load_all_eitems(eitems_raw = NULL, alteration_type = "MUT", origin = "Somatic")
```

## Arguments

- eitems_raw:

  raw data frame with evidence items

- alteration_type:

  type of alteration ('MUT','MUT_LOF','CNA')

- origin:

  variant origin ('Germline','Somatic')

## Value

all_eitems
