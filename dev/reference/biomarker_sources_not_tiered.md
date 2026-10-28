# Biomarker sources not used for tier assignment

When OncoKB is run, CGI evidence items are shown, but not used for tier
assignment: the CGI biomarker set is no longer updated, and is partly
mapped to generic cancer types (e.g. 'Cancer'), which would otherwise be
considered pan-cancer evidence. CIViC and OncoKB are used on equal
terms.

## Usage

``` r
biomarker_sources_not_tiered(settings = NULL)
```

## Arguments

- settings:

  PCGR run/configuration settings

## Value

character vector with sources (BM_SOURCE_DB), or NULL
