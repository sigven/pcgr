# Function that retrieves prevalent signatures for a given tumor type/primary site Data is collected from COSMIC v3.4.

Function that retrieves prevalent signatures for a given tumor
type/primary site Data is collected from COSMIC v3.4.

## Usage

``` r
get_prevalent_site_signatures(
  site = "Any",
  custom_collection = NULL,
  ref_data = NULL,
  min_prevalence_pct = 0.1,
  incl_poss_artifacts = T
)
```

## Arguments

- site:

  Primary tumor site

- custom_collection:

  Custom collection of signatures from COSMIC

- ref_data:

  PCGR reference data object

- min_prevalence_pct:

  Minimum prevalence (pct) of signature in cohorts associated with
  primary site - used to select reference signatures for inclusion in
  signature reconstruction

- incl_poss_artifacts:

  logical indicating if artefact signatures are to be included
