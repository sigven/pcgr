# Biomarker databases whose evidence items define actionability tiers - when OncoKB is run, CGI evidence is shown for tiered aberrations, but not used for tier assignment (cf. biomarker_sources_not_tiered()). Used in messages on the absence of biomarker evidence.

Biomarker databases whose evidence items define actionability tiers -
when OncoKB is run, CGI evidence is shown for tiered aberrations, but
not used for tier assignment (cf. biomarker_sources_not_tiered()). Used
in messages on the absence of biomarker evidence.

## Usage

``` r
bm_db_tier_label(report, fusion = FALSE)
```

## Arguments

- report:

  PCGR report object

- fusion:

  logical, RNA fusions (CGI evidence is never used for fusions)

## Value

character, e.g. 'CIViC and OncoKB'
