# Rank RNA fusions for the PDF report: actionability tier (Tier 1 first, no tier last), then recurrence in the tumor site of the sample (Mitelman DB), then the number of split reads

Rank RNA fusions for the PDF report: actionability tier (Tier 1 first,
no tier last), then recurrence in the tumor site of the sample (Mitelman
DB), then the number of split reads

## Usage

``` r
rank_fusions(fusions, site = "Any")
```

## Arguments

- fusions:

  data frame with fusions (fusion callset - variant_display)

- site:

  primary tumor site of the sample

## Value

ordered data frame, with column SITE_RECURRENCE
