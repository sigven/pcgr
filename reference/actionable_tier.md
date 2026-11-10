# Keep only the therapeutically actionable tiers (1 and 2) of an actionability tier vector

Used for the oncogenicity tables, where tier 1 and 2 alterations are
listed first and flagged. Lower tiers (uncertain significance) are not
treated as actionable here.

## Usage

``` r
actionable_tier(x)
```

## Arguments

- x:

  actionability tier

## Value

integer vector with values 1, 2 or NA
