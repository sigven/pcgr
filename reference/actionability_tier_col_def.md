# Column definition for the tier badge of the oncogenicity tables

The column is hidden if no alteration in the table is actionable.

## Usage

``` r
actionability_tier_col_def(data, color_palette = pcgrr::color_palette)
```

## Arguments

- data:

  table with column ACTIONABILITY_TIER, as returned by
  [`actionable_tier()`](https://sigven.github.io/pcgr/reference/actionable_tier.md)

- color_palette:

  PCGR color palette (tier colors)

## Value

a reactable column definition
