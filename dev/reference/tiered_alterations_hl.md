# Somatic alterations with an actionability tier (SNVs/InDels, CNAs, RNA fusions) as display labels

Somatic alterations with an actionability tier (SNVs/InDels, CNAs, RNA
fusions) as display labels

## Usage

``` r
tiered_alterations_hl(report, tier = c(1, 2))
```

## Arguments

- report:

  PCGR report object

- tier:

  actionability tiers

## Value

data frame with LABEL, SYMBOL, ALTERATION and ACTIONABILITY_TIER
