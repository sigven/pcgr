# Actionable RNA fusions as records for the PDF report tables (gene pair as 'SYMBOL', as for SNVs/InDels and CNAs)

Actionable RNA fusions as records for the PDF report tables (gene pair
as 'SYMBOL', as for SNVs/InDels and CNAs)

## Usage

``` r
fusion_tier_records(fusions)
```

## Arguments

- fusions:

  data frame with fusions (fusion callset 'variant')

## Value

data frame with SYMBOL, ALTERATION, VARIANT_CLASS and ACTIONABILITY_TIER
(integer), and ONCOGENICITY_OKB when OncoKB annotations are present
