# Order the variants of an oncogenicity table (SNVs/InDels)

Therapeutically actionable variants (tier 1, then tier 2) are listed
first, followed by all other variants. Within each of these groups,
variants are ranked by oncogenicity score and cancer association ranks.
VUS are further ranked by the presence of a pathogenic ClinVar
classification, within the positions held by VUS. The table is limited
after the actionable variants are put first, so that these are not lost
for samples with many variants.

## Usage

``` r
order_oncogenicity_table(d, n_max)
```

## Arguments

- d:

  data frame with variants

- n_max:

  maximum number of variants to keep

## Value

ordered data frame, where ACTIONABILITY_TIER holds tiers 1 and 2 only
(see
[`actionable_tier()`](https://sigven.github.io/pcgr/reference/actionable_tier.md)),
and is placed next to ONCOGENICITY
