# Labels of actionable variants/mutation hotspots removed by filtering

Labels of actionable variants/mutation hotspots removed by filtering

## Usage

``` r
filtered_variant_labels(d, markdown = TRUE)
```

## Arguments

- d:

  data frame with removed variants (SYMBOL, ALTERATION/HGVSP,
  ACTIONABILITY_TIER, MUTATION_HOTSPOT, SOMATIC_CLASSIFICATION,
  VAF_TUMOR, DP_TUMOR - all optional except SYMBOL)

- markdown:

  bold gene symbols (markdown)

## Value

character vector, e.g. "**BRAF** p.V600E (Tier 1; AF 0.04, DP 12)"
