# Function that filters variant set on depth and allelic fraction according to settings provided by user (tumor and control)

The filtering summary is attached to the result as attribute
'read_support_filtering' (see
[`attr()`](https://rdrr.io/r/base/attr.html)): a list with 'n_input'
(number of variants before filtering), 'n_retained', 'steps' (data frame
with one row per configured filter: STEP, THRESHOLD, APPLIED, NOTE,
N_REMOVED, N_REMAINING) and 'actionable_removed' (data frame with Tier
1/2 variants and mutation hotspots removed by the filters: SYMBOL,
ALTERATION, CONSEQUENCE, HGVSP, ACTIONABILITY_TIER, MUTATION_HOTSPOT,
DP_TUMOR, VAF_TUMOR)

## Usage

``` r
filter_read_support(vcf_df, config = NULL)
```

## Arguments

- vcf_df:

  data frame with variants

- config:

  list with PCGR configuration settings

## Value

vcf_df (filtered), with attribute 'read_support_filtering'
