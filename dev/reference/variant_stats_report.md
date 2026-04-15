# Function that generate stats for a given variant set, considering number of variants/genes affected across tiers, types of variants ()

Function that generate stats for a given variant set, considering number
of variants/genes affected across tiers, types of variants ()

## Usage

``` r
variant_stats_report(callset = NULL, name = "vstats", vartype = "snv_indel")
```

## Arguments

- callset:

  list object with callset data (CNA or SNVs/InDels)

- name:

  type of variant statistic

- vartype:

  type of variant ('snv_indel', 'cna')
