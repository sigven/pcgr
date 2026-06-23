# Function that sorts chromosomal segments according to chromosome and chromosomal start/end position

Function that sorts chromosomal segments according to chromosome and
chromosomal start/end position

## Usage

``` r
sort_chromosomal_segments(
  df,
  chromosome_column = "CHROM",
  start_segment = "POS",
  end_segment = "POS"
)
```

## Arguments

- df:

  data frame with chromosome and start + end segment

- chromosome_column:

  name of column for chromosome name is sigven

- start_segment:

  name of column that indicates start of chromosomal segment

- end_segment:

  name of column that indicates end of chromosomal segment

## Value

df_final data frame with sorted chromosomal segments
