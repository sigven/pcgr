# Extract a scalar value from a list or vector returned by the OncoKB API

Safely coerces whatever the OncoKB JSON parser returns for a single
field (NULL, a length-1 list/vector, or a multi-element list/vector)
into a single character scalar suitable for storing in a data frame
column.

## Usage

``` r
oncokb_extract_field(x)
```

## Arguments

- x:

  Object returned by the OncoKB API for a single field.

## Value

A length-1 character scalar, or `NA` when `x` is NULL.
