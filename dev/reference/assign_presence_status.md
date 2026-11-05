# Function that assigns a logical status column reflecting whether a variant has a (non-missing) value in a given annotation column

Function that assigns a logical status column reflecting whether a
variant has a (non-missing) value in a given annotation column

## Usage

``` r
assign_presence_status(sample_calls, annotation_col, status_col)
```

## Arguments

- sample_calls:

  data frame with sample variant calls

- annotation_col:

  annotation column (e.g. 'COSMIC_ID'); if not present, sample_calls is
  returned unchanged

- status_col:

  name of the status column to assign
