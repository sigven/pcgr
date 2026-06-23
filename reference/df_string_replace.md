# Function that performs stringr::str_replace on strings of multiple string columns of a dataframe

Function that performs stringr::str_replace on strings of multiple
string columns of a dataframe

## Usage

``` r
df_string_replace(df, strings, pattern, replacement, replace_all = FALSE)
```

## Arguments

- df:

  data frame

- strings:

  name of columns for which string replace is to be performed

- pattern:

  pattern to replace

- replacement:

  string to replace

- replace_all:

  logical - replace all occurrences

## Value

df
