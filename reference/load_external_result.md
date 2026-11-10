# Function that loads the normalised result of an external assay/tool

Reads a (single-row) tab-separated file with the result of an external
assay/tool (e.g. MSI, HRD), as written by the PCGR Python workflow. All
columns are read as character, and the given numeric/integer columns are
converted. Missing values are given as '.', 'NA' or empty.

## Usage

``` r
load_external_result(
  fname = NULL,
  required_cols = NULL,
  num_cols = NULL,
  int_cols = NULL,
  label = "external"
)
```

## Arguments

- fname:

  name of the file with the normalised external result

- required_cols:

  columns that must be present

- num_cols:

  columns converted to numeric

- int_cols:

  columns converted to integer

- label:

  type of result, used in error messages (e.g. "MSI", "HRD")

## Value

data frame with one row
