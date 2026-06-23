# Build column definitions for oncogenicity table with category-aware styling

This function generates a list of column definitions for a reactable
table that displays oncogenicity annotations. It applies specific
styling to the SYMBOL and SAMPLE_ALTERATION columns based on cancer
association rank and ClinVar classification, respectively. The
ONCOGENICITY column is styled with background colors corresponding to
oncogenicity categories. The function also ensures that only relevant
columns are shown, while others are hidden, and that styling lookup
columns are not displayed in the table.

## Usage

``` r
build_oncogenicity_col_defs(data, primary_cols, symbol_rank_col)
```

## Arguments

- data:

  The data frame containing the oncogenicity annotations and styling
  lookup columns.

- primary_cols:

  Character vector of column names that should be considered primary and
  kept visible in the table.

- symbol_rank_col:

  Name of the column in `data` that contains the cancer association rank
  used for styling the SYMBOL column (default: "TISSUE_ASSOC_RANK").

## Value

A list of column definitions that can be passed to the `colDef`
