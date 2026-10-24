# Oncogenicity classification as a sortable table value

reactable sorts a column by its values, which would order the
oncogenicity classes alphabetically (Benign first). Prefixes each class
with its rank in `color_palette$oncogenicity$levels` (e.g.
"1\|Oncogenic", "5\|Benign"), so that sorting (ascending) lists the most
oncogenic variants first. The prefix is stripped for display by
[`render_oncogenicity_cell()`](https://sigven.github.io/pcgr/dev/reference/render_oncogenicity_cell.md),
and ignored by the column filter of
[`build_oncogenicity_col_defs()`](https://sigven.github.io/pcgr/dev/reference/build_oncogenicity_col_defs.md).

## Usage

``` r
oncogenicity_sort_value(x, color_palette = pcgrr::color_palette)
```

## Arguments

- x:

  oncogenicity classification (character)

- color_palette:

  PCGR color palette (oncogenicity levels)

## Value

character vector - values not among the levels (and NA) are returned
unchanged
