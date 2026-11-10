# Render the gnomAD allele frequency of a germline variant

Renders GNOMAD_AF (exome AF, falling back to genome AF - see
germline.qmd) in scientific notation, or 'absent' when the variant is
not found in gnomAD (AF missing or 0), with both the exome (gnomADe_AF)
and genome (gnomADg_AF) frequencies in the hover tooltip.

## Usage

``` r
render_gnomad_af_cell(data)
```

## Arguments

- data:

  The germline variant data frame (passed via closure)

## Value

A function that can be used as a reactable cell renderer for the
GNOMAD_AF column
