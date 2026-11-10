# Append flag for segments with imputed allele-specific copy number

Ensures that column 'CN_MAJOR_MINOR_IMPUTED' is present (logical) in a
data frame with copy number segments - segments with an observed total
copy number, but imputed major/minor split (optional input column
'nMajor_nMinor_Imputed'). A missing column (or NA) means not imputed.

## Usage

``` r
append_cna_imputed_flag(cna_segment = NULL)
```

## Arguments

- cna_segment:

  data frame with copy number segments
