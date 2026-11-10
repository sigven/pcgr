# Classification values as colored pills in a PDF report table - colors as in the HTML report

Classification values as colored pills in a PDF report table - colors as
in the HTML report

## Usage

``` r
pdf_pill_class(
  x,
  type = c("tier", "oncogenicity", "pathogenicity", "elevel"),
  clnsig = "therapeutic_sensitivity"
)
```

## Arguments

- x:

  values per table row: actionability tier (1-3) for 'tier',
  classification for 'oncogenicity' (e.g. 'Likely_Oncogenic') and
  'pathogenicity' (e.g. 'Likely Pathogenic'), evidence level code for
  'elevel' (e.g. 'A', '3B')

- type:

  'tier' (color_palette\$tier_sensitivity / tier_resistance),
  'oncogenicity', 'pathogenicity' or 'elevel' (color_palette\$elevel)

- clnsig:

  for type 'tier': 'therapeutic_sensitivity' or 'therapeutic_resistance'

## Value

character vector with Typst markup (see pdf_pill); unknown values in
grey
