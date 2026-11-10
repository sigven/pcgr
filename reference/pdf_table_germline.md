# Germline variants in cancer predisposition genes (CPSR classification) for the PDF report - genotype, clinical significance (colored pills), ACMG evidence codes and mode of inheritance

Germline variants in cancer predisposition genes (CPSR classification)
for the PDF report - genotype, clinical significance (colored pills),
ACMG evidence codes and mode of inheritance

## Usage

``` r
pdf_table_germline(
  variants,
  low_depth_threshold = 10,
  color = pcgrr::color_palette$report_color$values[1]
)
```

## Arguments

- variants:

  data frame with classified germline variants (content -
  germline_classified - callset - variant_display)

- low_depth_threshold:

  variants with a control (germline) sequencing depth below this value
  are flagged (only when depth is known)

- color:

  main report color (hex)

## Value

tinytable object
