# Kataegis events (localized hypermutation) for the PDF report - region, chromosome arm, length, number of mutations, proportion of C\>X mutations and confidence

Kataegis events (localized hypermutation) for the PDF report - region,
chromosome arm, length, number of mutations, proportion of C\>X
mutations and confidence

## Usage

``` r
pdf_table_kataegis(events, color = pcgrr::color_palette$report_color$values[1])
```

## Arguments

- events:

  data frame with kataegis events (content - kataegis - events)

- color:

  main report color (hex)

## Value

tinytable object, or NULL when there are no events
