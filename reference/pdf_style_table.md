# Style the header row of a PDF report table

Dark header with white, bold text, white separator lines between the
header columns, and light zebra striping of the body rows; cell contents
left-aligned and centered vertically

## Usage

``` r
pdf_style_table(
  tt,
  color = pcgrr::color_palette$report_color$values[1],
  header_color = NULL,
  stripe_color = NULL
)
```

## Arguments

- tt:

  tinytable object

- color:

  main report color (hex), see pdf_report_color() - header and stripe
  colors are derived from it (darkened/lightened) unless given
  explicitly

- header_color:

  background color of the header row (hex); default: main report color,
  darkened

- stripe_color:

  background color of every second body row (hex); default: main report
  color, lightened

## Value

styled tinytable object
