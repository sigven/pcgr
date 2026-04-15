# Plotly Pie Chart - variant statistics

Function that generates a pie chart using plotly for a given category in
a data frame

## Usage

``` r
plotly_pie_chart(
  df_variant_stats = NULL,
  category = "CODING_STATUS",
  plot_margin_top = 50,
  plot_margin_bottom = 20,
  plot_margin_left = 20,
  plot_margin_right = 20,
  font_family = "Helvetica",
  font_size = 15,
  pie_line_width = 3,
  opacity_filtered_categories = 0.4,
  hole_size_pie = 0.4
)
```

## Arguments

- df_variant_stats:

  Data frame with variant statistics

- category:

  Category for pie chart (e.g. CODING_STATUS)

- plot_margin_top:

  Top margin

- plot_margin_bottom:

  Bottom margin

- plot_margin_left:

  Left margin

- plot_margin_right:

  Right margin

- font_family:

  Font family

- font_size:

  Font size

- pie_line_width:

  Line width for pie chart segments

- opacity_filtered_categories:

  Opacity for filtered categories

- hole_size_pie:

  Hole size for pie chart (0 to 1)

## Value

Plotly pie chart object
