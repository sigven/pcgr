# Function that generates a pie chart for germline filtering statistics for callsets coming from tumor-only sequencing

Function that generates a pie chart for germline filtering statistics
for callsets coming from tumor-only sequencing

## Usage

``` r
plot_filtering_stats_germline(
  report = NULL,
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

- report:

  list object with PCGR report content

- plot_margin_top:

  top margin of the plot

- plot_margin_bottom:

  bottom margin of the plot

- plot_margin_left:

  left margin of the plot

- plot_margin_right:

  right margin of the plot

- font_family:

  font family for plot text

- font_size:

  font size for plot text

- pie_line_width:

  line width for pie chart segments

- opacity_filtered_categories:

  opacity for filtered categories

- hole_size_pie:

  size of the hole in the pie chart

## Value

filtering_stats list with data frame and plotly pie chart
