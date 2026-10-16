# Donut chart of variant statistics (static, ggplot2) for the PDF report - counterpart of plotly_pie_chart() in the HTML report

Donut chart of variant statistics (static, ggplot2) for the PDF report -
counterpart of plotly_pie_chart() in the HTML report

## Usage

``` r
pdf_plot_donut(
  df_variant_stats,
  category = "CODING_STATUS",
  title = NULL,
  center_label = NULL,
  show_counts = TRUE,
  color_palette = pcgrr::color_palette
)
```

## Arguments

- df_variant_stats:

  data frame with 'N', 'Pct' and the category column, from
  stats_type_snv_indel()

- category:

  category column (e.g. 'CODING_STATUS')

- title:

  chart title

- center_label:

  text in the center of the donut (default: total N)

- show_counts:

  show counts (N) with the legend labels

- color_palette:

  PCGR color palette

## Value

ggplot2 object
