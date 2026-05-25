# Function that writes contents of PCGR object to various output formats (Rmarkdown/flexdashboard HTML reports, JSON, tab-separated etc)

param report List object with all report data (PCGR/CPSR), settings etc.
param tier_model type of tier model param output_format contents/file
format of output (html/json/tsv/cna_tsv etc) param flexdb logical
indicating if HTML output should be dashboard

## Usage

``` r
write_report_tsv(report = NULL, output_type = "snv_indel")
```

## Arguments

- report:

  List object with all report data, settings etc.

- output_type:

  character indicating output type for TSV, i.e. 'snv_indel',
  'snv_indel_unfiltered', 'cna_gene', 'fusion', or 'msigs'

## Details

export Function that writes contents of PCGR object to flexdashboard
HTML reports

param report List object with all report data, settings etc.

export Function that writes contents of PCGR object to a TSV file
