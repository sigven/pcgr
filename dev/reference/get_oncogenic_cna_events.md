# Get oncogenic copy number events

Function that extracts oncogenic copy number events from a data frame of
annotated transcripts (within copy number segments), utilizing
oncogene/tumor suppressor status and copy number variant class
(gain/loss). This set is used to highlight any copy-number altered
transcripts that are presumably oncogenic (e.g. oncogene amplifications,
tumor suppressor deletions)

## Usage

``` r
get_oncogenic_cna_events(cna_df_display = NULL)
```

## Arguments

- cna_df_display:

  data frame with transcript annotations per copy number segment
