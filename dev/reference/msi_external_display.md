# Function that prepares display strings for an external MSI result

Function that prepares display strings for an external MSI result

## Usage

``` r
msi_external_display(msi_content = NULL, classifier_ref = "below")
```

## Arguments

- msi_content:

  MSI report element

- classifier_ref:

  where the details of the PCGR MSI classifier prediction can be found
  in the report (e.g. 'below', or 'the tab "PCGR MSI classifier"'), used
  in the notes on concordance. If NULL, no reference is made.

## Value

list with elements: 'status', 'caller', 'caller_label' (caller, with the
threshold on the score if the status is derived), 'score_title',
'score_value', 'loci_value' (NA if not available), 'basis_note' (how the
status was obtained), 'zero_loci' (logical), 'concordance' (see
`msi_concordance`), 'discordance_note' (NA unless discordant) and
'classifier_class' and 'classifier_detail' (class and probability of the
classifier prediction)
