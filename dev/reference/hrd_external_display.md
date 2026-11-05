# Function that prepares display strings for an external HRD result

Function that prepares display strings for an external HRD result

## Usage

``` r
hrd_external_display(hrd_content = NULL)
```

## Arguments

- hrd_content:

  HRD report element

## Value

list with elements: 'status', 'caller', 'caller_label' (caller, with the
threshold on the score if the status is derived), 'score_type',
'score_type_label', 'score_title', 'score_value', 'components' (named
vector with the genomic scar components given, possibly empty),
'basis_note' (how the status was obtained) and 'internal_note' (NA
unless the PCGR genomic scar score is also available)
