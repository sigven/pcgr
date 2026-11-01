# Function that compares the external MSI result with the PCGR MSI classifier

Function that compares the external MSI result with the PCGR MSI
classifier

## Usage

``` r
msi_concordance(msi_content = NULL)
```

## Arguments

- msi_content:

  MSI report element

## Value

"concordant", "discordant", "not_evaluable" (the external result is
indeterminate), or NA (no external result and/or no classifier
prediction)
