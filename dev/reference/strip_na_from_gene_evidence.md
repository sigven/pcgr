# Remove missing-value placeholders from gene evidence strings

Entries in the bundle's cancer gene evidence strings can be literal "NA"
tokens (e.g. CancerMine citations that lack a link). These are removed,
along with any evidence label left with no entries. Text inside citation
links is left untouched.

## Usage

``` r
strip_na_from_gene_evidence(evidence = NULL)
```

## Arguments

- evidence:

  character vector of evidence strings

## Value

character vector of cleaned strings (NA if nothing remains)
