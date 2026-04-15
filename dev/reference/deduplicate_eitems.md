# Function that removes redundancy in variant evidence items (i.e. if a variant is assicated with evidence at the codon level, evidence at the exon/gene level is ignored)

Function that removes redundancy in variant evidence items (i.e. if a
variant is assicated with evidence at the codon level, evidence at the
exon/gene level is ignored)

## Usage

``` r
deduplicate_eitems(
  var_eitems = NULL,
  target_type = "exact",
  target_other = c("codon", "exon", "gene")
)
```

## Arguments

- var_eitems:

  data frame with variant evidence items

- target_type:

  which resolution level should be used as the "best" level ('exact' or
  'codon)

- target_other:

  resolution levels for other evidence items that should be ignored if
  evidence is found at the target_type level
