# Whether there is prognostic or diagnostic biomarker evidence to show for a variant category (tiers as given by progn_diagn_tiers()) - based on the evidence items, as the tier classification also holds variants with a Tier 3 by gene role only (no evidence of this type)

Whether there is prognostic or diagnostic biomarker evidence to show for
a variant category (tiers as given by progn_diagn_tiers()) - based on
the evidence items, as the tier classification also holds variants with
a Tier 3 by gene role only (no evidence of this type)

## Usage

``` r
has_progn_diagn_evidence(report, variant_category, clnsig)
```

## Arguments

- report:

  PCGR report object

- variant_category:

  'snv_indel', 'cna' or 'fusion'

- clnsig:

  'prognostic_poor', 'prognostic_better', 'diagnostic_positive' or
  'diagnostic_negative'

## Value

logical
