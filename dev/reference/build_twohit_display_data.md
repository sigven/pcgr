# Build display data for potential two-hit events (nested reactable)

For CNA records flagged with somatic and/or germline loss-of-function
SNV/InDel candidates, this function parses the comma-separated
`TWOHIT_CANDIDATE_SOMATIC` / `TWOHIT_CANDIDATE_GERMLINE` columns
(format: `VAR_ID;CONSEQUENCE` per entry) and joins them against the
somatic and germline variant callsets to produce two data frames
suitable for a nested reactable:

## Usage

``` r
build_twohit_display_data(
  cna_variant = NULL,
  snv_somatic = NULL,
  snv_germline = NULL
)
```

## Arguments

- cna_variant:

  data frame - gene-level CNA callset
  (`pcg_report$content$cna$callset$variant`)

- snv_somatic:

  data frame - somatic SNV/InDel callset
  (`pcg_report$content$snv_indel$callset$variant`)

- snv_germline:

  data frame - germline classified callset
  (`pcg_report$content$germline_classified$callset$variant`)

## Value

Named list with elements `main` and `nested` (both data frames, empty if
no two-hit candidates found).

## Details

- `main`:

  One row per two-hit CNA gene: SYMBOL, GENENAME, VARIANT_CLASS,
  CN_TOTAL, LOH, and a hidden `.row_id` key.

- `nested`:

  One row per overlapping LoF variant: `.row_id` (FK), ORIGIN (Somatic /
  Germline), ALTERATION, CONSEQUENCE, VAF_TUMOR (somatic only),
  CLASSIFICATION.
