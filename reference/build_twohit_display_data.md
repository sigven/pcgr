# Build display data for potential two-hit events (nested reactable)

For CNA records flagged with somatic and/or germline loss-of-function
SNV/InDel candidates, this function parses the comma-separated
`TWOHIT_CANDIDATE_SOMATIC` / `TWOHIT_CANDIDATE_GERMLINE` columns.

## Usage

``` r
build_twohit_display_data(
  cna_variant = NULL,
  snv_somatic = NULL,
  snv_germline = NULL,
  settings = NULL
)
```

## Arguments

- cna_variant:

  data frame - gene-level CNA callset
  (`pcg_report$content$cna$callset$variant`)

- snv_somatic:

  data frame - somatic SNV/InDel callset (retained for backward
  compatibility; no longer used for the somatic join)

- snv_germline:

  data frame - germline classified callset
  (`pcg_report$content$germline_classified$callset$variant`)

- settings:

  PCGR settings list (`pcg_report$settings`); used to apply
  `tumor_dp_min` and `tumor_af_min` thresholds

## Value

Named list with elements `main` and `nested` (both data frames, empty if
no two-hit candidates found).

## Details

Somatic entry format (semicolon-separated):
`VAR_ID;CONSEQUENCE;VAF_FLAG;ALTERATION;VAF_TUMOR_PCT;ONCOGENICITY`

Display fields for somatic variants are embedded by the Python pipeline
from the unfiltered somatic callset, so variants below the depth filter
still render correctly without a secondary join.

- `main`:

  One row per two-hit CNA gene: SYMBOL, GENENAME, VARIANT_CLASS,
  CN_TOTAL, LOH, and a hidden `.row_id` key.

- `nested`:

  One row per overlapping LoF variant: `.row_id` (FK), ORIGIN (Somatic /
  Germline), ALTERATION, CONSEQUENCE, VAF_GENOTYPE, VAF_FLAG,
  CLASSIFICATION.
