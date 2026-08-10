# Re-synchronize biomarker evidence items/classifications with a (possibly further-filtered) variant set

Biomarker evidence matching
([`map_biomarker_data()`](https://sigven.github.io/pcgr/reference/map_biomarker_data.md),
[`assign_amp_asco_cap_tiers()`](https://sigven.github.io/pcgr/reference/assign_amp_asco_cap_tiers.md))
is performed once, early, against the variant set as it existed right
after annotation. Callers may go on to remove variants from that set
afterwards (e.g. allelic depth/fraction filtering via
[`filter_read_support()`](https://sigven.github.io/pcgr/reference/filter_read_support.md),
or germline/non-exonic filtering for tumor-only input). Without
re-syncing, evidence items for variants removed by such downstream
filtering remain "orphaned" - present in `bm_evidence` but absent from
the variant set - which surfaces as biomarker-matched variants missing
from the variant listing in the report while still appearing in the
biomarker evidence listing (see
<https://github.com/sigven/pcgr/issues/302>).

## Usage

``` r
sync_biomarker_evidence(bm_evidence = NULL, var_df = NULL)
```

## Arguments

- bm_evidence:

  list with biomarker evidence data (as initialized by
  [`init_biomarker_content()`](https://sigven.github.io/pcgr/reference/init_biomarker_content.md)),
  i.e. top-level 'eitems'/'classification' data frames plus one sub-list
  per clinical significance category, each with its own
  'eitems'/'classification' data frames

- var_df:

  data frame with the final (filtered) variant set

## Value

bm_evidence list, with all 'eitems'/'classification' data frames limited
to records matching a variant in var_df
