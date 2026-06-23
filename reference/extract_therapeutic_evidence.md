# Extract therapeutic evidence items from OncoKB annotation

Extract therapeutic evidence items from OncoKB annotation

## Usage

``` r
extract_therapeutic_evidence(
  oncokb_annotation,
  gene = NA,
  alteration = NA,
  vartype = NA,
  oncotree_code = NA,
  variant_id = NA,
  match_by = "hgvsp"
)
```

## Arguments

- oncokb_annotation:

  OncoKB annotation list (from API)

- gene:

  Gene symbol (e.g., "BRAF")

- alteration:

  Alteration description (e.g., HGVSp short format for SNVs/InDels)

- vartype:

  Variant type (e.g., "snv_indel", "fusion", "cna")

- oncotree_code:

  OncoTree code for tumor type (e.g., "BRCA" for breast cancer)

- variant_id:

  Variant identifier for tracking

- match_by:

  Matching strategy for SNVs/InDels (e.g., "hgvsp", "genomic")

## Value

Tibble with one row per treatment-cancer type combination
