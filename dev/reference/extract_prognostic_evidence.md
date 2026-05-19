# Extract prognostic implications from OncoKB annotation

Extract prognostic implications from OncoKB annotation

## Usage

``` r
extract_prognostic_evidence(
  oncokb_annotation = NULL,
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

  OncoKB annotation list

- gene:

  Gene symbol (e.g., "BRAF")

- alteration:

  Alteration description

- vartype:

  Variant type (e.g., "snv_indel", "fusion", "cna") (e.g., HGVSp short
  format for SNVs/InDels)

- oncotree_code:

  OncoTree code for tumor type (e.g., "BRCA" for breast cancer)

- variant_id:

  Variant identifier

- match_by:

  Matching strategy for SNVs/InDels (e.g., "hgvsp", "genomic")

## Value

Tibble with prognostic evidence items
