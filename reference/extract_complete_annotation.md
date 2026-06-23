# Extract mutation effect information

param oncokb_annotation OncoKB annotation list param variant_id Variant
identifier return Tibble with mutation effect data export Extract
oncogenicity classification

## Usage

``` r
extract_complete_annotation(
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

  OncoKB annotation list

- gene:

  gene symbol / gene fusion

- alteration:

  alteration description (e.g., HGVSp short format for SNVs/InDels)

- vartype:

  variant type (e.g., "snv_indel", "fusion", "cna")

- oncotree_code:

  OncoTree code for tumor type

- variant_id:

  Variant identifier

- match_by:

  Matching strategy for SNVs/InDels (e.g., "hgvsp", "genomic")

## Value

Named list with all extracted components

## Details

param oncokb_annotation OncoKB annotation list param variant_id Variant
identifier return Tibble with oncogenicity data export Extract clinical
summaries

param oncokb_annotation OncoKB annotation list param variant_id Variant
identifier return Tibble with summary text fields export Extract
complete structured annotation from OncoKB JSON
