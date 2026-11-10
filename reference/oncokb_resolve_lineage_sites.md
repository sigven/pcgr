# Primary site of OncoKB evidence for the queried OncoTree code

OncoKB evidence items are tiered by their primary site (OncoTree tissue
of the cancer type). When the OncoTree code used for OncoKB queries has
another primary site than the query tumor (e.g. 'MPNST' (Peripheral
Nervous System) for tumor site 'Soft Tissue'), the evidence for the
code's own cancer type would be considered evidence from another tumor
site. Items whose cancer type is in the lineage of the OncoTree code
(the code itself, or one of its ancestors, e.g. 'Nerve Sheath Tumor')
are therefore assigned the primary site of the query tumor.

## Usage

``` r
oncokb_resolve_lineage_sites(
  biomarker_items = NULL,
  primary_site = "Any",
  oncotree_lineage = NULL
)
```

## Arguments

- biomarker_items:

  data.frame with biomarker evidence items

- primary_site:

  primary site of query tumor

- oncotree_lineage:

  character vector with OncoTree names of the queried OncoTree code and
  its ancestors

## Value

data.frame with biomarker evidence items
