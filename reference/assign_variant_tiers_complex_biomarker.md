# AMP/ASCO/CAP tier classification for complex biomarkers (MSI-H, TMB-H)

Assigns a tier of clinical significance to a complex biomarker (MSI-high
or TMB-high), based exclusively on the strength and tumor type
specificity of the biomarker evidence items (CIViC, OncoKB). The
biomarker is not tied to a gene, and there are no variant properties to
consider (as for e.g. oncogenes in fusion partners) - biomarkers without
matching evidence are assigned tier 5.

## Usage

``` r
assign_variant_tiers_complex_biomarker(
  primary_site = "Any",
  biomarker_mapping_confidence = "medium",
  var_df = NULL,
  etype_for_tiering = c("predictive"),
  biomarker_items = NULL
)
```

## Arguments

- primary_site:

  primary tumor site

- biomarker_mapping_confidence:

  confidence level of biomarker mapping (e.g. 'high' or 'medium')

- var_df:

  data frame with the biomarker record (VAR_ID, VARIANT_CLASS,
  ENTREZGENE)

- etype_for_tiering:

  evidence type(s) used for tiering (e.g. 'predictive')

- biomarker_items:

  data frame with biomarker evidence items
