# AMP/ASCO/CAP tier classification for somatic variants in cancer

Function that assigns tier classifications to genes subject to somatic
CNA (amplifications and deletions), somatic SNVs/InDels, or RNA fusions,
based on the presence and strength of biomarker evidence associated with
records in the variant set. Tier classification is based on the
AMP/ASCO/CAP guidelines for somatic variant interpretation in cancer (Li
et al., J Mol Diagn. 2017). The function also considers variant
properties associated with oncogenes and tumor suppressor genes (e.g.
low MAF, coding status), which are used to assign tier 3 status to
variants with uncertain clinical significance.

## Usage

``` r
assign_amp_asco_cap_tiers(
  vartype = "snv_indel",
  clinical_significance = "therapeutic_sensitivity",
  primary_site = "Any",
  biomarker_mapping_confidence = "medium",
  var_df = NULL,
  biomarker_items = NULL
)
```

## Arguments

- vartype:

  variant type ('snv_indel', 'cna', 'fusion')

- clinical_significance:

  character indicate the clinical significance types of evidence used
  for tiering of predictive evidence. Possible values:
  "therapeutic_sensitivity", "therapeutic_resistance",
  "prognostic_poor", "prognostic_better", "diagnostic_positive",
  "diagnostic_negative"

- primary_site:

  primary tumor site

- biomarker_mapping_confidence:

  confidence level of variant-biomarker mapping resolution (e.g. 'high'
  or 'medium'). 'High' indicates matches at genomic/hgvsp/hgvsc level,
  'medium' indicates matches at gene/exon level Note that that matches
  at gene/exon level also consider variant properties (loss-of-function,
  hotspot overlap etc) to avoid noise from presumably less impactful
  variants that are likely not relevant for the biomarker relationship

- var_df:

  data frame with variants (SNVs/InDels, CNAs, fusions)

- biomarker_items:

  data frame with biomarker evidence items
