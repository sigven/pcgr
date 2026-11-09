The *Key findings* summary is generated from a fixed template, so that its structure and ordering are
identical for all samples. Only findings present in the data are listed, and each list shows at most
five items (remaining items are indicated as *+n more*). Genes within a list are always ordered by the
same rules, with ties resolved alphabetically by gene symbol.

**Leading findings**

- __Therapeutically actionable alterations__ (AMP/ASCO/CAP Tier 1, then Tier 2) among SNVs/InDels, copy number
  aberrations and RNA fusions. Tiers refer to _therapeutic_ evidence only (as intended by the
  AMP/ASCO/CAP guidelines); prognostic and diagnostic evidence is not tiered. Within a tier:
  tissue-specific cancer association rank (highest first), then global cancer association rank.
  If none are found, this is stated explicitly. When a primary tumor site is specified, therapeutic
  evidence is matched against that site: strong evidence (levels A/B) in the same tumor type or in a pan-cancer
  biomarker is Tier 1, whereas strong evidence in other tumor types, or weak evidence (levels C/D/E) in the same
  tumor type, is Tier 2. When no primary tumor site is specified
  (site code 0, "Any"), evidence from _all_ tumor types is considered rather than pan-cancer evidence only:
  strong evidence (levels A/B) in a pan-cancer biomarker is Tier 1, strong evidence in a tumor-type-specific
  biomarker (any tumor type) is Tier 2, and weak evidence (levels C/D/E) is Tier 3 and therefore not listed here.
- __Diagnostic evidence__: when a primary tumor site is specified, alterations with _strong_ diagnostic
  evidence (CIViC levels A/B; OncoKB level Dx1) matching that site (e.g. an _EWSR1::FLI1_ fusion in
  soft tissue tumors) are listed in a section of their own, directly below the therapeutically actionable
  alterations, one item per hit and with the cancer types of the evidence (at most two shown per item, and
  at most three items; further hits are listed among the supporting findings). Hits are ranked by the
  resolution of the biomarker match: exact matches (e.g. fusions, copy number aberrations, variants matched
  at the HGVS or codon level) before gene- or region-level matches; ties are resolved by cancer association
  ranks, then gene symbol.
  A hit shown here is not repeated among the supporting findings. No diagnostic evidence is shown when no primary
  tumor site is specified (site code 0, "Any"), as evidence from any tumor type is then considered -
  such alterations are listed among the supporting findings. Diagnostic evidence is not tiered
- __Composite genomic biomarkers__: a tumor classified as *TMB-high* (with its TMB estimate, in
  mutations/Mb), *MSI-high* (high-confidence predictions only, P(MSI-H) >= 0.9) and/or
  *HR-deficient* (external HRD result only, `--input_hrd`) is marked with a highlighted label
  directly below the therapeutically actionable alterations. Nothing is shown when the tumor is neither TMB-high,
  MSI-high nor HR-deficient (TMB, MSI and HRD status are reported in full in the *Composite genomic
  biomarkers* section)

- __Filtered actionable variants__: actionable SNVs/InDels (Tier 1/2) or mutation hotspots that were
  removed by depth/allelic support filtering, and - for tumor-only input - actionable variants filtered
  as likely germline, are flagged with a warning (they are not considered elsewhere in the report; see
  *Variant filtering* in the SNV/InDel section)

**Supporting findings** (in this order)

- __Prognostic biomarkers__: alterations (SNVs/InDels, copy number aberrations, RNA fusions) with
  evidence of poor or better outcome, each with the cancer types of its evidence (CIViC disease or
  OncoKB tumor type, at most three). The cancer types are the main information element: prognostic
  and diagnostic evidence is not tiered (tiers refer to therapeutic evidence only)
- __Diagnostic biomarkers__: alterations with evidence supporting a diagnosis, each with the cancer
  types of its evidence (shown and ordered as for prognostic biomarkers). Alterations are listed
  here also when they have no therapeutic relevance at all. Strong hits matching a specified tumor site
  are shown as diagnostic labels in the leading findings instead (see above), and are not repeated here
- With a primary tumor site specified, only prognostic and diagnostic evidence matching that site
  is considered (strong or weak). When no primary tumor site is specified (site code 0, "Any"),
  evidence from any tumor type is considered, including weak evidence. Strong evidence is defined as
  CIViC levels A/B and OncoKB levels Dx1/Px1; weak evidence as CIViC levels C/D/E and OncoKB levels
  Dx2/Dx3/Px2/Px3
- __SNV/InDel variants__: number of oncogenic / likely oncogenic SNVs/InDels (all, including those
  listed as actionable - their number is stated), with the first five of the variants *not* already
  listed as actionable (oncogenicity score, highest first, then tissue-specific and global cancer
  association rank)
- __Copy number variants__: amplified proto-oncogenes and deep deletions of tumor suppressors (ordered by
  cancer association rank), whole-genome doubling, HRD score
- __Mutational signatures__: the (at most two) proposed aetiology groups with the largest
  contributions, each contributing at least 10% (signatures of *Other* aetiologies are not listed)
- __Kataegis__: number of potential events
- __Germline (CPSR)__: pathogenic before likely pathogenic variants, then cancer association rank
- __Transcriptome findings__: number of RNA fusions, with proto-oncogene partners listed by the number of split
  reads; overexpressed proto-oncogenes and underexpressed tumor suppressors (expression outliers),
  (selected among the expression outliers, i.e. genes in the extreme percentiles of the reference
  cohort) ordered by cancer association rank, then absolute Z-score; the sample's percentile in the
  reference cohort is given in parentheses (rounded, e.g. *99th pctile*)
