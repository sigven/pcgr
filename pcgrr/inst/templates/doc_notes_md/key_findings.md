The *Key findings* summary is generated from a fixed template, so that its structure and ordering are
identical for all samples. Only findings present in the data are listed, and each list shows at most
five items (remaining items are indicated as *+n more*). Genes within a list are always ordered by the
same rules, with ties resolved alphabetically by gene symbol.

**Leading findings**

- __Actionable alterations__ (AMP/ASCO/CAP Tier 1, then Tier 2) among SNVs/InDels, copy number
  aberrations and RNA fusions. Within a tier: tissue-specific cancer association rank (highest first),
  then global cancer association rank. If none are found, this is stated explicitly.
- __Composite genomic biomarkers__: a tumor classified as *TMB-high* (with its TMB estimate, in
  mutations/Mb) and/or *MSI-high* is marked with a highlighted label directly below the actionable
  alterations. Nothing is shown when the tumor is neither TMB-high nor MSI-high (TMB and MSI status
  are reported in full in the *Composite genomic biomarkers* section)

- __Filtered actionable variants__: actionable SNVs/InDels (Tier 1/2) or mutation hotspots that were
  removed by depth/allelic support filtering, and - for tumor-only input - actionable variants filtered
  as likely germline, are flagged with a warning (they are not considered elsewhere in the report; see
  *Variant filtering* in the SNV/InDel section)

**Supporting findings** (in this order)

- __Prognostic biomarkers__: alterations (SNVs/InDels, copy number aberrations, RNA fusions) with
  AMP/ASCO/CAP Tier 1/2 evidence of poor or better outcome, each with its tier (Tier 1 first, then
  cancer association rank)
- __Diagnostic biomarkers__: alterations with Tier 1/2 evidence supporting a diagnosis, each with
  its tier (ordered as for prognostic biomarkers)
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
