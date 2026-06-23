# FAQ

Frequently asked questions regarding PCGR usage and functionality:

**1. I do not see any data related to allelic depth/support in my
report. I thought that PCGR can grab this information automatically from
my VCF?**

*Answer: VCF variant genotype data (i.e. AD/DP) is something that you as
a user need to specify explicitly when running PCGR. In our experience,
there is currently no uniform way that variant callers format these
types of data (allelic fraction/depth, tumor/normal) in the VCF, and
this makes it very challenging for PCGR to automatically grab this
information from any VCF. Please take a careful look at the example VCF
files (`examples` folder) that comes with PCGR for how PCGR expects this
information to be formatted, and make sure your VCF is formatted
accordingly. There is also an in-depth explanation on the matter
[described
here](https://sigven.github.io/pcgr/articles/input.html#formatting-of-allelic-depthsupport-dpad)*

**2. Is it possible to utilize PCGR for analysis of multiple samples?**

*Answer: As the name of the tool implies, PCGR was developed for the
detailed analysis of individual tumor samples. However, if you take
advantage of the different outputs from PCGR, it can also be utilized
for analysis of multiple samples. First, make sure your input files are
organized per sample (i.e. one VCF file per sample, one CNA file per
sample), so that they can be fed directly to PCGR. Now, once all samples
have been processed with PCGR, note that all the tab-separated output
files (i.e. annotated SNVs, gene copy numbers, fusions) contain the
sample identifier, which enable them to be aggregated and suitable for a
downstream multi-sample analysis. Also note the multi-sheet Excel
workbook, which contains numerous outputs from PCGR, and can be
processed to aggregate findings across samples.*

**3. I do not see the expected transcript-specific consequence for a
particular variant. In what way is the primary variant consequence
established?**

*Answer: PCGR relies upon*
[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) *for
consequence prioritization, in which a specific transcript-specific
consequence is chosen as the primary variant consequence. In the PCGR
configuration file, you may customise how this is chosen by changing the
order of criteria applied when choosing a primary consequence block -
parameter*
[vep_pick_order](https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options)

**4. Is it possible to use RefSeq as the underlying gene transcript
model in PCGR?**

*Answer: PCGR uses GENCODE as the primary gene transcript model, but we
provide cross-references to corresponding RefSeq transcripts when this
is available.*

**5. I have a VCF with structural variants (SVs) detected in my tumor
sample, can PCGR process those as well?**

*Answer: This is currently not supported as input for PCGR, but is
something we want to incorporate in the future. We have a skeleton of SV
support working for CPSR, focusing to support large, multi-exon
deletions*

**6. Is it possible to see all the individual cancer subtypes that
belong to each of the 30 different tumor sites?**

*Answer: Yes, see* [an overview of phenotypes associated with primary
tumor
sites](https://sigven.github.io/pcgr/articles/primary_tumor_sites.md).
See also the related GitHub repository
[phenOncoX](https://github.com/sigven/phenOncoX)

**7. Is it possible for the users to update the data bundle to get the
most recent versions of all underlying data sources?**

*Answer: As of now, the data bundle is updated only with each release of
PCGR. The data harmonization pipeline of knowledge databases in PCGR
contain numerous and complex procedures, with several cleaning, quality
control, and re-formatting steps, and is semi-automated in its present
form. The versions of all databases and key software elements are
outlined in each PCGR report.*

**8. When OncoKB is enabled, I sometimes see variants where PCGR’s
internal `LOSS_OF_FUNCTION` flag disagrees with OncoKB’s mutation effect
(e.g. OncoKB says “Likely Loss-of-function” but `LOSS_OF_FUNCTION` is
FALSE, or vice versa). Why?**

*Answer: The two annotations are derived from independent evidence
sources and should be interpreted as complementary rather than
contradictory. PCGR’s internal `LOSS_OF_FUNCTION` flag is rule-based: it
fires on variant consequence types that are mechanistically expected to
disrupt gene function (stop-gained, frameshift, splice site disruption
assessed by MaxEntScan, etc.), irrespective of any curated knowledge
about the specific variant. OncoKB’s mutation effect annotation, by
contrast, is manually curated and may draw on functional assay data,
structural evidence, or recurrence in cancer cohorts — and therefore can
assign loss-of-function status to variants (e.g. certain missense
changes) that PCGR’s consequence-based logic would not flag, or
conversely may lack an entry for a variant that PCGR’s rules would
classify as LoF. When OncoKB is active and provides a mutation effect,
that annotation takes precedence in the two-hit candidate detection
logic alongside the internal flag (either source can qualify a variant
as a LoF hit). For detailed interpretation, cross-check both columns in
the SNV/InDel table.*

**9. When running PCGR with OncoKB enabled but without specifying
`--oncokb_oncotree_code`, I see a variant assigned OncoKB level 3B
(e.g. for Bladder Cancer). If I re-run with a more specific OncoTree
code (e.g. Urethral Urothelial Carcinoma), the same variant is assigned
level 1. Why does a less specific tumor type yield a lower evidence
level?**

*Answer: OncoKB evidence levels are tied to the specific tumor type
context in which a biomarker–drug association has been validated. When
PCGR derives the OncoTree code from `--tumor_site` alone (without
`--oncokb_oncotree_code`), it maps to a broad, site-level code
(e.g. `BLCA` for Bladder Cancer). If the highest-confidence evidence in
OncoKB is annotated under a more specific subtype (e.g. Urethral
Urothelial Carcinoma), the broad code will not match that entry and
OncoKB instead returns a lower level reflecting the closest match it can
find. Providing a more specific OncoTree code via
`--oncokb_oncotree_code` allows OncoKB to resolve the exact subtype
match and return the correct, higher evidence level. For tumor types
where subtype granularity matters clinically, we therefore recommend
explicitly setting `--oncokb_oncotree_code` rather than relying on the
site-level default. A full list of OncoTree codes is available at
[oncotree.mskcc.org](https://oncotree.mskcc.org).*

**10. Does PCGR support the ESMO Scale for Clinical Actionability of
Molecular Targets (ESCAT)?**

*Answer: ESCAT is not currently implemented in PCGR, primarily since an
automated implementation is considerably non-trivial. The ESCAT
guidelines have relatively low specificity in certain areas - often
leaving room for subjective judgements (e.g. “clinically meaningful
improvement of a survival endpoint in prospective, randomised clinical
trials”). This observation is reflected in a relatively low inter-rater
institutional agreement of ESCAT-based variant rankings (see [Lebedeva
et al., Ann Oncol 2024](https://pubmed.ncbi.nlm.nih.gov/39368036/)). We
continue to monitor developments in this space (e.g. [Kordes et al.,
medRxiv
2026](https://www.medrxiv.org/content/10.64898/2026.05.16.26353390v1.full-text))
and hope to offer ESCAT support in forthcoming releases.*

**11. I notice that PCGR’s internal oncogenicity classification
sometimes assigns a variant “Likely Oncogenic” while OncoKB labels the
same variant “Oncogenic”. Why does this discrepancy occur, given that
PCGR uses OncoKB data as one of its inputs?**

*Answer: PCGR’s internal oncogenicity classification implements the
joint VICC/CGC/ClinGen guidelines, which combine evidence from multiple
sources, and weigh each source according to a defined rule-based scoring
framework. Crucially, when available, OncoKB data feeds into that
framework among several other sources (e.g. hotspot status, functional
impact, population frequency), and the resulting score may not reach the
threshold required for a definitive “Oncogenic” call even when OncoKB
itself has curated the variant as such. OncoKB’s own classifications, on
the other hand, are based on manual expert curation that can integrate
functional assay results, structural biology evidence, and broader
literature context in ways that the rule-based VICC/CGC/ClinGen
algorithm cannot fully capture algorithmically. PCGR’s classification
thus tends to be more conservative: a variant that a curator has
confidently labelled “Oncogenic” in OncoKB may accumulate evidence to
reach only “Likely Oncogenic” under the VICC/CGC/ClinGen scoring scheme.
Both annotations are reported in PCGR (columns `ONCOGENICITY` for the
internal call and `ONCOGENICITY_OKB` when OncoKB is enabled), and users
are encouraged to treat them as complementary.*

**12. In the Excel workbook, the biomarker sheets contain a
`BM_ACTIONABILITY_SUPPORT` column with values such as `tier-defining` or
`additional`. Does this mean that prognostic, diagnostic, and resistance
markers also influence variant tiering?**

*Answer: No. Variant tiering in PCGR is (as of v2.3.0) driven
exclusively by **treatment sensitivity** evidence — the `TIER` column
values `T1`, `T2`, `T3`, and `T4` reflect only drug sensitivity
biomarkers. The `BM_ACTIONABILITY_SUPPORT` column (`tier-defining` /
`additional`) describes the role of a given biomarker **within its own
evidence category**: for example, a prognostic biomarker may be the
strongest (`tier-defining`) or a corroborating (`additional`) piece of
evidence for that category. Prognostic (PP1/PB1), diagnostic (D1/D2),
and resistance (R1/R2) markers are listed in the biomarker sheets for
completeness and transparency, but they do not contribute to the
`T1`–`T4` variant rank. To identify which biomarkers actually drove the
tier assignment of a given variant, focus on rows where `TIER` is
`T1`–`T3` and the evidence type is treatment sensitivity.*
