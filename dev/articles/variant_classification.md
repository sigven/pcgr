# Variant classification

## Oncogenicity

Somatic aberrations (SNV/InDels) are evaluated for oncogenicity through
an implementation of standard operating procedures proposed by
ClinGen/CGC/VICC (Horak et al. 2022). Here, various properties of the
variants and genes affected are assigned specific scores according to
several criteria, both negative and positive, pending on whether the
properties support an oncogenic or benign variant type. These scores are
in turn aggregated towards an overall oncogenicity score.

Note that all properties/criteria provided in the SOP’s are *not*
readily implemented in PCGR, specifically the ones requiring manual
curation or expert review (i.e. experimental oncogenic variant evidence,
requiring support from *in vitro* or *in vivo* functional studies
(criteria *OS2*)). Also, the current source for existing oncogenicity
classifications (ClinVar), which is needed for implementation of
criteria *OS1*, *OM3*, is also considerably limited. Considering the
current limitations, some oncogenic variants are likely to be missed and
classified with uncertain significance (VUS) by PCGR. We highlight
conventional ClinVar classifications (i.e. with respect to
pathogenicity) alongside the current oncogenicity classifications, this
may add value to the interpretation for uncertain cases. Furthermore,
taking into the account the nature of the current implementation, we
have adopted slightly different score thresholds for variant
classifications to those proposed originally by (Horak et al. 2022). We
are working to further improve the oncogenicity classification in PCGR,
and welcome feedback on this matter.

Note also that for somatic copy number aberrations, we showcase
potential oncogenic events as **proto-oncogenes subject to
amplifications** (where level of amplification is configurable by the
user), as well as **tumor suppressor genes subject to homozygous
deletions**.

The following criteria/codes are currently used for variant oncogenicity
classification in PCGR (key resources/tools used for implementation
indicated in parentheses):

- *ONCG_OVS1* - Null variant - predicted as LoF - in bona fide tumor
  suppressor gene (VEP;CGC;CancerMine)
- *ONCG_OVS1_A* - Null variant - annotated by OncoKB as
  Loss-of-function - in bona fide tumor suppressor gene (OncoKB)
- *ONCG_OVS1_B* - Null variant - annotated by OncoKB as Likely
  Loss-of-function - in bona fide tumor suppressor gene (OncoKB)
- *ONCG_OS1* - Same amino acid change as previously established
  oncogenic variant - regardless of nucleotide change (ClinVar)
- *ONCG_OS2_A* - Well established in vitro/in vivo functional studies
  (OncoKB-curated) show oncogenic effect of variant (OncoKB)
- *ONCG_OS2_B* - Well established in vitro/in vivo functional studies
  (OncoKB-curated) show likely oncogenic effect of variant (OncoKB)
- *ONCG_OS3* - Located in a mutation hotspot with \>= 50 samples with
  variant at AA position, \>= 10 samples with same AA change
  (cancerhotspots.org)
- *ONCG_OM1* - Presumably critical site of functional domain (CIViC)
- *ONCG_OM2* - Protein length changes from in-frame dels/ins in known
  oncogene/tumor suppressor genes or stop-loss variants in a tumor
  suppressor gene (VEP;CGC;CancerMine)
- *ONCG_OM3* - Missense variant at an amino acid residue where a
  different missense variant determined to be oncogenic (using this
  standard) has been documented (ClinVar)
- *ONCG_OM4* - Located in a mutation hotspot with \< 50 samples with
  variant at AA position, \>= 10 samples with same AA change
  (cancerhotspots.org)
- *ONCG_OP1* - Multiple lines of computational evidence support of a
  damaging variant effect on the gene or gene product (dbNSFP)
- *ONCG_OP3* - Located in a mutation hotspot with \< 10 samples with the
  same amino acid change (cancerhotspots.org)
- *ONCG_OP4* - Absent from controls (gnomAD) / very low MAF (any five
  major gnomAD subpopulations) (gnomAD)
- *ONCG_SBVS1* - Very high MAF (any five major gnomAD subpopulations)
  (gnomAD)
- *ONCG_SBS1* - High MAF (any five major gnomAD subpopulations) (gnomAD)
- *ONCG_SBS2_A* - Well established in vitro/in vivo functional studies
  (OncoKB-curated) show neutral effect of variant (OncoKB)
- *ONCG_SBS2_B* - Well established in vitro/in vivo functional studies
  (OncoKB-curated) show likely neutral effect of variant (OncoKB)
- *ONCG_SBP1* - Multiple lines of computational evidence support a
  benign variant effect on the gene or gene product (dbNSFP)
- *ONCG_SBP2* - Silent and intronic changes outside of the consensus
  splice site (VEP)

## Actionability

Clinical actionability assessment of somatic SNVs/InDels, gene copy
number aberrations, and RNA fusions found in the tumor sample implements
recommendation guidelines by AMP/ASCO/CAP (Li et al. 2017).
Specifically, different levels of actionability are implemented in the
following manner:

- **Tier I: Variants of strong clinical significance** — aberrations
  linked to predictive, prognostic, or diagnostic biomarkers in the
  [CIViC database](https://civicdb.org) and the [Cancer Biomarkers
  Database](https://www.cancergenomeinterpreter.org/2021/biomarkers)
  (and optionally [OncoKB](https://oncokb.org)) that show
  - Strong clinical evidence (approved therapies, guidelines, late-phase
    trials; [CIViC evidence levels
    A/B](https://civic.readthedocs.io/en/latest/model/evidence/level.html);
    [OncoKB therapeutic levels 1, 2, 3A,
    3B](https://www.oncokb.org/therapeutic-levels)) **in the same tumor
    type** as specified by the user, **OR**
  - Strong clinical evidence in a **pan-cancer** (`Any`) context
    (treated as equivalent to a matching-site strong item)
- **Tier II: Variants of potential clinical significance** — aberrations
  linked to predictive, prognostic, or diagnostic biomarkers in the
  [CIViC database](https://civicdb.org) and the [Cancer Biomarkers
  Database](https://www.cancergenomeinterpreter.org/2021/biomarkers)
  (and optionally [OncoKB](https://oncokb.org)) that show either
  - Strong clinical evidence ([CIViC
    A/B](https://civic.readthedocs.io/en/latest/model/evidence/level.html);
    [OncoKB therapeutic levels 1, 2, 3A,
    3B](https://www.oncokb.org/therapeutic-levels)) in a
    **non-matching** tumor type, **OR**
  - Weak clinical evidence (early trials, case reports; [CIViC evidence
    levels
    C/D/E](https://civic.readthedocs.io/en/latest/model/evidence/level.html);
    [OncoKB therapeutic level
    4](https://www.oncokb.org/therapeutic-levels)) in the **same tumor
    type** as specified by the user
- **Tier III: Variants of unknown/uncertain clinical significance** —
  Tier I and II assignment is alteration-type agnostic (driven solely by
  biomarker evidence strength and tumor-type matching). Tier III,
  however, is assigned through alteration-type specific routes:
  1.  **Biomarker-linked, weak evidence only** (all alteration types) —
      variants linked to known biomarkers but only through weak clinical
      evidence ([CIViC
      C/D/E](https://civic.readthedocs.io/en/latest/model/evidence/level.html);
      [OncoKB therapeutic level
      4](https://www.oncokb.org/therapeutic-levels)) in a non-matching
      tumor type or pan-cancer (`Any`) context
  2.  **No biomarker link** — alteration-type specific criteria:
      - *SNVs/InDels*: coding variants in oncogenes or tumor suppressor
        genes with a population allele frequency below 0.001 (gnomAD),
        not linked to any known biomarker in underlying databases
      - *Copy number aberrations*: amplified oncogenes, or deleted
        (homozygous, hemizygous, or heterozygous) tumor suppressor
        genes, yet not linked to any known biomarker in underlying
        databases
      - *RNA fusions*: fusions where at least one partner gene (5’ or
        3’) is an oncogene and where there is a record in the Mitelman
        database, yet not linked to any known biomarker in underlying
        databases

In PCGR, we skip the classification of variants into the
AMP/ASCO/CAP-specified *Tier IV* (benign/likely benign variants), but
rather take a more cautious approach. Specifically, for SNVs/indels that
do not fall into tier I, II, or III, we classify them into *Tier V:
Other coding variants*, which includes protein-coding variants in
non-cancer related genes, as well as *Tier VI: Other non-coding
variants*, which includes synonymous variants, intronic variants, and
other variants in non-coding regions.

### Tier I/II — alteration-type agnostic biomarker matching

Tier I and II assignment applies identically to all somatic alteration
types (SNVs/InDels, copy number aberrations, and RNA fusions).
Assignment depends on whether the user has specified a primary tumor
site (tumor-type specific mode) or the generic `Any` site (tumor-type
agnostic mode):

**Tumor-type specific** (a primary tumor site is specified, e.g. *Lung*,
*Breast*):

| Biomarker site | Evidence strength | Tier |
|----|----|----|
| Matching site **or** pan-cancer (`Any`) | Strong (CIViC A/B; [OncoKB therapeutic 1/2/3A/3B](https://www.oncokb.org/therapeutic-levels)) | **I** |
| Non-matching tumor-specific site | Strong (CIViC A/B; [OncoKB therapeutic 1/2/3A/3B](https://www.oncokb.org/therapeutic-levels)) | **II** |
| Matching site | Weak (CIViC C/D/E; OncoKB therapeutic 4) | **II** |
| Non-matching tumor-specific site **or** pan-cancer (`Any`) | Weak (CIViC C/D/E; OncoKB therapeutic 4) | **III** |

**Tumor-type agnostic** (primary site set to `Any`):

| Biomarker site | Evidence strength | Tier |
|----|----|----|
| Pan-cancer (`Any`) | Strong (CIViC A/B; [OncoKB therapeutic 1/2/3A/3B](https://www.oncokb.org/therapeutic-levels)) | **I** |
| Non-pan-cancer (any tumor-specific site) | Strong (CIViC A/B; [OncoKB therapeutic 1/2/3A/3B](https://www.oncokb.org/therapeutic-levels)) | **II** |
| Pan-cancer (`Any`) **or** any tumor-specific site | Weak (CIViC C/D/E; OncoKB therapeutic 4) | **III** |

Pan-cancer strong evidence is treated as equivalent to a matching-site
strong item and always promotes a variant to Tier I, regardless of query
mode. Pan-cancer weak evidence contributes to Tier III in both modes.

### Tier III — alteration-type specific routes

Unlike Tier I/II, Tier III assignment is alteration-type specific. Two
distinct routes lead to Tier III and it is important to distinguish
between them:

**Route 1 — Biomarker-linked, weak evidence only** (all alteration
types)

Variants that *are* linked to known biomarkers but only through weak
clinical evidence (CIViC C/D/E; [OncoKB therapeutic level
4](https://www.oncokb.org/therapeutic-levels)) in non-matching tumor
types or pan-cancer contexts. These variants have some clinical signal
but insufficient evidence to reach Tier II. This route applies equally
to SNVs/InDels, copy number aberrations, and RNA fusions.

**Route 2 — No biomarker link** (alteration-type specific criteria)

Variants not linked to any known biomarker in the underlying databases,
but flagged as potentially significant based on gene/alteration
properties:

| Alteration type | Tier III criterion |
|----|----|
| SNVs/InDels | Coding variant in an oncogene or tumor suppressor gene and population allele frequency \< 0.001 (gnomAD) |
| Copy number aberrations | Amplified oncogene, or deleted (homozygous, hemizygous, or heterozygous) tumor suppressor gene |
| RNA fusions | At least one fusion partner gene (5’ or 3’) is an oncogene and at least one record for the fusion event in the Mitelman database |

### Evidence items: tier-defining vs. additional

Each biomarker evidence item associated with a variant is classified as
either *tier-defining* or *additional* support:

- **Tier-defining** — the item directly justifies the variant’s assigned
  tier (e.g. a matching-site strong item on a Tier I variant, or a weak
  item on a Tier III variant)
- **Additional** — the item provides supporting context but did not
  drive the tier assignment (e.g. a non-matching-site strong item on a
  Tier I variant, or pan-cancer weak evidence on a Tier I or II variant)

## References

Horak, Peter, Malachi Griffith, Arpad M Danos, et al. 2022. “Standards
for the Classification of Pathogenicity of Somatic Variants in Cancer
(Oncogenicity): Joint Recommendations of Clinical Genome Resource
(ClinGen), Cancer Genomics Consortium (CGC), and Variant Interpretation
for Cancer Consortium (VICC).” *Genet. Med.*, January.
<http://dx.doi.org/10.1016/j.gim.2022.01.001>.

Li, Marilyn M, Michael Datto, Eric J Duncavage, et al. 2017. “Standards
and Guidelines for the Interpretation and Reporting of Sequence Variants
in Cancer: A Joint Consensus Recommendation of the Association for
Molecular Pathology, American Society of Clinical Oncology, and College
of American Pathologists.” *J. Mol. Diagn.* 19 (1): 4–23.
<https://doi.org/10.1016/j.jmoldx.2016.10.002>.
