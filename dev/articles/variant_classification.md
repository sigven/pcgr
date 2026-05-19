# Variant classification

## Oncogenicity

PCGR evaluates somatic aberrations (SNV/InDels) for oncogenic potential
through an implementation of standard operating procedures proposed by
[ClinGen/CGC/VICC](https://pubmed.ncbi.nlm.nih.gov/35101336/). Here,
various properties of the variants and genes affected are assigned
criteria-specific scores, both negative and positive, pending on whether
the properties support an oncogenic or benign variant nature.
Criteria-specific scores are in turn aggregated towards an overall
oncogenicity score per variant.

Note that all properties/criteria provided in the SOP’s are *not*
readily implemented in PCGR, specifically the ones requiring manual
curation or expert review (i.e. experimental oncogenic variant evidence,
requiring support from *in vitro* or *in vivo* functional studies
(criteria *OS2*)). Also, the current source for existing oncogenicity
classifications (ClinVar), which is needed for implementation of
criteria *OS1* and *OM3*, is considerably limited. Considering the
current limitations, some oncogenic variants are likely to be missed and
classified with uncertain significance (VUS) by PCGR. We highlight
conventional ClinVar classifications (i.e. with respect to
pathogenicity) alongside the current oncogenicity classifications, this
may likely add value to the interpretation for uncertain cases.
Furthermore, taking into the account the nature of the current
implementation, we have adopted slightly different score thresholds for
variant classifications to those proposed originally by [Horak et al.,
2022](https://pubmed.ncbi.nlm.nih.gov/35101336/). We are working to
further improve the oncogenicity classification in PCGR, and welcome
feedback on this matter.

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

PCGR prioritizes and evaluates variants according to clinical
actionability. Currently, PCGR implements its tier classification
framework along the proposed AMP/ASCO/CAP guidelines, as outlined in [Li
et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/):

- **Tier I: Variants of strong clinical significance** - constitutes
  aberrations linked to therapeutic (drug sensitivity/response)
  biomarkers in [OncoKB](https://www.oncokb.org/),
  [CIViC](https://civicdb.org), or the [Cancer Biomarkers
  Database](https://www.cancergenomeinterpreter.org/biomarkers). These
  biomarkers are further
  - reported within the same tumor type/class as the tumor in question
    (query tumor sample), **AND**
  - Of strong clinical evidence (i.e. FDA-approved/standard-of-care
    therapies, part of guidelines, validated or discovered in late
    clinical trials ([OncoKB evidence levels
    1/2/3A](https://www.oncokb.org/therapeutic-levels), [CIViC evidence
    levels
    A/B](https://civic.readthedocs.io/en/latest/model/evidence/level.html)))
- **Tier II: Variants of potential clinical significance** - constitutes
  other aberrations linked to therapeutic (drug sensitivity/response)
  biomarkers in [OncoKB](https://www.oncokb.org/),
  [CIViC](https://civicdb.org), or [Cancer Biomarkers
  Database](https://www.cancergenomeinterpreter.org/biomarkers) that are
  either
  - Of strong clinical evidence in other tumor types/classes than the
    one specified by the user (i.e. ([OncoKB evidence level
    3B](https://www.oncokb.org/therapeutic-levels)), **OR**
  - Of weak clinical evidence (early trials, case reports etc. ([CIViC
    evidence levels
    C/D/E](https://civic.readthedocs.io/en/latest/model/evidence/level.html)),
    [OncoKB evidence level
    4](https://www.oncokb.org/therapeutic-levels)) in the same tumor
    type/class as specified by the user
- **Tier III: Variants of uncertain clinical significance (SNVs/InDels
  only)** -
  - Other coding variants, not observed at significant allele
    frequencies (gnomAD MAF \< 0.001), found in oncogenes or tumor
    suppressor genes, yet *not* linked to any known predictive,
    prognostic, or diagnostic biomarkers in the [CIViC
    database](https://civicdb.org) and the [Cancer Biomarkers
    Database](https://www.cancergenomeinterpreter.org/biomarkers)

In PCGR, we skip the classification of variants into the
AMP/ASCO/CAP-specified *Tier IV* (benign/likely benign variants), but
rather take a more cautious approach. Specifically, for SNVs/indels that
do not fall into tier I, II, or III, we classify them into *Tier V:
Other coding variants*, which includes protein-coding variants in
non-cancer related genes, as well as *Tier VI: Other non-coding
variants*, which includes synonymous variants, intronic variants, and
other variants in non-coding regions.
