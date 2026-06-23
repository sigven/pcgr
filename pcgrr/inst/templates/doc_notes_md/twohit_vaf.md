The **potential two-hit events** table identifies tumor suppressor genes (TSGs) subject to
allele-specific loss-of-heterozygosity (LOH) — either a deletion (nMinor=0, nMajor < ploidy
baseline) or copy-neutral LOH (nMinor=0, nMajor = ploidy baseline) — that also carry a
co-occurring somatic or germline loss-of-function (LOF) SNV/InDel. Together these may
constitute biallelic inactivation (Knudson two-hit hypothesis).

__VAF consistency assessment__

For somatic candidates, the **VAF consistency** column assesses whether the observed tumor
VAF is consistent with the variant residing on the _retained_ (major) allele of the LOH
segment. The expected VAF of a clonal mutation is derived from the purity-adjusted mixture model
([Carter et al., Nat Biotechnol, 2012](https://doi.org/10.1038/nbt.2203)):

$$\text{VAF}_{\text{expected}} = \frac{p \times n_{\text{mut}}}{p \times n_{\text{total}} + 2 \times (1 - p)}$$

where $p$ is tumor purity, $n_{\text{mut}}$ is the number of mutant allele copies on the
retained allele, and $n_{\text{total}}$ is the total segment copy number. A variant is considered consistent if the observed tumor VAF falls within a tolerance of ±0.15 of the expected VAF.

| Flag | Interpretation |
|---|---|
| **VAF_CONSISTENT** | Observed VAF ≥ expected − tolerance. The variant is likely on the retained allele — counter-evidence for the mutation residing on the deleted allele. |
| **VAF_LOW** | Observed VAF below purity-adjusted expectation. The variant may reside on the _deleted_ allele (in which case there may be cells retaining a functioning copy), or may be subclonal relative to the LOH event (mutation and LOH in different clonal lineages). Interpret with caution. |
| **VAF_UNKNOWN** | Tumor purity not provided; or VAF not available for the variant - VAF consistency cannot be assessed. |

__Interpretation notes__

- `VAF_CONSISTENT` is strong counter-evidence for the variant residing on the deleted allele
  (scenario 1 below), since a variant on the deleted allele would be expected to show VAF ≈ 0.
  It does **not**, however, rule out subclonality scenarios (2–4).
- A `VAF_LOW` flag may reflect any of the following:
  1. The LOF variant resides on the allele that was lost by LOH — a functioning WT allele
     is retained on the major allele, and this is **not** a true two-hit event.
  2. The LOH is clonal but the mutation arose later in a subclone — cells lacking the
     mutation still carry one functioning allele.
  3. The mutation is clonal but LOH is subclonal — cells without LOH carry both a mutant
     and a WT allele.
  4. Mutation and LOH arose independently in separate clonal lineages.
- For **copy-neutral LOH** the lower-bound expectation (n_mut = 1) is used, as the phase
  of the variant relative to the duplicated haplotype cannot be determined without allele
  phasing data.
- Germline candidates (pathogenic/likely-pathogenic variants from CPSR) are matched by
  gene symbol only and do not carry a VAF flag, as the relevant comparison is the germline
  genotype (heterozygous/homozygous) rather than tumor VAF.
