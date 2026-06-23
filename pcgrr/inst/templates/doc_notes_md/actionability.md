Clinical actionability assessment of somatic SNVs/InDels, gene copy number aberrations, and
RNA fusions found in the tumor sample implements recommendation guidelines by AMP/ASCO/CAP [@Li2017-ew].
Specifically, different levels of actionability are implemented in the following manner:

- **Tier I: Variants of strong clinical significance** — aberrations linked to
  predictive, prognostic, or diagnostic biomarkers in the
  [CIViC database]({civic_url}) and the [Cancer Biomarkers Database]({cgi_url})
  (and optionally [OncoKB]({oncokb_url})) that show
    - Strong clinical evidence (approved therapies, guidelines, late-phase trials;
      [CIViC evidence levels A/B]({civic_docs_url}model/evidence/level.html);
      [OncoKB therapeutic levels 1, 2, 3A, 3B]({oncokb_levels_url})) **in the same tumor type** as
      specified by the user, **OR**
    - Strong clinical evidence in a **pan-cancer** (`Any`) context (treated as
      equivalent to a matching-site strong item)

- **Tier II: Variants of potential clinical significance** — aberrations linked
  to predictive, prognostic, or diagnostic biomarkers in the
  [CIViC database]({civic_url}) and the [Cancer Biomarkers Database]({cgi_url})
  (and optionally [OncoKB]({oncokb_url})) that show either
    - Strong clinical evidence ([CIViC A/B]({civic_docs_url}model/evidence/level.html);
      [OncoKB therapeutic levels 1, 2, 3A, 3B]({oncokb_levels_url})) in a **non-matching** tumor type, **OR**
    - Weak clinical evidence (early trials, case reports;
      [CIViC evidence levels C/D/E]({civic_docs_url}model/evidence/level.html);
      [OncoKB therapeutic level 4]({oncokb_levels_url})) in the **same tumor type** as specified by the user

- **Tier III: Variants of unknown/uncertain clinical significance** — Tier I and II
  assignment is alteration-type agnostic (driven solely by biomarker evidence
  strength and tumor-type matching). Tier III, however, is assigned through
  alteration-type specific routes:
    1. **Biomarker-linked, weak evidence only** (all alteration types) — variants
       linked to known biomarkers but only through weak clinical evidence
       ([CIViC C/D/E]({civic_docs_url}model/evidence/level.html);
       [OncoKB therapeutic level 4]({oncokb_levels_url})) in a non-matching tumor type or pan-cancer
       (`Any`) context
    2. **No biomarker link** — alteration-type specific criteria:
       - _SNVs/InDels_: coding variants in oncogenes or tumor suppressor genes
         with a population allele frequency below 0.001 (gnomAD), not linked to
         any known biomarker in underlying databases
       - _Copy number aberrations_: amplified oncogenes, or deleted (homozygous,
         hemizygous, or heterozygous) tumor suppressor genes, yet not linked to any
         known biomarker in underlying databases
       - _RNA fusions_: fusions where at least one partner gene (5' or 3') is an
         oncogene and where there is a record in the Mitelman database, yet 
         not linked to any known biomarker in underlying databases

In PCGR, we skip the classification of variants into the AMP/ASCO/CAP-specified
*Tier IV* (benign/likely benign variants), but rather take a more cautious approach.
Specifically, for SNVs/indels that do not fall into tier I, II, or III, we classify them into
*Tier V: Other coding variants*, which includes protein-coding variants in non-cancer
related genes, as well as *Tier VI: Other non-coding variants*, which includes
synonymous variants, intronic variants, and other variants in non-coding regions.
