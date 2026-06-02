Clinical actionability assessment of somatic SNVs/InDels, gene copy number aberrations, and
RNA fusions found in the tumor sample implements recommendation guidelines by AMP/ASCO/CAP [@Li2017-ew].
Specifically, different levels  of actionability are implemented in the following manner:

- **Tier I: Variants of strong clinical significance** - constitutes aberrations
  linked to predictive, prognostic, or diagnostic biomarkers in the
  [CIViC database]({civic_url}) and the [Cancer Biomarkers Database]({cgi_url})
  (and optionally [OncoKB]({oncokb_url})) that are
    - Found within the same tumor type/class as specified by the user, **AND**
    - Of strong clinical evidence (i.e. approved therapies, part of guidelines,
      validated or discovered in late clinical trials
      ([CIViC evidence levels A/B]({civic_docs_url}model/evidence/level.html);
      OncoKB therapeutic levels 1, 2, 3A, 3B, R1, R2; OncoKB diagnostic/prognostic
      level Dx1/Px1))

- **Tier II: Variants of potential clinical significance** - constitutes other
  aberrations linked to predictive, prognostic, or diagnostic biomarkers in the
  [CIViC database]({civic_url}) and the [Cancer Biomarkers Database]({cgi_url})
  (and optionally [OncoKB]({oncokb_url})) that are either
    - Of strong clinical evidence in other tumor types/classes than the one 
      specified by the user, **OR**
    - Of weak clinical evidence (early trials, case reports etc.
      ([CIViC evidence levels C/D/E]({civic_docs_url}model/evidence/level.html);
      OncoKB therapeutic level 4)
      in the same tumor type/class as specified by the user

- **Tier III: Variants of uncertain clinical significance (SNVs/InDels only)**
    - Other coding variants, not observed at significant allele frequencies
      (gnomAD MAF < 0.001), found in oncogenes or tumor suppressor genes,
      yet _not_ linked to any known predictive, prognostic, or diagnostic
      biomarkers in the [CIViC database]({civic_url}), the
      [Cancer Biomarkers Database]({cgi_url}), or [OncoKB]({oncokb_url})
      (if OncoKB integration is enabled)

In PCGR, we skip the classification of variants into the AMP/ASCO/CAP-specified 
*Tier IV* (benign/likely benign variants), but rather take a more cautious approach. 
Specifically, for SNVs/indels that do not fall into tier I, II, or III, we classify them into
*Tier V: Other coding variants*, which includes protein-coding variants in non-cancer 
related genes, as well as *Tier VI: Other non-coding variants*, which includes 
synonymous variants, intronic variants, and other variants in non-coding regions.
