CHANGELOG
---------

0.8.0 - May 20th 2019
^^^^^^^^^^^^^^^^^^^^^

Fixed
'''''

-  Bug in value box for Tier 2 variants (new line carriage) `Issue
   #73 <https://github.com/sigven/pcgr/issues/73>`__

Added
'''''

-  Upgraded VEP to v96

   -  Skipping the *–regulatory* VEP option to avoid forking issues and
      to improve speed (See `this
      issue <https://github.com/Ensembl/ensembl-vep/issues/384>`__)
   -  Added option to configure *pick-order* for choice of primary
      transcript in configuration file

-  Pre-made configuration files for each tumor type in *conf* folder
-  Possibility to append a CNA plot file (.png format) to the section of
   the report with *Somatic CNAs* `previous feature
   request <https://github.com/sigven/pcgr/issues/58>`__
-  Added possibility to input estimates of **tumor purity** and
   **ploidy**

   -  shown as value boxes in *Main results*

-  Tumor mutational burden is now compared with the distribution of TMB
   observed for TCGA’s cohorts (organized by primary site)

   -  Default target size is now 34Mb (approx. estimate from exome-wide
      calculation of protein-coding parts of GENCODE)

-  Added flexibility for variant filtering in tumor-only input callsets

   -  Added additional options to exclude likely germline variants (both
      requires the tumor VAF tag to be correctly specified in the input
      VCF)

      -  **exclude_likely_hom_germline** - removes any variant with an
         allelic fraction of 1 (100%) - very unlikely somatic event

   -  **exclude_likely_het_germline** - removes any variant with

      -  an allelic fraction between 0.4 and 0.6, and
      -  presence in dbSNP + gnomAD, and
      -  no presence as somatic event in COSMIC/TCGA

   -  Added possibility to input *PANEL-OF-NORMALS* VCF - this to
      support the many labs that have sequenced a database/pool of
      healthy controls. This set of variants are utilized in PCGR to
      improve the variant filtering when running in tumor-only mode. The
      *PANEL-OF-NORMALS* annotation work as follows:

      -  all variants in the tumor that coincide with any variant listed
         in the *PANEL-OF-NORMALS* VCF is appended with a
         **PANEL_OF_NORMALS** flag in the query VCF with tumor variants.

   -  If configuration parameter **exclude_pon** is set to True in
      **tumor_only** runs, all variants with a **PANEL_OF_NORMALS** flag
      are filtered/excluded

-  For tumor-only runs, added an `UpSet
   plot <https://github.com/hms-dbmi/UpSetR#Demo>`__ showing how
   different filtering sources (gnomAD, 1KG Project, panel-of-normals
   etc) contribute in the germline filtering procedure
-  Variants in *Tier 3 / Tier 4 / Noncoding* are now sorted (and
   color-coded) according to the target (gene) association score to the
   cancer phenotype, as provided by the `OpenTargets
   Platform <https://docs.targetvalidation.org/getting-started/scoring>`__
-  Added annotation of TCGA’s ten oncogenic signaling pathways
-  Added *EXONIC_STATUS* annotation tag (VCF and TSV)

   -  *exonic* denotes all protein-altering AND cannonical splicesite
      altering AND synonymous variants, *nonexonic* denotes the
      complement

-  Added *CODING_STATUS* annotation tag (VCF and TSV)

   -  *coding* denotes all protein-altering AND cannonical splicesite
      altering, *noncoding* denotes the complement

-  Added *SYMBOL_ENTREZ* annotation tag (VCF)

   -  Official gene symbol from NCBI EntreZ (SYMBOL provided by VEP can
      sometimes be non-official/alias (i.e. for GENCODE v19/grch37))

-  Added *SIMPLEREPEATS_HIT* annotation tag (VCF and TSV)

   -  Variant overlaps UCSC *simpleRepeat* sequence repeat track - used
      for MSI prediction

-  Added *WINMASKER_HIT* annotation tag (VCF and TSV)

   -  Variant overlaps UCSC *windowmaskerSdust* sequence repeat track -
      used for MSI prediction

-  Added *PUTATIVE_DRIVER_MUTATION* annotation tag (VCF and TSV)

   -  Putative cancer driver mutation discovered by multiple approaches
      from 9,423 tumor exomes in TCGA. Format:
      symbol:hgvsp:ensembl_transcript_id:discovery_approaches

-  Added *OPENTARGETS_DISEASE_ASSOCS* annotation tag (VCF and TSV)

   -  Associations between protein targets and disease based on multiple
      lines of evidence (mutations,affected pathways,GWAS, literature
      etc). Format: CUI:EFO_ID:IS_DIRECT:OVERALL_SCORE

-  Added *OPENTARGETS_TRACTABILITY_COMPOUND* annotation tag (VCF and
   TSV)

   -  Confidence for the existence of a modulator (small molecule) that
      interacts with the target (protein) to elicit a desired biological
      effect

-  Added *OPENTARGTES_TRACTABILITY_ANTIBODY* annotation tag (VCF and
   TSV)

   -  Confidence for the existence of a modulator (antibody) that
      interacts with the target (protein) to elicit a desired biological
      effect

-  Added *CLINVAR_REVIEW_STATUS_STARS* annotation tag

   -  Rating of the ClinVar variant (0-4 stars) with respect to level of
      review

Changed
'''''''

-  Moved from `IntoGen’s driver mutation
   resource <https://www.intogen.org/>`__ to `TCGA’s putative driver
   mutation list <https://doi.org/10.1016/j.cell.2018.02.060>`__ in
   display of driver mutation status
-  Moved option for vcf_validation from configuration file to run script
   (``--no_vcf_validate``)

Removed
'''''''

-  Original tier model ‘pcgr’

0.7.0 - Nov 27th 2018
^^^^^^^^^^^^^^^^^^^^^

.. _fixed-1:

Fixed
'''''

-  Bug in assignment of variants to tier1/tier2 `Issue
   #61 <https://github.com/sigven/pcgr/issues/61>`__
-  Missing config option for *maf_gnomad_asj* in TOML file (also setting
   operator to ``<=``) `Issue
   #60 <https://github.com/sigven/pcgr/issues/60>`__
-  Bug in new CancerMine oncogene/tumor suppressor annotation `Issue
   #53 <https://github.com/sigven/pcgr/issues/53>`__
-  vcfanno fix for empty Description (upgrade to vcfanno v0.3.1 `Issue
   #49 <https://github.com/sigven/pcgr/issues/49>`__)
-  Bug in message showing too few variants for MSI prediction, `Issue
   #55 <https://github.com/sigven/pcgr/issues/55>`__
-  Bug in appending of custom VCF tags

   -  Still unsolved: how to disambiguate identical FORMAT and INFO tags
      in vcf2tsv

-  Bug in SCNA value box display for multiple copy number hits (`Issue
   #47 <https://github.com/sigven/pcgr/issues/47>`__)
-  Bug in vcf2tsv (handling INFO tags encoded with ‘Type = String’,
   `Issue #39 <https://github.com/sigven/pcgr/issues/39>`__)
-  Bug in search of UniProt functional features (BED feature regions
   spanning exons are now handled)
-  Stripped off HTML elements (TCGA_FREQUENCY, DBSNP) in TSV output
-  Some effect predictions from dbNSFP were not properly parsed
   (e.g. multiple prediction entries from multiple transcript isoforms),
   these should now be retrieved correctly
-  Removed ‘COSM’ prefix in COSMIC mutation links
-  Bug in retrieval of splice site predictions from dbscSNV

.. _added-1:

Added
'''''

-  Possibility to run PCGR in a non-Docker environment (e.g. using the
   *–no-docker* option). Thanks to an excellent contribution by `Vlad
   Saveliev <https://github.com/vladsaveliev>`__, `Issue
   #35 <https://github.com/sigven/pcgr/issues/35>`__

   -  Added possibility to add docker user-id

-  Possibility for MAF file output (converted with vcf2maf), must be
   configured by the user in the TOML file (i.e. *vcf2maf = true*,
   `Issue #17 <https://github.com/sigven/pcgr/issues/17>`__)
-  Possibility for adding custom VCF INFO tags to PCGR output files
   (JSON/TSV), must be configured by the user in the TOML file (i.e.
   *custom_tags*)
-  Added MUTATION_HOTSPOT_CANCERTYPE in data tables (i.e. listing tumor
   types in which hotspot mutations have been found)
-  Included the ‘rs’ prefix for dbSNP identifiers (HTML and TSV output)
-  Individual entries/columns for variant effect predictions:

   -  Individual algorithms: SIFT_DBNSFP, M_CAP_DBNSFP, MUTPRED_DBNSFP,
      MUTATIONTASTER_DBNSFP, MUTATIONASSESSOR_DBNSFP, FATHMM_DBNSFP,
      FATHMM_MKL_DBNSFP, PROVEAN_DBNSFP
   -  Ensemble predictions (META_LR_DBNSFP), dbscSNV splice site
      predictions (SPLICE_SITE_RF_DBNSFP, SPLICE_SITE_ADA_DBNSFP)

-  Upgraded samtools to v1.9 (makes vcf2maf work properly)
-  Added Ensembl gene/transcript id and corresponding RefSeq mRNA id to
   TSV/JSON
-  Added for future implementation:

   -  SeqKat + karyoploteR for exploration of *kataegis/hypermutation*
   -  CELLector - genomics-guided selection of cancer cell lines

-  Upgraded VEP to v94

.. _changed-1:

Changed
'''''''

-  Changed CANCER_MUTATION_HOTSPOT to MUTATION_HOTSPOT
-  Moved from `TSGene 2.0 <https://bioinfo.uth.edu/TSGene/>`__ to
   `CancerMine <https://zenodo.org/record/1336650#.W9QMdRMzaL4>`__ for
   annotation of tumor suppressor genes and proto-oncogenes

   -  A minimum of n=3 citations were required to include
      literatured-mined tumor suppressor genes and proto-oncogenes from
      CancerMine

0.6.2.1 - May 14th 2018
^^^^^^^^^^^^^^^^^^^^^^^

.. _fixed-2:

Fixed
'''''

-  Bug in copy number annotation (broad/focal)

0.6.2 - May 9th 2018
^^^^^^^^^^^^^^^^^^^^

.. _fixed-3:

Fixed
'''''

-  Bug in copy number segment display (missing variable initalization,
   `Issue #34 <https://github.com/sigven/pcgr/issues/34>`__))
-  Typo in gnomAD filter statistic (fraction, `Issue
   #31 <https://github.com/sigven/pcgr/issues/31>`__)
-  Bug in mutational signature analysis for grch38 (forgot to pass
   BSgenome object, `Issue
   #27 <https://github.com/sigven/pcgr/issues/27>`__)
-  Missing proper ASCII-encoding in vcf2tsv conversion, `Issue
   # <https://github.com/sigven/pcgr/issues/35>`__
-  Removed ‘Noncoding mutations’ section when no input VCF is present
-  Bug in annotation of copy number event type (focal/broad)
-  Bug in copy number annotation (missing protein-coding transcripts)
-  Updated MSI prediction (variable importance, performance measures)

.. _added-2:

Added
'''''

-  Genome assembly is appended to every output file
-  Issue warning for copy number segment that goes beyond chromosomal
   lengths of specified assembly (segments will be skipped)
-  Added missing subtypes for ‘Skin_Cancer_NOS’ in the cancer phenotype
   dataset

0.6.1 - May 2nd 2018
^^^^^^^^^^^^^^^^^^^^

.. _fixed-4:

Fixed
'''''

-  Bug in tier assignment ‘pcgr_acmg’ (case for no variants in
   tier1,2,3)
-  Bug in tier assignment ‘pcgr_acmg’ (no tumor type specified, evidence
   items with weak support detected)
-  Bug: duplicated variants in ‘Tier 3’ resulting from genes encoded
   with dual roles as tumor suppressor genes/oncogenes
-  Bug: duplicated variants in ‘Tier 1/Noncoding variants’ resulting
   from rare cases of noncoding variants occurring in Tier 1 (synonymous
   variants with biomarker role)

0.6.0 - April 25th 2018
^^^^^^^^^^^^^^^^^^^^^^^

.. _added-3:

Added
'''''

-  New argument in pcgr.py

   -  *assembly* (grch37/grch38)

-  New option in pcgr.py

   -  *–basic* - run comprehensive VCF annotation only, skip report
      generation and additional analyses

-  New sections in HTML report

   -  *Settings and annotation sources* - now also listing key PCGR
      configuration settings
   -  *Main findings* - Six value boxes indicating the main findings of
      clinical relevance

-  New configuration options

   -  [tier_model](string) - choice between *pcgr_acmg* and *pcgr*
   -  [mutational_burden] - set TMB tertile limits

      -  *tmb_low_limit (float)*
      -  *tmb_intermediate_limit (float)*

   -  [tumor_type] - choose between 34 tumor types/classes:

      -  *Adrenal_Gland_Cancer_NOS (logical)*
      -  *Ampullary_Carcinoma_NOS (logical)*
      -  *Biliary_Tract_Cancer_NOS (logical)*
      -  *Bladder_Urinary_Tract_Cancer_NOS (logical)*
      -  *Blood_Cancer_NOS (logical)*
      -  *Bone_Cancer_NOS (logical)*
      -  *Breast_Cancer_NOS (logical)*
      -  *CNS_Brain_Cancer_NOS (logical)*
      -  *Colorectal_Cancer_NOS (logical)*
      -  *Cervical_Cancer_NOS (logical)*
      -  *Esophageal_Stomach_Cancer_NOS (logical)*
      -  *Head_And_Neck_Cancer_NOS (logical)*
      -  *Hereditary_Cancer_NOS (logical)*
      -  *Kidney_Cancer_NOS (logical)*
      -  *Leukemia_NOS (logical)*
      -  *Liver_Cancer_NOS (logical)*
      -  *Lung_Cancer_NOS (logical)*
      -  *Lymphoma_Hodgkin_NOS (logical)*
      -  *Lymphoma_Non_Hodgkin_NOS (logical)*
      -  *Ovarian_Fallopian_Tube_Cancer_NOS (logical)*
      -  *Pancreatic_Cancer_NOS (logical)*
      -  *Penile_Cancer_NOS (logical)*
      -  *Peripheral_Nervous_System_Cancer_NOS (logical)*
      -  *Peritoneal_Cancer_NOS (logical)*
      -  *Pleural_Cancer_NOS (logical)*
      -  *Prostate_Cancer_NOS (logical)*
      -  *Skin_Cancer_NOS (logical)*
      -  *Soft_Tissue_Cancer_NOS (logical)*
      -  *Stomach_Cancer_NOS (logical)*
      -  *Testicular_Cancer_NOS (logical)*
      -  *Thymic_Cancer_NOS (logical)*
      -  *Thyroid_Cancer_NOS (logical)*
      -  *Uterine_Cancer_NOS (logical)*
      -  *Vulvar_Vaginal_Cancer_NOS (logical)*

   -  [mutational_signatures]

      -  *mutsignatures_cutoff (float)* - discard any signature
         contributions with a weight less than the cutoff

   -  [cna]

      -  *transcript_cna_overlap (float)* - minimum percent overlap
         between copy number segment and transcripts (average) for tumor
         suppressor gene/proto-oncogene to be reported

   -  [allelic_support]

      -  If input VCF has correctly formatted depth/allelic fraction as
         INFO tags, users can add thresholds on depth/support that are
         applied prior to report generation

         -  *tumor_dp_min (integer)* - minimum sequencing depth for
            variant in tumor sample
         -  *tumor_af_min (float)* - minimum allelic fraction for
            variant in tumor sample
         -  *normal_dp_min (integer)* - minimum sequencing depth for
            variant in normal sample
         -  *normal_af_max (float)* - maximum allelic fraction for
            variant in normal sample

   -  [visual]

      -  *report_theme (string)* - visual theme of report (Bootstrap)

   -  [other]

      -  *vcf_validation (logical)* - keep/skip VCF validation by
         `vcf-validator <https://github.com/EBIvariation/vcf-validator>`__

-  New output file - JSON output of HTML report content
-  New INFO tags of PCGR-annotated VCF

   -  *CANCER_PREDISPOSITION*
   -  *PFAM_DOMAIN*
   -  *TCGA_FREQUENCY*
   -  *TCGA_PANCANCER_COUNT*
   -  *ICGC_PCAWG_OCCURRENCE*
   -  *ICGC_PCAWG_AFFECTED_DONORS*
   -  *CLINVAR_MEDGEN_CUI*

-  New column entries in annotated SNV/InDel TSV file:

   -  *CANCER_PREDISPOSITION*
   -  *ICGC_PCAWG_OCCURRENCE*
   -  *TCGA_FREQUENCY*

-  New column in CNA output

   -  *TRANSCRIPTS* - aberration-overlapping transcripts (Ensembl
      transcript IDs)
   -  *MEAN_TRANSCRIPT_CNA_OVERLAP* - Mean overlap (%) betweeen gene
      transcripts and aberration segment

.. _removed-1:

Removed
'''''''

-  Elements of databundle (now annotated directly through VEP):

   -  dbsnp
   -  gnomad/exac
   -  1000G project

-  INFO tags of PCGR-annotated VCF

   -  *DBSNPBUILDID*
   -  *DBSNP_VALIDATION*
   -  *DBSNP_SUBMISSIONS*
   -  *DBSNP_MAPPINGSTATUS*
   -  *GWAS_CATALOG_PMID*
   -  *GWAS_CATALOG_TRAIT_URI*
   -  *DOCM_DISEASE*

-  Output files

   -  TSV files with mutational signature results and biomarkers (i.e.
      *sample_id.pcgr.snvs_indels.biomarkers.tsv* and
      *sample_id.pcgr.mutational_signatures.tsv*)

      -  Data can still be retrieved - now from the JSON dump

   -  MAF file

      -  The previous MAF output was generated in a custom fashion, a
         more accurate MAF output based on
         https://github.com/mskcc/vcf2maf will be incorporated in the
         next release

.. _changed-2:

Changed
'''''''

-  HTML report sections

   -  *Tier statistics* and *Variant statistics* are now grouped into
      the section *Tier and variant statistics*
   -  *Tier 5* is now *Noncoding mutations* (i.e. not considered a tier
      per se)
   -  Sliders for allelic fraction in the *Global variant browser* are
      now fixed from 0 to 1 (0.05 intervals)
