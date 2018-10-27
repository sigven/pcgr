CHANGELOG
---------

dev/unreleased - Oct 27th 2018
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fixed
'''''

-  Bug in appending of custom VCF tags

   -  Still unsolved: how to disambiguate identical FORMAT and INFO tags
      in vcf2tsv

-  Bug in SCNA value box display for multiple copy number hits (paste
   error)
-  Bug in vcf2tsv (handling INFO tags encoded with 'Type = String')
-  Bug in search of UniProt functional features (BED feature regions
   spanning exons are now handled)
-  Stripped off HTML elements (TCGA\_FREQUENCY, DBSNP) in TSV output
-  Some effect predictions from dbNSFP were not properly parsed (e.g.
   multiple prediction entries from multiple transcript isoforms), these
   should now be retrieved correctly
-  Removed 'COSM' prefix in COSMIC mutation links
-  Bug in retrieval of splice site predictions from dbscSNV

Added
'''''

-  Possibility to run PCGR in a non-Docker environment (e.g. using the
   *--no-docker* option). Thanks to an excellent contribution by `Vlad
   Saveliev <https://github.com/vladsaveliev>`__

   -  Added possibility to add docker user-id

-  Possibility for MAF file output (converted with vcf2maf), must be
   configured by the user in the TOML file (i.e. *vcf2maf = true*)
-  Possibility for adding custom VCF INFO tags to PCGR output files
   (JSON/TSV), must be configured by the user in the TOML file (i.e.
   *custom\_tags*)
-  Addded MUTATION\_HOTSPOT\_CANCERTYPES in data tables (i.e. listing
   tumor types in which hotspot mutations have been found)
-  Included the 'rs' prefix for dbSNP identifiers (HTML and TSV output)
-  Individual entries/columns for variant effect predictions:

   -  Individual algorithms: SIFT\_DBNSFP, M\_CAP\_DBNSFP,
      MUTPRED\_DBNSFP, MUTATIONTASTER\_DBNSFP, MUTATIONASSESSOR\_DBNSFP,
      FATHMM\_DBNSFP, FATHMM\_MKL\_DBNSFP, PROVEAN\_DBNSFP
   -  Ensemble predictions (META\_LR\_DBNSFP), dbscSNV splice site
      predictions (SPLICE\_SITE\_RF\_DBNSFP, SPLICE\_SITE\_ADA\_DBNSFP)

-  Upgraded samtools to v1.9 (makes vcf2maf work properly)
-  Added Ensembl gene/transcript id and corresponding RefSeq mRNA id to
   TSV/JSON
-  Added for future implementation:

   -  SeqKat + karyoploteR for exploration of *kataegis/hypermutation*
   -  CELLector - genomics-guided selection of cancer cell lines

-  Upgraded VEP to v94

Changed
'''''''

-  Changed CANCER\_MUTATION\_HOTSPOT to MUTATION\_HOTSPOT
-  Moved from `TSGene 2.0 <https://bioinfo.uth.edu/TSGene/>`__ to
   `CancerMine <https://zenodo.org/record/1336650#.W9QMdRMzaL4>`__ for
   annotation of tumor suppressor genes and proto-oncogenes

   -  A minimum of n=3 citations were required to include
      literatured-mined tumor suppressor genes and proto-oncogenes from
      CancerMine

0.6.2.1 - May 14th 2018
^^^^^^^^^^^^^^^^^^^^^^^

Fixed
'''''

-  Bug in copy number annotation (broad/focal)

0.6.2 - May 9th 2018
^^^^^^^^^^^^^^^^^^^^

Fixed
'''''

-  Bug in copy number segment display (missing variable initalization)
-  Typo in gnomAD filter statistic (fraction)
-  Bug in mutational signature analysis for grch38 (forgot to pass
   BSgenome object)
-  Missing proper ASCII-encoding in vcf2tsv conversion
-  Removed 'Noncoding mutations' section when no input VCF is present
-  Bug in annotation of copy number event type (focal/broad)
-  Bug in copy number annotation (missing protein-coding transcripts)
-  Updated MSI prediction (variable importance, performance measures)

Added
'''''

-  Genome assembly is appended to every output file
-  Issue warning for copy number segment that goes beyond chromosomal
   lengths of specified assembly (segments will be skipped)
-  Added missing subtypes for 'Skin\_Cancer\_NOS' in the cancer
   phenotype dataset

0.6.1 - May 2nd 2018
^^^^^^^^^^^^^^^^^^^^

Fixed
'''''

-  Bug in tier assignment 'pcgr\_acmg' (case for no variants in
   tier1,2,3)
-  Bug in tier assignment 'pcgr\_acmg' (no tumor type specified,
   evidence items with weak support detected)
-  Bug: duplicated variants in 'Tier 3' resulting from genes encoded
   with dual roles as tumor suppressor genes/oncogenes
-  Bug: duplicated variants in 'Tier 1/Noncoding variants' resulting
   from rare cases of noncoding variants occurring in Tier 1 (synonymous
   variants with biomarker role)

0.6.0 - April 25th 2018
^^^^^^^^^^^^^^^^^^^^^^^

Added
'''''

-  New argument in pcgr.py

   -  *assembly* (grch37/grch38)

-  New option in pcgr.py

   -  *--basic* - run comprehensive VCF annotation only, skip report
      generation and additional analyses

-  New sections in HTML report

   -  *Settings and annotation sources* - now also listing key PCGR
      configuration settings
   -  *Main findings* - Six value boxes indicating the main findings of
      clinical relevance

-  New configuration options

   -  [tier\_model](string) - choice between *pcgr\_acmg* and *pcgr*
   -  [mutational\_burden] - set TMB tertile limits

      -  *tmb\_low\_limit (float)*
      -  *tmb\_intermediate\_limit (float)*

   -  [tumor\_type] - choose between 34 tumor types/classes:

      -  *Adrenal\_Gland\_Cancer\_NOS (logical)*
      -  *Ampullary\_Carcinoma\_NOS (logical)*
      -  *Biliary\_Tract\_Cancer\_NOS (logical)*
      -  *Bladder\_Urinary\_Tract\_Cancer\_NOS (logical)*
      -  *Blood\_Cancer\_NOS (logical)*
      -  *Bone\_Cancer\_NOS (logical)*
      -  *Breast\_Cancer\_NOS (logical)*
      -  *CNS\_Brain\_Cancer\_NOS (logical)*
      -  *Colorectal\_Cancer\_NOS (logical)*
      -  *Cervical\_Cancer\_NOS (logical)*
      -  *Esophageal\_Stomach\_Cancer\_NOS (logical)*
      -  *Head\_And\_Neck\_Cancer\_NOS (logical)*
      -  *Hereditary\_Cancer\_NOS (logical)*
      -  *Kidney\_Cancer\_NOS (logical)*
      -  *Leukemia\_NOS (logical)*
      -  *Liver\_Cancer\_NOS (logical)*
      -  *Lung\_Cancer\_NOS (logical)*
      -  *Lymphoma\_Hodgkin\_NOS (logical)*
      -  *Lymphoma\_Non\_Hodgkin\_NOS (logical)*
      -  *Ovarian\_Fallopian\_Tube\_Cancer\_NOS (logical)*
      -  *Pancreatic\_Cancer\_NOS (logical)*
      -  *Penile\_Cancer\_NOS (logical)*
      -  *Peripheral\_Nervous\_System\_Cancer\_NOS (logical)*
      -  *Peritoneal\_Cancer\_NOS (logical)*
      -  *Pleural\_Cancer\_NOS (logical)*
      -  *Prostate\_Cancer\_NOS (logical)*
      -  *Skin\_Cancer\_NOS (logical)*
      -  *Soft\_Tissue\_Cancer\_NOS (logical)*
      -  *Stomach\_Cancer\_NOS (logical)*
      -  *Testicular\_Cancer\_NOS (logical)*
      -  *Thymic\_Cancer\_NOS (logical)*
      -  *Thyroid\_Cancer\_NOS (logical)*
      -  *Uterine\_Cancer\_NOS (logical)*
      -  *Vulvar\_Vaginal\_Cancer\_NOS (logical)*

   -  [mutational\_signatures]

      -  *mutsignatures\_cutoff (float)* - discard any signature
         contributions with a weight less than the cutoff

   -  [cna]

      -  *transcript\_cna\_overlap (float)* - minimum percent overlap
         between copy number segment and transcripts (average) for tumor
         suppressor gene/proto-oncogene to be reported

   -  [allelic\_support]

      -  If input VCF has correctly formatted depth/allelic fraction as
         INFO tags, users can add thresholds on depth/support that are
         applied prior to report generation

         -  *tumor\_dp\_min (integer)* - minimum sequencing depth for
            variant in tumor sample
         -  *tumor\_af\_min (float)* - minimum allelic fraction for
            variant in tumor sample
         -  *normal\_dp\_min (integer)* - minimum sequencing depth for
            variant in normal sample
         -  *normal\_af\_max (float)* - maximum allelic fraction for
            variant in normal sample

   -  [visual]

      -  *report\_theme (string)* - visual theme of report (Bootstrap)

   -  [other]

      -  *vcf\_validation (logical)* - keep/skip VCF validation by
         `vcf-validator <https://github.com/EBIvariation/vcf-validator>`__

-  New output file - JSON output of HTML report content
-  New INFO tags of PCGR-annotated VCF

   -  *CANCER\_PREDISPOSITION*
   -  *PFAM\_DOMAIN*
   -  *TCGA\_FREQUENCY*
   -  *TCGA\_PANCANCER\_COUNT*
   -  *ICGC\_PCAWG\_OCCURRENCE*
   -  *ICGC\_PCAWG\_AFFECTED\_DONORS*
   -  *CLINVAR\_MEDGEN\_CUI*

-  New column entries in annotated SNV/InDel TSV file:

   -  *CANCER\_PREDISPOSITION*
   -  *ICGC\_PCAWG\_OCCURRENCE*
   -  *TCGA\_FREQUENCY*

-  New column in CNA output

   -  *TRANSCRIPTS* - aberration-overlapping transcripts (Ensembl
      transcript IDs)
   -  *MEAN\_TRANSCRIPT\_CNA\_OVERLAP* - Mean overlap (%) betweeen gene
      transcripts and aberration segment

Removed
'''''''

-  Elements of databundle (now annotated directly through VEP):

   -  dbsnp
   -  gnomad/exac
   -  1000G project

-  INFO tags of PCGR-annotated VCF

   -  *DBSNPBUILDID*
   -  *DBSNP\_VALIDATION*
   -  *DBSNP\_SUBMISSIONS*
   -  *DBSNP\_MAPPINGSTATUS*
   -  *GWAS\_CATALOG\_PMID*
   -  *GWAS\_CATALOG\_TRAIT\_URI*
   -  *DOCM\_DISEASE*

-  Output files

   -  TSV files with mutational signature results and biomarkers (i.e.
      *sample\_id.pcgr.snvs\_indels.biomarkers.tsv* and
      *sample\_id.pcgr.mutational\_signatures.tsv*)

      -  Data can still be retrieved - now from the JSON dump

   -  MAF file

      -  The previous MAF output was generated in a custom fashion, a
         more accurate MAF output based on
         https://github.com/mskcc/vcf2maf will be incorporated in the
         next release

Changed
'''''''

-  HTML report sections

   -  *Tier statistics* and *Variant statistics* are now grouped into
      the section *Tier and variant statistics*
   -  *Tier 5* is now *Noncoding mutations* (i.e. not considered a tier
      per se)
   -  Sliders for allelic fraction in the *Global variant browser* are
      now fixed from 0 to 1 (0.05 intervals)
