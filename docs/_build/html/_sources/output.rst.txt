Output
------

Interactive HTML report
~~~~~~~~~~~~~~~~~~~~~~~

An interactive and tier-structured HTML report that shows the most
relevant findings in the query cancer genome is provided with the
following naming convention:

**sample_id**.__tier_model__.__genome_assembly__.html

-  The **sample_id** is provided as input by the user, and reflects a
   unique identifier of the tumor-normal sample pair to be analyzed.
-  The **tier_model** is provided by the user and will be either
   *‘pcgr’* or *‘pcgr_acmg’*

The report is structured in seven main sections, described in more
detail below:

1. **Report and query settings**

   -  Lists key configurations provided by user

2. **Main results**

   -  Eight value boxes that highlight the main properties of the tumor
      sample:

      1. Tumor purity - as provided by the user
      2. Tumor ploidy - as provided by the user
      3. Mutational signatures - two most prevalent signatures (other
         than aging)
      4. Tier 1 variants (top four)
      5. Tier 2 variants (top four)
      6. Tumor mutational burden (TMB)
      7. Microsatellite instability (MSI) status
      8. Somatic copy number aberrations of clinical significance

3. **Somatic SNVs/InDels**

   -  *Tumor mutational burden (TMB)*

      -  given a coding target region size specified by the user
         (ideally the **callable target size**), an estimate of the
         mutational burden is provided
      -  is presently only computed for tumor-normal input (e.g.
         *vcf_tumor_only = false*)
      -  The estimated mutational burden is assigned a descriptive
         *tertile* based on thresholds defined by the user (these should
         reflect thresholds of clinical significance, and may vary for
         different tumor types)
      -  The estimated TMB is shown in the context of TMB distributions
         from different primary sites in TCGA

   -  *Tier & variant statistics*

      -  indicate total variant numbers across variant types, coding
         types and tiers

   -  *Global distribution - allelic support*

      -  distribution (histogram) of variant allelic support for somatic
         variants (will only be present in the report if specific fields
         in input VCF is defined and specified by the user)

   -  *Global variant browser*

      -  permits exploration of the whole SNV/InDel dataset by filtering
         along several dimensions (call confidence, variant sequencing
         depth/support, variant consequence etc.)

   -  *Tier tables*

      -  Variants are organized into five (tier 1-4 + noncoding)
         interactive datatables) according to clinical utility
      -  Contents of the tier tables are outlined below

4. **Somatic CNAs**

   -  *Segments - amplifications and homozygous deletions*

      -  Based on user-defined/default log-ratio thresholds of
         gains/losses, the whole CNA dataset can be navigated further
         through filters:
      -  cytoband
      -  type of CNA event - *focal* (less than 25% of chromosome arm
         affected) or *broad*
      -  log ratio

   -  *Proto-oncogenes subject to copy number amplifications*

      -  Datatable listing known proto-oncogenes covered by
         user-defined/default amplifications and potential targeted
         therapies

   -  *Tumor suppressor genes subject to homozygous deletions*

      -  Datatable listing known tumor suppressor genes covered by
         user-defined/default losses and potential targeted therapies

   -  *Copy number aberrations as biomarkers for prognosis, diagnosis,
      and drug response*

      -  Interactive data table where the user can navigate aberrations
         acting as biomarkers across therapeutic contexts, tumor
         subtypes, evidence levels etc.

5. **MSI status**

   -  Indicates predicted microsatellite stability from the somatic
      mutation profile and supporting evidence (details of the
      underlying MSI statistical classifier can be found
      `here <http://rpubs.com/sigven/msi_classification_v3>`__)
   -  The MSI classifier was trained on TCGA exome samples.
   -  Will only be present in the report if specified by the user in the
      configuration file (i.e. *msi = true*) and if the input is
      tumor-normal (i.e. *vcf_tumor_only = false*)

6. **Mutational signatures**

   -  Estimation of relative contribution of `30 known mutational
      signatures <http://cancer.sanger.ac.uk/cosmic/signatures>`__ in
      tumor sample (using
      `deconstructSigs <https://github.com/raerose01/deconstructSigs>`__
      as the underlying framework)
   -  Datatable with signatures and proposed underlying etiologies
   -  Will only be present in the report if specified by the user in the
      configuration file (i.e. *mutsignatures = true*) and if the input
      is tumor-normal (i.e. *vcf_tumor_only = false*)
   -  `Trimer (i.e. DNA 3-mer)
      normalization <https://github.com/raerose01/deconstructSigs>`__
      can be configured according to sequencing approach used (WES, WXS
      etc.) using the ‘mutsignatures_normalization’ option, as can the

      -  minimum number of mutations required for analysis (option
         ‘mutsignatures_mutation_limit’)
      -  the maximum number of mutational signatures in the search space
         (option ‘mutsignatures_signature_limit’)
      -  possibility to discard any signature contributions with a
         weight less than a given cutoff (option ‘mutsignatures_cutoff’)

7. **Documentation**

   -  Annotation resources - software, databases and tools with version
      information
   -  References - supporting scientific literature (key report
      elements)

Interactive datatables
^^^^^^^^^^^^^^^^^^^^^^

The interactive datatables contain a number of hyperlinked annotations
similar to those defined for the annotated VCF file, including the
following:

-  SYMBOL - Gene symbol (Entrez/NCBI)
-  PROTEIN_CHANGE - Amino acid change (VEP)
-  CANCER_TYPE - Biomarker (tier 1/2): associated cancer type
-  EVIDENCE_LEVEL - Biomarker (tier 1/2): evidence level (A,B,C,D,E)
-  CLINICAL_SIGNIFICANCE - Biomarker (tier 1/2): drug sensitivity,
   poor/better outcome etc
-  EVIDENCE_TYPE - Biomarker (tier 1/2):
   predictive/diagnostic/therapeutic
-  DISEASE_ONTOLOGY_ID - Biomarker (tier 1/2): associated cancer type
   (Disease Ontology)
-  EVIDENCE_DIRECTION - Biomarker (tier 1/2): supports/does not support
-  DESCRIPTION - Biomarker (tier 1/2): description
-  VARIANT_ORIGIN - Biomarker (tier 1/2): variant origin
   (germline/somatic)
-  BIOMARKER_MAPPING - Biomarker (tier 1/2): accuracy of genomic mapping
   (exact,codon,exon)
-  CITATION - Biomarker (tier 1/2): supporting literature
-  THERAPEUTIC_CONTEXT - Biomarker (tier 1/2): associated drugs
-  RATING - Biomarker (tier 1/2): trust rating from 1 to 5 (CIVIC)
-  GENE_NAME - gene name description (Entrez/NCBI)
-  PROTEIN_DOMAIN - PFAM protein domain
-  PROTEIN_FEATURE - UniProt feature overlapping variant site
-  CDS_CHANGE - Coding sequence change
-  MUTATION_HOTSPOT - Known cancer mutation hotspot
-  MUTATION_HOTSPOT_CANCERTYPE - Hotspot-associated cancer types
-  TCGA_FREQUENCY - Frequency of variant in TCGA cohorts
-  ICGC_PCAWG_OCCURRENCE - Frequency of variant in ICGC-PCAWG cohorts
-  DOCM_LITERATURE - Literature links - DoCM
-  DOCM_DISEASE - Associated diseases - DoCM
-  OPENTARGETS_RANK - Strength of gene-phenotype associatino according
   to the Open Targets Platform
-  OPENTARGETS_ASSOCIATIONS - Phenotype associations with the gene
   retrieved from the Open Targets Platform
-  INTOGEN_DRIVER_MUT - predicted driver mutation - IntOGen
-  CONSEQUENCE - VEP consequence (primary transcript)
-  HGVSc - from VEP
-  HGVSp - from VEP
-  ONCOGENE - Predicted as proto-oncogene from literature mining
-  TUMOR_SUPPRESSOR - Predicted as tumor suppressor gene from literature
   mining
-  ONCOSCORE - Literature-derived score for oncogenic potential (gene
   level)
-  PREDICTED_EFFECT - Effect predictions from dbNSFP
-  VEP_ALL_CSQ - All VEP transcript block consequences
-  DBSNP - dbSNP rsID
-  COSMIC - Cosmic mutation IDs
-  CLINVAR - ClinVar variant origin and associated phenotypes
-  CANCER_ASSOCIATIONS - Gene-associated cancer types from DisGenet
-  TARGETED_DRUGS - Targeted drugs from DGIdb-ChEMBL
-  KEGG_PATHWAY - Gene-associated pathways from KEGG
-  CALL_CONFIDENCE - Variant confidence (as set by user in input VCF)
-  DP_TUMOR - Variant sequencing depth in tumor (as set by user in input
   VCF)
-  AF_TUMOR - Variant allelic fraction in tumor (as set by user in input
   VCF)
-  DP_CONTROL - Variant sequencing depth in control sample (as set by
   user in input VCF)
-  AF_CONTROL - Variant allelic fraction in control sample (as set by
   user in input VCF)
-  GENOMIC_CHANGE - Variant ID
-  GENOME_VERSION - Genome assembly

Example reports:

-  `View an example report for a breast tumor sample
   (TCGA) <http://folk.uio.no/sigven/tumor_sample.BRCA.pcgr_acmg.grch37.0.8.0.html>`__
-  `View an example report for a colon adenocarcinoma sample
   (TCGA) <http://folk.uio.no/sigven/tumor_sample.COAD.pcgr_acmg.grch37.0.8.0.html>`__

The HTML reports have been tested using the following browsers:

-  Safari (Version 12.1 (14607.1.40.1.4))
-  Mozilla Firefox (52.0.2)
-  Google Chrome (Version 74.0.3729.131 )

JSON (beta)
~~~~~~~~~~~

A JSON file that stores the HTML report content is provided. This file
will easen the process of extracting particular parts of the report for
further analysis. The JSON contains two main objects, *metadata* and
*content*, where the former contains information about the settings,
data versions, and the latter contains the various sections of the
report. Examples (using R) on how to extract information from the JSON
file will soon be posted here.

Output files - somatic SNVs/InDels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Variant call format - VCF
^^^^^^^^^^^^^^^^^^^^^^^^^

A VCF file containing annotated, somatic calls (single nucleotide
variants and insertion/deletions) is generated with the following naming
convention:

**sample_id**.__tier_model__.__genome_assembly__.vcf.gz

Here, the **sample_id** is provided as input by the user, and reflects a
unique identifier of the tumor-normal sample pair to be analyzed.
Following common standards, the annotated VCF file is compressed with
`bgzip <http://www.htslib.org/doc/tabix.html>`__ and indexed with
`tabix <http://www.htslib.org/doc/tabix.html>`__. Below follows a
description of all annotations/tags present in the VCF INFO column after
processing with the PCGR annotation pipeline:

*VEP consequence annotations*
'''''''''''''''''''''''''''''

-  CSQ - Complete consequence annotations from VEP. Format:
   Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON\|
   INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids\|
   Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS\|
   SYMBOL_SOURCE|HGNC_ID|CANONICAL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC\|
   RefSeq|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF\|
   gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF\|
   gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO\|
   MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE
-  Consequence - Impact modifier for the consequence type (picked by
   VEP’s –flag_pick_allele option)
-  Gene - Ensembl stable ID of affected gene (picked by VEP’s
   –flag_pick_allele option)
-  Feature_type - Type of feature. Currently one of Transcript,
   RegulatoryFeature, MotifFeature (picked by VEP’s –flag_pick_allele
   option)
-  Feature - Ensembl stable ID of feature (picked by VEP’s
   –flag_pick_allele option)
-  cDNA_position - Relative position of base pair in cDNA sequence
   (picked by VEP’s –flag_pick_allele option)
-  CDS_position - Relative position of base pair in coding sequence
   (picked by VEP’s –flag_pick_allele option)
-  CDS_CHANGE - Coding, transcript-specific sequence annotation (picked
   by VEP’s –flag_pick_allele option)
-  AMINO_ACID_START - Protein position indicating absolute start of
   amino acid altered (fetched from Protein_position)
-  AMINO_ACID_END - Protein position indicating absolute end of amino
   acid altered (fetched from Protein_position)
-  Protein_position - Relative position of amino acid in protein (picked
   by VEP’s –flag_pick_allele option)
-  Amino_acids - Only given if the variant affects the protein-coding
   sequence (picked by VEP’s –flag_pick_allele option)
-  Codons - The alternative codons with the variant base in upper case
   (picked by VEP’s –flag_pick_allele option)
-  IMPACT - Impact modifier for the consequence type (picked by VEP’s
   –flag_pick_allele option)
-  VARIANT_CLASS - Sequence Ontology variant class (picked by VEP’s
   –flag_pick_allele option)
-  SYMBOL - Gene symbol (picked by VEP’s –flag_pick_allele option)
-  SYMBOL_ENTREZ - Official gene symbol as provided by NCBI’s Entrez
   gene
-  SYMBOL_SOURCE - The source of the gene symbol (picked by VEP’s
   –flag_pick_allele option)
-  STRAND - The DNA strand (1 or -1) on which the transcript/feature
   lies (picked by VEP’s –flag_pick_allele option)
-  ENSP - The Ensembl protein identifier of the affected transcript
   (picked by VEP’s –flag_pick_allele option)
-  FLAGS - Transcript quality flags: cds_start_NF: CDS 5’, incomplete
   cds_end_NF: CDS 3’ incomplete (picked by VEP’s –flag_pick_allele
   option)
-  SWISSPROT - Best match UniProtKB/Swiss-Prot accession of protein
   product (picked by VEP’s –flag_pick_allele option)
-  TREMBL - Best match UniProtKB/TrEMBL accession of protein product
   (picked by VEP’s –flag_pick_allele option)
-  UNIPARC - Best match UniParc accession of protein product (picked by
   VEP’s –flag_pick_allele option)
-  HGVSc - The HGVS coding sequence name (picked by VEP’s
   –flag_pick_allele option)
-  HGVSp - The HGVS protein sequence name (picked by VEP’s
   –flag_pick_allele option)
-  HGVSp_short - The HGVS protein sequence name, short version (picked
   by VEP’s –flag_pick_allele option)
-  HGVS_OFFSET - Indicates by how many bases the HGVS notations for this
   variant have been shifted (picked by VEP’s –flag_pick_allele option)
-  MOTIF_NAME - The source and identifier of a transcription factor
   binding profile aligned at this position (picked by VEP’s
   –flag_pick_allele option)
-  MOTIF_POS - The relative position of the variation in the aligned
   TFBP (picked by VEP’s –flag_pick_allele option)
-  HIGH_INF_POS - A flag indicating if the variant falls in a high
   information position of a transcription factor binding profile (TFBP)
   (picked by VEP’s –flag_pick_allele option)
-  MOTIF_SCORE_CHANGE - The difference in motif score of the reference
   and variant sequences for the TFBP (picked by VEP’s –flag_pick_allele
   option)
-  CELL_TYPE - List of cell types and classifications for regulatory
   feature (picked by VEP’s –flag_pick_allele option)
-  CANONICAL - A flag indicating if the transcript is denoted as the
   canonical transcript for this gene (picked by VEP’s –flag_pick_allele
   option)
-  CCDS - The CCDS identifier for this transcript, where applicable
   (picked by VEP’s –flag_pick_allele option)
-  INTRON - The intron number (out of total number) (picked by VEP’s
   –flag_pick_allele option)
-  EXON - The exon number (out of total number) (picked by VEP’s
   –flag_pick_allele option)
-  LAST_EXON - Logical indicator for last exon of transcript (picked by
   VEP’s –flag_pick_allele option)
-  LAST_INTRON - Logical indicator for last intron of transcript (picked
   by VEP’s –flag_pick_allele option)
-  DISTANCE - Shortest distance from variant to transcript (picked by
   VEP’s –flag_pick_allele option)
-  BIOTYPE - Biotype of transcript or regulatory feature (picked by
   VEP’s –flag_pick_allele option)
-  TSL - Transcript support level (picked by VEP’s –flag_pick_allele
   option)>
-  PUBMED - PubMed ID(s) of publications that cite existing variant -
   VEP
-  PHENO - Indicates if existing variant is associated with a phenotype,
   disease or trait - VEP
-  GENE_PHENO - Indicates if overlapped gene is associated with a
   phenotype, disease or trait - VEP
-  ALLELE_NUM - Allele number from input; 0 is reference, 1 is first
   alternate etc - VEP
-  REFSEQ_MATCH - The RefSeq transcript match status; contains a number
   of flags indicating whether this RefSeq transcript matches the
   underlying reference sequence and/or an Ensembl transcript (picked by
   VEP’s –flag_pick_allele option)
-  PICK - Indicates if this block of consequence data was picked by
   VEP’s –flag_pick_allele option
-  VEP_ALL_CONSEQUENCE - All transcript consequences
   (Consequence:SYMBOL:Feature_type:Feature:BIOTYPE) - VEP
-  EXONIC_STATUS - Indicates if variant consequence type is ‘exonic’ or
   ‘nonexonic’. We here define ‘exonic’ as any variant with either of
   the following consequence:

   -  stop_gained / stop_lost
   -  start_lost
   -  frameshift_variant
   -  missense_variant
   -  splice_donor_variant
   -  splice_acceptor_variant
   -  inframe_insertion / inframe_deletion
   -  synonymous_variant
   -  protein_altering

-  CODING_STATUS - Indicates if primary variant consequence type is
   ‘coding’ or ‘noncoding’ (wrt. protein-alteration). ‘coding’ variants
   are here defined as those with an ‘exonic’ status, with the exception
   of synonymous variants

*Gene information*
''''''''''''''''''

-  ENTREZ_ID - `Entrez <http://www.ncbi.nlm.nih.gov/gene>`__ gene
   identifier
-  APPRIS - Principal isoform flags according to the `APPRIS principal
   isoform database <http://appris.bioinfo.cnio.es/#/downloads>`__
-  UNIPROT_ID - `UniProt <http://www.uniprot.org>`__ identifier
-  UNIPROT_ACC - `UniProt <http://www.uniprot.org>`__ accession(s)
-  ENSEMBL_GENE_ID - Ensembl gene identifier for VEP’s picked transcript
   (*ENSGXXXXXXX*)
-  ENSEMBL_TRANSCRIPT_ID - Ensembl transcript identifier for VEP’s
   picked transcript (*ENSTXXXXXX*)
-  REFSEQ_MRNA - Corresponding RefSeq transcript(s) identifier for VEP’s
   picked transcript (*NM_XXXXX*)
-  CORUM_ID - Associated protein complexes (identifiers) from
   `CORUM <http://mips.helmholtz-muenchen.de/corum/>`__
-  DISGENET_CUI - Tumor types associated with gene, as found in
   DisGeNET. Tumor types are listed as unique
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept IDs
   (*CUIs*)
-  TUMOR_SUPPRESSOR - Gene is predicted as tumor suppressor candidate
   according to
   (`CancerMine <https://zenodo.org/record/2662509#.XNM4VtMzaL5>`__)
-  ONCOGENE - Gene is predicted as an oncogene according to
   (`CancerMine <https://zenodo.org/record/2662509#.XNM4VtMzaL5>`__)
-  ONCOSCORE - Literature-derived score for cancer gene relevance
   `Bioconductor/OncoScore <http://bioconductor.org/packages/release/bioc/html/OncoScore.html>`__,
   range from 0 (low oncogenic potential) to 1 (high oncogenic
   potential)
-  INTOGEN_DRIVER - Gene is predicted as a cancer driver in the `IntoGen
   Cancer Drivers Database <https://www.intogen.org/downloads>`__
-  TCGA_DRIVER - Gene is predicted as a cancer driver in the `TCGA
   pan-cancer analysis of cancer driver genes and
   mutations <https://www.ncbi.nlm.nih.gov/pubmed/29625053>`__

*Variant effect and protein-coding information*
'''''''''''''''''''''''''''''''''''''''''''''''

-  MUTATION_HOTSPOT - mutation hotspot codon in
   `cancerhotspots.org <http://cancerhotspots.org/>`__. Format:
   gene_symbol \| codon \| q-value
-  MUTATION_HOTSPOT_TRANSCRIPT - hotspot-associated transcripts (Ensembl
   transcript ID)
-  MUTATION_HOTSPOT_CANCERTYPE - hotspot-associated cancer types (from
   cancerhotspots.org)
-  UNIPROT_FEATURE - Overlapping protein annotations from `UniProt
   KB <http://www.uniprot.org>`__
-  PFAM_DOMAIN - Pfam domain identifier (from VEP)
-  INTOGEN_DRIVER_MUT - Indicates if existing variant is predicted as
   driver mutation from IntoGen Catalog of Driver Mutations
-  PUTATIVE_DRIVER_MUTATION - Variant is predicted as driver mutation in
   the `TCGA pan-cancer analysis of cancer driver genes and
   mutations <https://www.ncbi.nlm.nih.gov/pubmed/29625053>`__
-  EFFECT_PREDICTIONS - All predictions of effect of variant on protein
   function and pre-mRNA splicing from `database of non-synonymous
   functional predictions -
   dbNSFP <https://sites.google.com/site/jpopgen/dbNSFP>`__. Predicted
   effects are provided by different sources/algorithms (separated by
   ‘&’):

   1.  `SIFT <https://sift.bii.a-star.edu.sg/>`__
   2.  `SIFT4G <https://sift.bii.a-star.edu.sg/sift4g/>`__
   3.  `LRT <http://www.genetics.wustl.edu/jflab/lrt_query.html>`__
       (2009)
   4.  `MutationTaster <http://www.mutationtaster.org/>`__ (data release
       Nov 2015)
   5.  `MutationAssessor <http://mutationassessor.org/>`__ (release 3)
   6.  `FATHMM <http://fathmm.biocompute.org.uk>`__ (v2.3)
   7.  `PROVEAN <http://provean.jcvi.org/index.php>`__ (v1.1 Jan 2015)
   8.  `FATHMM_MKL <http://fathmm.biocompute.org.uk/fathmmMKL.htm>`__
   9.  `PRIMATEAI <https://www.nature.com/articles/s41588-018-0167-z>`__
   10. `DEOGEN2 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5570203/>`__
   11. `DBNSFP_CONSENSUS_SVM <https://www.ncbi.nlm.nih.gov/pubmed/25552646>`__
       (Ensembl/consensus prediction, based on support vector machines)
   12. `DBNSFP_CONSENSUS_LR <https://www.ncbi.nlm.nih.gov/pubmed/25552646>`__
       (Ensembl/consensus prediction, logistic regression based)
   13. `SPLICE_SITE_EFFECT_ADA <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       adaptive boosting)
   14. `SPLICE_SITE_EFFECT_RF <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       random forest)
   15. `M-CAP <http://bejerano.stanford.edu/MCAP>`__
   16. `MutPred <http://mutpred.mutdb.org>`__
   17. `GERP <http://mendel.stanford.edu/SidowLab/downloads/gerp/>`__

-  SIFT_DBNSFP - predicted effect from SIFT (dbNSFP)
-  SIFT4G_DBNSFP - predicted effect from SIFT4G (dbNSFP)
-  PROVEAN_DBNSFP - predicted effect from PROVEAN (dbNSFP)
-  MUTATIONTASTER_DBNSFP - predicted effect from MUTATIONTASTER (dbNSFP)
-  MUTATIONASSESSOR_DBNSFP - predicted effect from MUTATIONASSESSOR
   (dbNSFP)
-  M_CAP_DBNSFP - predicted effect from M-CAP (dbNSFP)
-  MUTPRED_DBNSFP - score from MUTPRED (dbNSFP)
-  FATHMM_DBNSFP - predicted effect from FATHMM (dbNSFP)
-  PRIMATEAI_DBNSFP - predicted effect from PRIMATEAI (dbNSFP)
-  DEOGEN2_DBNSFP - predicted effect from DEOGEN2 (dbNSFP)
-  FATHMM_MKL_DBNSFP - predicted effect from FATHMM-mkl (dbNSFP)
-  META_LR_DBNSFP - predicted effect from ensemble prediction (logistic
   regression - dbNSFP)
-  SPLICE_SITE_RF_DBNSFP - predicted effect of splice site disruption,
   using random forest (dbscSNV)
-  SPLICE_SITE_ADA_DBNSFP - predicted effect of splice site disruption,
   using boosting (dbscSNV)

*Variant frequencies/annotations in germline/somatic databases*
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

-  AFR_AF_GNOMAD - African/American germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  AMR_AF_GNOMAD - American germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  GLOBAL_AF_GNOMAD - Adjusted global germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  SAS_AF_GNOMAD - South Asian germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  EAS_AF_GNOMAD - East Asian germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  FIN_AF_GNOMAD - Finnish germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  NFE_AF_GNOMAD - Non-Finnish European germline allele frequency
   (`Genome Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  OTH_AF_GNOMAD - Other germline allele frequency (`Genome Aggregation
   Database release 2 <http://gnomad.broadinstitute.org/>`__)
-  ASJ_AF_GNOMAD - Ashkenazi Jewish allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  AFR_AF_1KG - `1000G Project - phase 3 <http://www.1000genomes.org>`__
   germline allele frequency for samples from AFR (African)
-  AMR_AF_1KG - `1000G Project - phase 3 <http://www.1000genomes.org>`__
   germline allele frequency for samples from AMR (Ad Mixed American)
-  EAS_AF_1KG - `1000G Project - phase 3 <http://www.1000genomes.org>`__
   germline allele frequency for samples from EAS (East Asian)
-  EUR_AF_1KG - `1000G Project - phase 3 <http://www.1000genomes.org>`__
   germline allele frequency for samples from EUR (European)
-  SAS_AF_1KG - `1000G Project - phase 3 <http://www.1000genomes.org>`__
   germline allele frequency for samples from SAS (South Asian)
-  GLOBAL_AF_1KG - `1000G Project - phase
   3 <http://www.1000genomes.org>`__ germline allele frequency for all
   1000G project samples (global)
-  DBSNPRSID - `dbSNP <http://www.ncbi.nlm.nih.gov/SNP/>`__ reference
   ID, as provided by VEP
-  COSMIC_MUTATION_ID - Mutation identifier in `Catalog of somatic
   mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__
   database, as provided by VEP
-  TCGA_PANCANCER_COUNT - Raw variant count across all TCGA tumor types
-  TCGA_FREQUENCY - Frequency of variant across TCGA tumor types.
   Format: tumortype\| percent affected|affected cases|total cases
-  ICGC_PCAWG_OCCURRENCE - Mutation occurrence in
   `ICGC-PCAWG <http://docs.icgc.org/pcawg/>`__. By project:
   project_code|affected_donors|tested_donors|frequency)
-  ICGC_PCAWG_AFFECTED_DONORS - Number of donors with the current
   mutation in `ICGC-PCAWG <http://docs.icgc.org/pcawg/>`__

*Clinical associations*
'''''''''''''''''''''''

-  CLINVAR_MSID - `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
   Measure Set/Variant ID
-  CLINVAR_ALLELE_ID - `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
   allele ID
-  CLINVAR_PMID - Associated Pubmed IDs for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - germline
   state-of-origin
-  CLINVAR_HGVSP - Protein variant expression using HGVS nomenclature
-  CLINVAR_PMID_SOMATIC - Associated Pubmed IDs for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - somatic
   state-of-origin
-  CLINVAR_CLNSIG - Clinical significance for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - germline
   state-of-origin
-  CLINVAR_CLNSIG_SOMATIC - Clinical significance for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - somatic
   state-of-origin
-  CLINVAR_MEDGEN_CUI - Associated
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept identifiers
   (*CUIs*) - germline state-of-origin
-  CLINVAR_MEDGEN_CUI_SOMATIC - Associated
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept identifiers
   (*CUIs*) - somatic state-of-origin
-  CLINVAR_VARIANT_ORIGIN - Origin of variant (somatic, germline, de
   novo etc.) for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
-  CLINVAR_REVIEW_STATUS_STARS - Rating of the
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ variant (0-4 stars)
   with respect to level of review
-  DOCM_PMID - Associated Pubmed IDs for variant in `Database of Curated
   Mutations <http://docm.genome.wustl.edu>`__
-  OPENTARGETS_DISEASE_ASSOCS - Associations between protein targets and
   disease based on multiple lines of evidence (mutations,affected
   pathways,GWAS, literature etc). Format:
   CUI:EFO_ID:IS_DIRECT:OVERALL_SCORE
-  OPENTARGETS_TRACTABILITY_COMPOUND - Confidence for the existence of a
   modulator (small molecule) that interacts with the target to elicit a
   desired biological effect
-  OPENTARGETS_TRACTABILITY_ANTIBODY - Confidence for the existence of a
   modulator (antibody) that interacts with the target to elicit a
   desired biological effect

*Other*
'''''''

-  CHEMBL_COMPOUND_ID - antineoplastic drugs targeting the encoded
   protein (from `Drug-Gene Interaction
   Database <http://dgidb.genome.wustl.edu/>`__, drugs are listed as
   `ChEMBL <https://www.ebi.ac.uk/chembl/>`__ compound identifiers)
-  CIVIC_ID, CIVIC_ID_2 - Variant identifiers in the `CIViC
   database <http://civic.genome.wustl.edu>`__, CIVIC_ID refers to
   markers mapped at variant level, CIVIC_ID_2 refers to region markers
   (codon, exon etc.)
-  CBMDB_ID - Variant identifier in the `Cancer Biomarkers
   database <https://www.cancergenomeinterpreter.org/biomarkers>`__

Tab-separated values (TSV)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotated List of all SNVs/InDels
'''''''''''''''''''''''''''''''''

We provide a tab-separated values file with most important annotations
for SNVs/InDels. The file has the following naming convention:

**sample_id**.__tier_model__.__genome_assembly__.snvs_indels.tiers.tsv

The SNVs/InDels are organized into different **tiers** (as defined above
for the HTML report)

The following variables are included in the tiered TSV file:

::

   1. CHROM - Chromosome
   2. POS - Position (VCF-based)
   3. REF - Reference allele
   4. ALT - Alternate allele
   5. GENOMIC_CHANGE - Identifier for variant at the genome (VCF) level, e.g. 1:g.152382569A>G
         Format: (<chrom>:g.<position><ref_allele>><alt_allele>)
   6. GENOME_VERSION - Assembly version, e.g. GRCh37
   7. VCF_SAMPLE_ID - Sample identifier
   8. VARIANT_CLASS - Variant type, e.g. SNV/insertion/deletion
   9. SYMBOL - Gene symbol
   10. GENE_NAME - Gene description
   11. CCDS - CCDS identifier
   12. CANONICAL - indication of canonical transcript
   13. ENTREZ_ID - Entrez gene identifier
   14. UNIPROT_ID - UniProt protein identifier
   15. ENSEMBL_TRANSCRIPT_ID - Ensembl transcript identifier
   16. ENSEMBL_GENE_ID - Ensembl gene identifier
   17. REFSEQ_MRNA - RefSeq mRNA identifier
   18. ONCOSCORE - Literature-derived score for cancer gene relevance
   19. ONCOGENE - Gene is predicted as an oncogene according to literature mining (CancerMine)
   20. TUMOR_SUPPRESSOR - Gene is predicted as tumor suppressor according to literature mining (CancerMine)
   21. DISGENET_CUI - Associated tumor types from DisGeNET (MedGen concept IDs)
   22. DISGENET_TERMS - Associated tumor types from DisGeNET (MedGen concept terms)
   23. CONSEQUENCE - Variant consequence (as defined above for VCF output:
       Consequence)
   24. PROTEIN_CHANGE - Protein change (HGVSp without reference accession)
   25. PROTEIN_DOMAIN - Protein domain
   26. CODING_STATUS - Coding variant status wrt. protein alteration ('coding' or 'noncoding')
   27. EXONIC_STATUS - Exonic variant status ('exonic' or 'nonexonic')
   28. CDS_CHANGE - composite VEP-based variable for coding change, format:
       Consequence:Feature:cDNA_position:EXON:HGVSp_short
   29. HGVSp
   30. HGVSc
   31. EFFECT_PREDICTIONS - as defined above for VCF
   32. MUTATION_HOTSPOT - mutation hotspot codon in
       cancerhotspots.org. Format: gene_symbol | codon | q-value
   33. MUTATION_HOTSPOT_TRANSCRIPT - hotspot-associated transcripts (Ensembl transcript ID)
   34. MUTATION_HOTSPOT_CANCERTYPE - hotspot-associated cancer types (from cancerhotspots.org)
   35. PUTATIVE_DRIVER_MUTATION - Indicates if variant is predicted as
       driver mutation from TCGA's PanCancer study of cancer driver mutation
   36. VEP_ALL_CSQ - all VEP transcript block consequences
   37. DBSNPRSID - dbSNP reference cluster ID
   38. COSMIC_MUTATION_ID - COSMIC mutation ID
   39. TCGA_PANCANCER_COUNT - Raw variant count across all TCGA tumor types
   40. TCGA_FREQUENCY - Frequency of variant across TCGA tumor types. Format: tumortype|
   percent affected|affected cases|total cases
   41. ICGC_PCAWG_OCCURRENCE - Mutation occurrence in ICGC-PCAWG by project:
   project_code|affected_donors|tested_donors|frequency
   42. CHEMBL_COMPOUND_ID - Compounds (as ChEMBL IDs) that target the encoded protein (from DGIdb)
   43. CHEMBL_COMPOUND_TERMS - Compounds (as drug names) that target the encoded protein (from DGIdb)
   44. SIMPLEREPEATS_HIT - Variant overlaps UCSC _simpleRepeat_ sequence repeat track
   45. WINMASKER_HIT - Variant overlaps UCSC _windowmaskerSdust_ sequence repeat track
   46. OPENTARGETS_RANK - OpenTargets association score (between 0 and 1) for gene (maximum across cancer phenotypes)
   47. CLINVAR - ClinVar association: variant origin and associated traits
   48. CLINVAR_CLNSIG - clinical significance of ClinVar variant
   49. GLOBAL_AF_GNOMAD - global germline allele frequency in gnomAD
   50. GLOBAL_AF_1KG - 1000G Project - phase 3, germline allele frequency
   51. CALL_CONFIDENCE - confidence indicator for somatic variant
   52. DP_TUMOR - sequencing depth at variant site (tumor sample)
   53. AF_TUMOR - allelic fraction of alternate allele (tumor sample)
   54. DP_CONTROL - sequencing depth at variant site (control sample)
   55. AF_CONTROL - allelic fraction of alternate allele (control sample)
   56. TIER
   57. TIER_DESCRIPTION

**NOTE**: The user has the possibility to append the TSV file with data
from other tags in the input VCF of interest (i.e. using the
*custom_tags* option in the TOML configuration file)

Output files - somatic copy number aberrations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _tab-separated-values-tsv-1:

1. Tab-separated values (TSV)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copy number segments are intersected with the genomic coordinates of all
transcripts from `GENCODE’s basic gene
annotation <https://www.gencodegenes.org/releases/current.html>`__. In
addition, PCGR attaches cancer-relevant annotations for the affected
transcripts. The naming convention of the compressed TSV file is as
follows:

**sample_id**.__tier_model__.__genome_assembly__.cna_segments.tsv.gz

The format of the compressed TSV file is the following:

::

   1. chrom - chromosome
   2. segment_start - start of copy number segment
   3. segment_end - end of copy number segment
   4. segment_length_Mb - length of segment in Mb
   5. event_type - focal or broad (covering more than 25% of chromosome arm)
   6. cytoband
   7. LogR - Copy log-ratio
   8. ensembl_gene_id
   9. symbol - gene symbol
   10. ensembl_transcript_id
   11. transcript_start
   12. transcript_end
   13. transcript_overlap_percent - percent of transcript length covered by CN segment
   14. name - gene name description
   15. biotype - type of gene
   16. disgenet_cui - tumor types associated with gene (from DisGeNET, tumor types
      are listed as MedGen concept IDs (CUI)
   17. tsgene - tumor suppressor gene status (CancerMine literature database)
   18. p_oncogene - oncogene status (CancerMine literature database)
   19. intogen_drivers - predicted driver gene status (IntoGen Cancer Drivers Database)
   20. chembl_compound_id - antineoplastic drugs targeting the encoded protein
      (from DGIdb, drugs are listed as ChEMBL compound identifiers)
   21. gencode_gene_biotype
   22. gencode_tag
   23. gencode_release
