Output
------

Interactive HTML report (flexdashboard)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An interactive and tier-structured HTML report (dashboard-like) that
shows the most relevant findings in the query cancer genome has the
following naming convention:

**sample_id**.pcgr_acmg.__genome_assembly__.flexdb.html

-  The **sample_id** is provided as input by the user, and reflects a
   unique identifier of the tumor-normal sample pair to be analyzed.

-  DOCUMENTATION COMING SOON

Interactive HTML report (rmarkdown)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An interactive and tier-structured HTML report that shows the most
relevant findings in the query cancer genome has the following naming
convention:

**sample_id**.pcgr_acmg.__genome_assembly__.html

-  The **sample_id** is provided as input by the user, and reflects a
   unique identifier of the tumor-normal sample pair to be analyzed.

The report is structured in seven main sections, described in more
detail below:

1. **Report and query settings**

   -  Lists key configurations provided by user

2. **Main results**

   -  Nine value boxes that highlight the main properties of the tumor
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
      9. Kataegis events

3. **Somatic SNVs/InDels**

   -  *Tumor mutational burden (TMB)*

      -  given a coding target region size specified by the user
         (ideally the **callable target size**), an estimate of the
         mutational burden is provided
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

   -  *Other drug targets subject to copy number aplifications*
   -  *Copy number aberrations as biomarkers for prognosis, diagnosis,
      and drug response (both of strong significance and potential
      significance)*

      -  Interactive data table where the user can navigate aberrations
         acting as biomarkers across therapeutic contexts, tumor
         subtypes, evidence levels etc

5. **MSI status**

   -  Indicates predicted microsatellite stability from the somatic
      mutation profile and supporting evidence (details of the
      underlying MSI statistical classifier can be found
      `here <http://rpubs.com/sigven/msi_classification_v3>`__)
   -  The MSI classifier was trained on TCGA exome samples.

6. **Mutational signatures**

   -  Estimation of relative contribution of `67 known mutational
      signatures <http://cancer.sanger.ac.uk/cosmic/signatures>`__ in
      tumor sample (using
      `MutationalPatterns <https://github.com/raerose01/deconstructSigs>`__
      as the underlying framework)
   -  Datatable with signatures and proposed underlying etiologies

7. **Kataegis events**

   -  Detection of potential kataegis events using the algorithm
      outlined in the `Kataegis
      Portal <https://github.com/MeichunCai/KataegisPortal>`__

8. **Clinial trials**

   -  Interactive table with annotated, molecularly targeted clinical
      trials from `clinicaltrials.gov <https//clinicaltrials.gov>`__ in
      the relevant tumor type

9. **Documentation**

   -  Annotation resources - software, databases and tools with version
      information
   -  Report content
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
-  OPENTARGETS_RANK - Strength (max across all phenotypes) of phenotype
   association according to the Open Targets Platform
-  OPENTARGETS_ASSOCIATIONS - Phenotype associations with the gene
   retrieved from the Open Targets Platform
-  INTOGEN_DRIVER_MUT - predicted driver mutation - IntOGen
-  CONSEQUENCE - VEP consequence (primary transcript)
-  HGVSc - from VEP
-  HGVSp - from VEP
-  NCBI_REFSEQ - Transcript accession ID(s) (NCBI RefSeq)
-  ONCOGENE - Predicted as proto-oncogene from CancerMine/NCG
-  CANCERGENE_SUPPORT - Links to underlying publications (CancerMine)
   that support oncogenic/tumor suppressive role of gene
-  TUMOR_SUPPRESSOR - Predicted as tumor suppressor gene from
   CancerMine/NCG
-  PREDICTED_EFFECT - Effect predictions from dbNSFP
-  VEP_ALL_CSQ - All VEP transcript block consequences
-  DBSNP - dbSNP rsID
-  COSMIC - Cosmic mutation IDs
-  CLINVAR - ClinVar variant origin and associated phenotypes
-  CANCER_ASSOCIATIONS - Gene-associated cancer types from DisGenet
-  TARGETED_DRUGS - Targeted drugs from Open Targets Platform /ChEMBL
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

-  `Cervical cancer sample
   (tumor-only) <http://insilico.hpc.uio.no/pcgr/example_reports/0.9.1/TCGA-EA-A410-01A_TO.pcgr_acmg.grch37.flexdb.html>`__
-  `Lung cancer sample
   (tumor-control) <http://insilico.hpc.uio.no/pcgr/example_reports/0.9.1/TCGA-05-4427-01A.pcgr_acmg.grch37.flexdb.html>`__
-  `Colorectal cancer sample
   (tumor-control) <http://insilico.hpc.uio.no/pcgr/example_reports/0.9.1/TCGA-AD-5900-01A.pcgr_acmg.grch37.flexdb.html>`__
-  `Brain cancer sample
   (tumor-control) <http://insilico.hpc.uio.no/pcgr/example_reports/0.9.1/TCGA-QH-A6CU-01A.pcgr_acmg.grch37.flexdb.html>`__

The HTML reports have been tested using the following browsers:

-  Safari (Version 14.0 (15610.1.28.1.9, 15610))
-  Mozilla Firefox (83.0)
-  Google Chrome (Version 87.0.4280.67 (Official Build) (x86_64))

JSON
~~~~

A compressed JSON file that stores all the essential content of the
report is provided.

This file will easen the process of extracting particular parts of the
report for further analysis or integration with other workflows. The
JSON contains two main objects, *metadata* and *content*, where the
former contains information about the settings, data versions, and the
latter contains the various sections of the report.

Output files - somatic SNVs/InDels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Variant call format - VCF
^^^^^^^^^^^^^^^^^^^^^^^^^

A VCF file containing annotated, somatic calls (single nucleotide
variants and insertion/deletions) is generated with the following naming
convention:

**sample_id**.pcgr_acmg.__genome_assembly__.vcf.gz

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
   MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|NearestExonJB
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
-  NearestExonJB - VEP plugin that finds nearest exon junction for a
   coding sequence variant. Format: Ensembl exon identifier+distanceto
   exon boundary+boundary type(start/end)+exon length
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
-  INTRON_POSITION - Relative position of intron variant to nearest
   exon/intron junction (NearestExonJB VEP plugin)
-  EXON_POSITION - Relative position of exon variant to nearest
   intron/exon junction (NearestExonJB VEP plugin)
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
-  TUMOR_SUPPRESSOR - Indicates whether gene is predicted as a tumor
   suppressor gene, from Network of Cancer Genes (NCG) & the CancerMine
   text-mining resource
-  TUMOR_SUPPRESSOR_EVIDENCE - Underlying evidence for gene being a
   tumor suppressor. Format:
   NCG:<TRUE|FALSE>&CancerMine:<LC|MC|HC>:num_citations
-  ONCOGENE - Indicates whether gene is predicted as an oncogene, from
   Network of Cancer Genes (NCG) & the CancerMine text-mining resource
-  ONCOGENE_EVIDENCE - Underlying evidence for gene being an oncogene.
   Format: NCG:<TRUE|FALSE>&CancerMine:<LC|MC|HC>:num_citations
-  INTOGEN_DRIVER - Gene is predicted as a cancer driver in the `IntoGen
   Cancer Drivers Database <https://www.intogen.org/downloads>`__
-  TCGA_DRIVER - Gene is predicted as a cancer driver in the `TCGA
   pan-cancer analysis of cancer driver genes and
   mutations <https://www.ncbi.nlm.nih.gov/pubmed/29625053>`__
-  PROB_EXAC_LOF_INTOLERANT - dbNSFP_gene: the probability of being
   loss-of-function intolerant (intolerant of both heterozygous and
   homozygous lof variants) based on ExAC r0.3 data
-  PROB_EXAC_LOF_INTOLERANT_HOM - dbNSFP_gene: the probability of being
   intolerant of homozygous, but not heterozygous lof variants based on
   ExAC r0.3 data
-  PROB_EXAC_LOF_TOLERANT_NULL - dbNSFP_gene: the probability of being
   tolerant of both heterozygous and homozygous lof variants based on
   ExAC r0.3 data
-  PROB_EXAC_NONTCGA_LOF_INTOLERANT - dbNSFP_gene: the probability of
   being loss-of-function intolerant (intolerant of both heterozygous
   and homozygous lof variants) based on ExAC r0.3 nonTCGA subset
-  PROB_EXAC_NONTCGA_LOF_INTOLERANT_HOM - dbNSFP_gene: the probability
   of being intolerant of homozygous, but not heterozygous lof variants
   based on ExAC r0.3 nonTCGA subset
-  PROB_EXAC_NONTCGA_LOF_TOLERANT_NULL - dbNSFP_gene: the probability of
   being tolerant of both heterozygous and homozygous lof variants based
   on ExAC r0.3 nonTCGA subset
-  PROB_GNOMAD_LOF_INTOLERANT - dbNSFP_gene: the probability of being
   loss-of-function intolerant (intolerant of both heterozygous and
   homozygous lof variants based on gnomAD 2.1 data
-  PROB_GNOMAD_LOF_INTOLERANT_HOM - dbNSFP_gene: the probability of
   being intolerant of homozygous, but not heterozygous lof variants
   based on gnomAD 2.1 data
-  PROB_GNOMAD_LOF_TOLERANT_NULL - dbNSFP_gene: the probability of being
   tolerant of both heterozygous and homozygous lof variants based on
   gnomAD 2.1 data
-  PROB_HAPLOINSUFFICIENCY - dbNSFP_gene: Estimated probability of
   haploinsufficiency of the gene (from
   http://dx.doi.org/10.1371/journal.pgen.1001154)
-  ESSENTIAL_GENE_CRISPR - dbNSFP_gene: Essential (E) or Non-essential
   phenotype-changing (N) based on large scale CRISPR experiments. from
   http://dx.doi.org/10.1126/science.aac7041
-  ESSENTIAL_GENE_CRISPR2 - dbNSFP_gene: Essential (E), context-Specific
   essential (S), or Non-essential phenotype-changing (N) based on large
   scale CRISPR experiments. from
   http://dx.doi.org/10.1016/j.cell.2015.11.015

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
   18. `BayesDel <https://doi.org/10.1002/humu.23158>`__
   19. `LIST-S2 <https://doi.org/10.1093/nar/gkaa288>`__
   20. `ALoFT <https://www.nature.com/articles/s41467-017-00443-5>`__

-  BAYESDEL_ADDAF_DBNSFP - predicted effect from BayesDel (dbNSFP)

-  LIST_S2_DBNSFP - predicted effect from LIST-S2 (dbNSFP)

-  SIFT_DBNSFP - predicted effect from SIFT (dbNSFP)

-  SIFT4G_DBNSFP - predicted effect from SIFT4G (dbNSFP)

-  PROVEAN_DBNSFP - predicted effect from PROVEAN (dbNSFP)

-  MUTATIONTASTER_DBNSFP - predicted effect from MUTATIONTASTER (dbNSFP)

-  MUTATIONASSESSOR_DBNSFP - predicted effect from MUTATIONASSESSOR
   (dbNSFP)

-  M_CAP_DBNSFP - predicted effect from M-CAP (dbNSFP)

-  ALOFT_DBNSFP - predicted effect from ALoFT (dbNSFP)

-  MUTPRED_DBNSFP - score from MUTPRED (dbNSFP)

-  FATHMM_DBNSFP - predicted effect from FATHMM (dbNSFP)

-  PRIMATEAI_DBNSFP - predicted effect from PRIMATEAI (dbNSFP)

-  DEOGEN2_DBNSFP - predicted effect from DEOGEN2 (dbNSFP)

-  GERP_DBNSFP - evolutionary constraint measure from GERP (dbNSFP)

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
   protein (from `Open Targets
   Platform <https://www.targetvalidation.org/>`__, drugs are listed as
   `ChEMBL <https://www.ebi.ac.uk/chembl/>`__ compound identifiers)
-  CIVIC_ID, CIVIC_ID_SEGMENT - Variant/segment (exon, codon)
   identifiers in the `CIViC database <http://civic.genome.wustl.edu>`__
-  CGI_ID, CGI_ID_SEGMENT - Variant/segment (exon, codon) identifier in
   the `Cancer Genome Interpreter Cancer Biomarkers
   Database <https://www.cancergenomeinterpreter.org/biomarkers>`__

Tab-separated values (TSV)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotated List of all SNVs/InDels
'''''''''''''''''''''''''''''''''

We provide a tab-separated values file with most important annotations
for SNVs/InDels. The file has the following naming convention:

**sample_id**.pcgr_acmg.__genome_assembly__.snvs_indels.tiers.tsv

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
   18. ONCOGENE - Gene is predicted/classified as an oncogene (CancerMine/NCG)
   19. TUMOR_SUPPRESSOR - Gene is predicted/classified as tumor suppressor (CancerMine/NCG)
   20. ONCOGENE_EVIDENCE - Underlying evidence for oncogene prediction/classification
   21. TUMOR_SUPPRESSOR_EVIDENCE - Underlying evidence for tumor suppressor prediction/classification
   22. CONSEQUENCE - Variant consequence (as defined above for VCF output:
       Consequence)
   23. PROTEIN_CHANGE - Protein change (HGVSp without reference accession)
   24. PROTEIN_DOMAIN - Protein domain description (Pfam)
   25. CODING_STATUS - Coding variant status wrt. protein alteration ('coding' or 'noncoding')
   26. EXONIC_STATUS - Exonic variant status ('exonic' or 'nonexonic')
   27. CDS_CHANGE - composite VEP-based variable for coding change, format:
       Consequence:Feature:cDNA_position:EXON:HGVSp_short
   28. HGVSp
   29. HGVSc
   30. EFFECT_PREDICTIONS - as defined above for VCF
   31. MUTATION_HOTSPOT - mutation hotspot codon in
       cancerhotspots.org. Format: gene_symbol | codon | q-value
   32. MUTATION_HOTSPOT_TRANSCRIPT - hotspot-associated transcripts (Ensembl transcript ID)
   33. MUTATION_HOTSPOT_CANCERTYPE - hotspot-associated cancer types (from cancerhotspots.org)
   34. PUTATIVE_DRIVER_MUTATION - Indicates if variant is predicted as
       driver mutation from TCGA's PanCancer study of cancer driver mutation
   35. CHASMPLUS_DRIVER - Driver mutation predicted by CHASMplus algorithm
   36. CHASMPLUS_TTYPE - Tumor type for which mutation is predicted as driver by CHASMplus
   37. VEP_ALL_CSQ - all VEP transcript block consequences
   38. DBSNPRSID - dbSNP reference cluster ID
   39. COSMIC_MUTATION_ID - COSMIC mutation ID
   40. TCGA_PANCANCER_COUNT - Raw variant count across all TCGA tumor types
   41. TCGA_FREQUENCY - Frequency of variant across TCGA tumor types. Format: tumortype|
   percent affected|affected cases|total cases
   42. ICGC_PCAWG_OCCURRENCE - Mutation occurrence in ICGC-PCAWG by project:
   project_code|affected_donors|tested_donors|frequency
   43. CHEMBL_COMPOUND_ID - Compounds (as ChEMBL IDs) that are targeted towards the encoded protein (from Open Targets Platform)
   44. CHEMBL_COMPOUND_TERMS - Compounds (as drug names) that are targeted towards the encoded protein (from Open Targets Platform)
   45. SIMPLEREPEATS_HIT - Variant overlaps UCSC _simpleRepeat_ sequence repeat track
   46. WINMASKER_HIT - Variant overlaps UCSC _windowmaskerSdust_ sequence repeat track
   47. OPENTARGETS_RANK - OpenTargets association score (between 0 and 1) for gene (maximum across cancer phenotypes)
   48. CLINVAR - ClinVar association: variant origin and associated traits
   49. CLINVAR_CLNSIG - clinical significance of ClinVar variant
   50. GLOBAL_AF_GNOMAD - global germline allele frequency in gnomAD
   51. GLOBAL_AF_1KG - 1000G Project - phase 3, germline allele frequency
   52. CALL_CONFIDENCE - confidence indicator for somatic variant
   53. DP_TUMOR - sequencing depth at variant site (tumor sample)
   54. AF_TUMOR - allelic fraction of alternate allele (tumor sample)
   55. DP_CONTROL - sequencing depth at variant site (control sample)
   56. AF_CONTROL - allelic fraction of alternate allele (control sample)
   57. TIER
   58. TIER_DESCRIPTION

**NOTE**: The user has the possibility to append the TSV file with data
from other tags in the input VCF of interest (i.e. using the
*custom_tags* option in the TOML configuration file)

Mutational signature contributions
''''''''''''''''''''''''''''''''''

We provide a tab-separated values file information about mutational
signatures detected in the tumor sample. The file has the following
naming convention:

**sample_id**.pcgr_acmg.__genome_assembly__.mutational_signatures.tsv

The format of the TSV file is the following:

::

   1. signature_id - identifier for signature
   2. sample_id - sample identifier
   3. prop_signature - relative contribution of mutational signature
   4. group - keyword for signature aetiology
   5. all_reference_signatures - logical indicating if all reference signatures were used for reconstruction/inference
   6. tumor_type - tumor type (used for retrieval of reference signatures)
   7. reference_collection - collection used for reference signatures
   8. reference_signatures - signatures present in reference collection
   9. fitting_accuracy - accuracy of mutational signature fitting

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

**sample_id**.pcgr_acmg.__genome_assembly__.cna_segments.tsv.gz

**NOTE**: This is file is organized according to the *affected
transcripts* (i.e. one line/record per affected transcript).

The format of the compressed TSV file is the following:

::

   1. chrom - chromosome
   2. segment_start - start of copy number segment
   3. segment_end - end of copy number segment
   4. segment_id - segment identifier
   5. segment_length_Mb - length of segment in Mb
   6. event_type - focal or broad (covering more than 25% of chromosome arm)
   7. cytoband(s) - cytobands covered by segment
   8. LogR - Copy log-ratio
   9. sample_id - Sample identifier
   10. ensembl_gene_id
   11. symbol - gene symbol
   12. ensembl_transcript_id
   13. transcript_start - chromosomal position of transcript start
   14. transcript_end - chromosomal position of transcript end
   15. transcript_overlap_percent - percent of transcript length covered by CN segment
   16. name - gene name description
   17. biotype - type of gene
   18. tumor_suppressor - tumor suppressor gene status (CancerMine literature database/Network of Cancer Genes)
   19. oncogene - oncogene status (CancerMine literature database/Network of Cancer Genes)
   20. intogen_drivers - predicted driver gene status (IntoGen Cancer Drivers Database)
   21. chembl_compound_id - anti-cancer drugs that are targeted towards the encoded protein (from Open Targets Platform, drugs are listed as ChEMBL compound identifiers)
   22. gencode_tag
   23. gencode_release
