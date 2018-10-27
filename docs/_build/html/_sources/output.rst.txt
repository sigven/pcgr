Output
------

Interactive HTML report
~~~~~~~~~~~~~~~~~~~~~~~

An interactive and tier-structured HTML report that shows the most
relevant findings in the query cancer genome is provided with the
following naming convention:

**sample\_id**.\ **tier\_model**.\ **genome\_assembly**.html

-  The **sample\_id** is provided as input by the user, and reflects a
   unique identifier of the tumor-normal sample pair to be analyzed.
-  The **tier\_model** is provided by the user and will be either
   *'pcgr'* or *'pcgr\_acmg'*

The report is structured in seven main sections, described in more
detail below:

1. **Settings & annotation sources**

   -  Lists underlying tools and annotation sources (links and versions)
   -  Lists key configurations provided by user

2. **Main results**

   -  Six value boxes that highlight the main findings of clinical
      relevance in the tumor:

      1. Mutational signatures - two most prevalent signatures (other
         than aging)
      2. Tier 1 variants (top four)
      3. Tier 2 variants (top four)
      4. Tumor mutational burden
      5. Microsatellite instability prediction
      6. Somatic copy number aberrations of clinical significance

3. **Somatic SNVs/InDels**

   -  *Mutational burden (TMB)*

      -  given a coding target region size specified by the user
         (ideally the **callable target size**), an estimate of the
         mutational burden is provided
      -  is presently only computed for tumor-normal input (e.g.
         *vcf\_tumor\_only = false*)
      -  The estimated mutational burden is assigned a descriptive
         *tertile* based on thresholds defined by the user (these should
         reflect thresholds of clinical significance, and may vary for
         different tumor types)

   -  *Variant & tier statistics*

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

      -  Variants are organized into four tiers (interactive datatables)
         according to clinical utility
      -  Users can choose between two `tier
         models <tier_systems.html>`__:

         -  the original model (*pcgr*) that do not take into account
            tumor type of input when assigning variants to (top) tiers
         -  the new model (*pcgr\_acmg*) that takes into account tumor
            type of input and strength of clinical evidence when
            assigning variants to (top) tiers

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
         acting as biomarkers across therapeutic contexts, tumor types,
         evidence levels etc.

5. **MSI status**

   -  Indicates predicted microsatellite stability from the somatic
      mutation profile and supporting evidence (details of the
      underlying MSI statistical classifier can be found
      `here <http://rpubs.com/sigven/msi2018>`__)
   -  The MSI classifier was trained on TCGA exome samples.
   -  Will only be present in the report if specified by the user in the
      configuration file (i.e. *msi = true*) and if the input is
      tumor-normal (i.e. *vcf\_tumor\_only = false*)

6. **Mutational signatures**

   -  Estimation of relative contribution of `30 known mutational
      signatures <http://cancer.sanger.ac.uk/cosmic/signatures>`__ in
      tumor sample (using
      `deconstructSigs <https://github.com/raerose01/deconstructSigs>`__
      as the underlying framework)
   -  Datatable with signatures and proposed underlying etiologies
   -  Will only be present in the report if specified by the user in the
      configuration file (i.e. *mutsignatures = true*) and if the input
      is tumor-normal (i.e. *vcf\_tumor\_only = false*)
   -  `Trimer (i.e. DNA 3-mer)
      normalization <https://github.com/raerose01/deconstructSigs>`__
      can be configured according to sequencing approach used (WES, WXS
      etc.) using the 'mutsignatures\_normalization' option, as can the

      -  minimum number of mutations required for analysis (option
         'mutsignatures\_mutation\_limit')
      -  the maximum number of mutational signatures in the search space
         (option 'mutsignatures\_signature\_limit')
      -  possibility to discard any signature contributions with a
         weight less than a given cutoff (option
         'mutsignatures\_cutoff')

7. **References**

   -  Supporting scientific literature (key report elements)

Interactive datatables
^^^^^^^^^^^^^^^^^^^^^^

The interactive datatables contain a number of hyperlinked annotations
similar to those defined for the annotated VCF file, including the
following:

-  SYMBOL - Gene symbol (Entrez/NCBI)
-  PROTEIN\_CHANGE - Amino acid change (VEP)
-  CANCER\_TYPE - Biomarker (tier 1/2): associated cancer type
-  EVIDENCE\_LEVEL - Biomarker (tier 1/2): evidence level (A,B,C,D,E)
-  CLINICAL\_SIGNIFICANCE - Biomarker (tier 1/2): drug sensitivity,
   poor/better outcome etc
-  EVIDENCE\_TYPE - Biomarker (tier 1/2):
   predictive/diagnostic/therapeutic
-  DISEASE\_ONTOLOGY\_ID - Biomarker (tier 1/2): associated cancer type
   (Disease Ontology)
-  EVIDENCE\_DIRECTION - Biomarker (tier 1/2): supports/does not support
-  DESCRIPTION - Biomarker (tier 1/2): description
-  VARIANT\_ORIGIN - Biomarker (tier 1/2): variant origin
   (germline/somatic)
-  BIOMARKER\_MAPPING - Biomarker (tier 1/2): accuracy of genomic
   mapping (exact,codon,exon)
-  CITATION - Biomarker (tier 1/2): supporting literature
-  THERAPEUTIC\_CONTEXT - Biomarker (tier 1/2): associated drugs
-  RATING - Biomarker (tier 1/2): trust rating from 1 to 5 (CIVIC)
-  GENE\_NAME - gene name description (Entrez/NCBI)
-  PROTEIN\_DOMAIN - PFAM protein domain
-  PROTEIN\_FEATURE - UniProt feature overlapping variant site
-  CDS\_CHANGE - Coding sequence change
-  MUTATION\_HOTSPOT - Known cancer mutation hotspot
-  MUTATION\_HOTSPOT\_CANCERTYPE - Hotspot-associated cancer types
-  TCGA\_FREQUENCY - Frequency of variant in TCGA cohorts
-  ICGC\_PCAWG\_OCCURRENCE - Frequency of variant in ICGC-PCAWG cohorts
-  DOCM\_LITERATURE - Literature links - DoCM
-  DOCM\_DISEASE - Associated diseases - DoCM
-  INTOGEN\_DRIVER\_MUT - predicted driver mutation - IntOGen
-  CONSEQUENCE - VEP consequence (primary transcript)
-  HGVSc - from VEP
-  HGVSp - from VEP
-  ONCOGENE - Predicted as proto-oncogene from literature mining
-  TUMOR\_SUPPRESSOR - Predicted as tumor suppressor gene from
   literature mining
-  ONCOSCORE - Literature-derived score for oncogenic potential (gene
   level)
-  PREDICTED\_EFFECT - Effect predictions from dbNSFP
-  VEP\_ALL\_CONSEQUENCE - All VEP consequences (multiple transcripts)
-  DBSNP - dbSNP rsID
-  COSMIC - Cosmic mutation IDs
-  CLINVAR - ClinVar variant origin and associated phenotypes
-  CANCER\_ASSOCIATIONS - Gene-associated cancer types from DisGenet
-  TARGETED\_DRUGS - Targeted drugs from DGIdb-ChEMBL
-  KEGG\_PATHWAY - Gene-associated pathways from KEGG
-  CALL\_CONFIDENCE - Variant confidence (as set by user in input VCF)
-  DP\_TUMOR - Variant sequencing depth in tumor (as set by user in
   input VCF)
-  AF\_TUMOR - Variant allelic fraction in tumor (as set by user in
   input VCF)
-  DP\_NORMAL - Variant sequencing depth in normal (as set by user in
   input VCF)
-  AF\_NORMAL - Variant allelic fraction in tumor (as set by user in
   input VCF)
-  GENOMIC\_CHANGE - Variant ID
-  GENOME\_VERSION - Genome assembly

Example reports:

-  `View an example report for a breast tumor sample
   (TCGA) <http://folk.uio.no/sigven/tumor_sample.BRCA.pcgr_acmg.grch37.0.6.2.html>`__
-  `View an example report for a colon adenocarcinoma sample
   (TCGA) <http://folk.uio.no/sigven/tumor_sample.COAD.pcgr_acmg.grch37.0.6.2.html>`__

The HTML reports have been tested using the following browsers:

-  Safari (10.0.3)
-  Mozilla Firefox (52.0.2)
-  Google Chrome (57.0.2987.110)

JSON (beta)
~~~~~~~~~~~

A JSON file that stores the HTML report content is provided. This file
will easen the process of extracting particular parts of the report for
further analysis. Presently, there is no detailed schema documented for
the PCGR JSON structure. Examples (using R) on how to extract
information from the JSON file will soon be posted here.

Output files - somatic SNVs/InDels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Variant call format - VCF
^^^^^^^^^^^^^^^^^^^^^^^^^

A VCF file containing annotated, somatic calls (single nucleotide
variants and insertion/deletions) is generated with the following naming
convention:

**sample\_id**.\ **tier\_model**.\ **genome\_assembly**.vcf.gz

Here, the **sample\_id** is provided as input by the user, and reflects
a unique identifier of the tumor-normal sample pair to be analyzed.
Following common standards, the annotated VCF file is compressed with
`bgzip <http://www.htslib.org/doc/tabix.html>`__ and indexed with
`tabix <http://www.htslib.org/doc/tabix.html>`__. Below follows a
description of all annotations/tags present in the VCF INFO column after
processing with the PCGR annotation pipeline:

*VEP consequence annotations*
'''''''''''''''''''''''''''''

-  CSQ - Complete consequence annotations from VEP. Format:
   Allele\|Consequence\|IMPACT\|SYMBOL\|Gene\|Feature\_type\|Feature\|BIOTYPE\|EXON\|
   INTRON\|HGVSc\|HGVSp\|cDNA\_position\|CDS\_position\|Protein\_position\|Amino\_acids\|
   Codons\|Existing\_variation\|ALLELE\_NUM\|DISTANCE\|STRAND\|FLAGS\|PICK\|VARIANT\_CLASS\|
   SYMBOL\_SOURCE\|HGNC\_ID\|CANONICAL\|APPRIS\|CCDS\|ENSP\|SWISSPROT\|TREMBL\|UNIPARC\|
   RefSeq\|DOMAINS\|HGVS\_OFFSET\|AF\|AFR\_AF\|AMR\_AF\|EAS\_AF\|EUR\_AF\|SAS\_AF\|gnomAD\_AF\|
   gnomAD\_AFR\_AF\|gnomAD\_AMR\_AF\|gnomAD\_ASJ\_AF\|gnomAD\_EAS\_AF\|gnomAD\_FIN\_AF\|
   gnomAD\_NFE\_AF\|gnomAD\_OTH\_AF\|gnomAD\_SAS\_AF\|CLIN\_SIG\|SOMATIC\|PHENO\|
   MOTIF\_NAME\|MOTIF\_POS\|HIGH\_INF\_POS\|MOTIF\_SCORE\_CHANGE
-  Consequence - Impact modifier for the consequence type (picked by
   VEP's --flag\_pick\_allele option)
-  Gene - Ensembl stable ID of affected gene (picked by VEP's
   --flag\_pick\_allele option)
-  Feature\_type - Type of feature. Currently one of Transcript,
   RegulatoryFeature, MotifFeature (picked by VEP's --flag\_pick\_allele
   option)
-  Feature - Ensembl stable ID of feature (picked by VEP's
   --flag\_pick\_allele option)
-  cDNA\_position - Relative position of base pair in cDNA sequence
   (picked by VEP's --flag\_pick\_allele option)
-  CDS\_position - Relative position of base pair in coding sequence
   (picked by VEP's --flag\_pick\_allele option)
-  CDS\_CHANGE - Coding, transcript-specific sequence annotation (picked
   by VEP's --flag\_pick\_allele option)
-  AMINO\_ACID\_START - Protein position indicating absolute start of
   amino acid altered (fetched from Protein\_position)
-  AMINO\_ACID\_END - Protein position indicating absolute end of amino
   acid altered (fetched from Protein\_position)
-  Protein\_position - Relative position of amino acid in protein
   (picked by VEP's --flag\_pick\_allele option)
-  Amino\_acids - Only given if the variant affects the protein-coding
   sequence (picked by VEP's --flag\_pick\_allele option)
-  Codons - The alternative codons with the variant base in upper case
   (picked by VEP's --flag\_pick\_allele option)
-  IMPACT - Impact modifier for the consequence type (picked by VEP's
   --flag\_pick\_allele option)
-  VARIANT\_CLASS - Sequence Ontology variant class (picked by VEP's
   --flag\_pick\_allele option)
-  SYMBOL - Gene symbol (picked by VEP's --flag\_pick\_allele option)
-  SYMBOL\_SOURCE - The source of the gene symbol (picked by VEP's
   --flag\_pick\_allele option)
-  STRAND - The DNA strand (1 or -1) on which the transcript/feature
   lies (picked by VEP's --flag\_pick\_allele option)
-  ENSP - The Ensembl protein identifier of the affected transcript
   (picked by VEP's --flag\_pick\_allele option)
-  FLAGS - Transcript quality flags: cds\_start\_NF: CDS 5', incomplete
   cds\_end\_NF: CDS 3' incomplete (picked by VEP's --flag\_pick\_allele
   option)
-  SWISSPROT - Best match UniProtKB/Swiss-Prot accession of protein
   product (picked by VEP's --flag\_pick\_allele option)
-  TREMBL - Best match UniProtKB/TrEMBL accession of protein product
   (picked by VEP's --flag\_pick\_allele option)
-  UNIPARC - Best match UniParc accession of protein product (picked by
   VEP's --flag\_pick\_allele option)
-  HGVSc - The HGVS coding sequence name (picked by VEP's
   --flag\_pick\_allele option)
-  HGVSp - The HGVS protein sequence name (picked by VEP's
   --flag\_pick\_allele option)
-  HGVSp\_short - The HGVS protein sequence name, short version (picked
   by VEP's --flag\_pick\_allele option)
-  HGVS\_OFFSET - Indicates by how many bases the HGVS notations for
   this variant have been shifted (picked by VEP's --flag\_pick\_allele
   option)
-  MOTIF\_NAME - The source and identifier of a transcription factor
   binding profile aligned at this position (picked by VEP's
   --flag\_pick\_allele option)
-  MOTIF\_POS - The relative position of the variation in the aligned
   TFBP (picked by VEP's --flag\_pick\_allele option)
-  HIGH\_INF\_POS - A flag indicating if the variant falls in a high
   information position of a transcription factor binding profile (TFBP)
   (picked by VEP's --flag\_pick\_allele option)
-  MOTIF\_SCORE\_CHANGE - The difference in motif score of the reference
   and variant sequences for the TFBP (picked by VEP's
   --flag\_pick\_allele option)
-  CELL\_TYPE - List of cell types and classifications for regulatory
   feature (picked by VEP's --flag\_pick\_allele option)
-  CANONICAL - A flag indicating if the transcript is denoted as the
   canonical transcript for this gene (picked by VEP's
   --flag\_pick\_allele option)
-  CCDS - The CCDS identifier for this transcript, where applicable
   (picked by VEP's --flag\_pick\_allele option)
-  INTRON - The intron number (out of total number) (picked by VEP's
   --flag\_pick\_allele option)
-  EXON - The exon number (out of total number) (picked by VEP's
   --flag\_pick\_allele option)
-  DOMAINS - The source and identifier of any overlapping protein
   domains (picked by VEP's --flag\_pick\_allele option)
-  DISTANCE - Shortest distance from variant to transcript (picked by
   VEP's --flag\_pick\_allele option)
-  BIOTYPE - Biotype of transcript or regulatory feature (picked by
   VEP's --flag\_pick\_allele option)
-  TSL - Transcript support level (picked by VEP's --flag\_pick\_allele
   option)>
-  PUBMED - PubMed ID(s) of publications that cite existing variant -
   VEP
-  PHENO - Indicates if existing variant is associated with a phenotype,
   disease or trait - VEP
-  GENE\_PHENO - Indicates if overlapped gene is associated with a
   phenotype, disease or trait - VEP
-  ALLELE\_NUM - Allele number from input; 0 is reference, 1 is first
   alternate etc - VEP
-  REFSEQ\_MATCH - The RefSeq transcript match status; contains a number
   of flags indicating whether this RefSeq transcript matches the
   underlying reference sequence and/or an Ensembl transcript (picked by
   VEP's --flag\_pick\_allele option)
-  PICK - Indicates if this block of consequence data was picked by
   VEP's --flag\_pick\_allele option
-  VEP\_ALL\_CONSEQUENCE - All transcript consequences
   (Consequence:SYMBOL:Feature\_type:Feature:BIOTYPE) - VEP

*Gene information*
''''''''''''''''''

-  ENTREZ\_ID - `Entrez <http://www.ncbi.nlm.nih.gov/gene>`__ gene
   identifier
-  APPRIS - Principal isoform flags according to the `APPRIS principal
   isoform database <http://appris.bioinfo.cnio.es/#/downloads>`__
-  UNIPROT\_ID - `UniProt <http://www.uniprot.org>`__ identifier
-  UNIPROT\_ACC - `UniProt <http://www.uniprot.org>`__ accession(s)
-  ENSEMBL\_GENE\_ID - Ensembl gene identifier for VEP's picked
   transcript (*ENSGXXXXXXX*)
-  ENSEMBL\_TRANSCRIPT\_ID - Ensembl transcript identifier for VEP's
   picked transcript (*ENSTXXXXXX*)
-  REFSEQ\_MRNA - Corresponding RefSeq transcript(s) identifier for
   VEP's picked transcript (*NM\_XXXXX*)
-  CORUM\_ID - Associated protein complexes (identifiers) from
   `CORUM <http://mips.helmholtz-muenchen.de/corum/>`__
-  DISGENET\_CUI - Tumor types associated with gene, as found in
   DisGeNET. Tumor types are listed as unique
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept IDs
   (*CUIs*)
-  TUMOR\_SUPPRESSOR - Gene is predicted as tumor suppressor candidate
   according to (`TSGene v2.0 <https://bioinfo.uth.edu/TSGene/>`__)
-  ONCOGENE - Gene is curated as an oncogene according to (`TSGene
   v2.0 <https://bioinfo.uth.edu/TSGene/>`__)
-  ONCOSCORE - Literature-derived score for cancer gene relevance
   `Bioconductor/OncoScore <http://bioconductor.org/packages/release/bioc/html/OncoScore.html>`__,
   range from 0 (low oncogenic potential) to 1 (high oncogenic
   potential)
-  INTOGEN\_DRIVER - Gene is predicted as a cancer driver in the
   `IntoGen Cancer Drivers Database -
   2014.12 <https://www.intogen.org/downloads>`__
-  TCGA\_DRIVER - Gene is predicted as a cancer driver in the
   `Pan-cancer analysis of cancer driver
   genes <https://www.ncbi.nlm.nih.gov/pubmed/29625053>`__

*Variant effect and protein-coding information*
'''''''''''''''''''''''''''''''''''''''''''''''

-  MUTATION\_HOTSPOT - mutation hotspot codon in
   `cancerhotspots.org <http://cancerhotspots.org/>`__. Format:
   gene\_symbol \| codon \| q-value
-  MUTATION\_HOTSPOT\_TRANSCRIPT - hotspot-associated transcripts
   (Ensembl transcript ID)
-  MUTATION\_HOTSPOT\_CANCERTYPE - hotspot-associated cancer types (from
   cancerhotspots.org)
-  UNIPROT\_FEATURE - Overlapping protein annotations from `UniProt
   KB <http://www.uniprot.org>`__
-  PFAM\_DOMAIN - Pfam domain identifier (from VEP)
-  INTOGEN\_DRIVER\_MUT - Indicates if existing variant is predicted as
   driver mutation from IntoGen Catalog of Driver Mutations
-  EFFECT\_PREDICTIONS - All predictions of effect of variant on protein
   function and pre-mRNA splicing from `database of non-synonymous
   functional predictions - dbNSFP
   v3.5 <https://sites.google.com/site/jpopgen/dbNSFP>`__. Predicted
   effects are provided by different sources/algorithms (separated by
   '&'):

   1.  `SIFT <http://provean.jcvi.org/index.php>`__ (Jan 2015)
   2.  `LRT <http://www.genetics.wustl.edu/jflab/lrt_query.html>`__
       (2009)
   3.  `MutationTaster <http://www.mutationtaster.org/>`__ (data release
       Nov 2015)
   4.  `MutationAssessor <http://mutationassessor.org/>`__ (release 3)
   5.  `FATHMM <http://fathmm.biocompute.org.uk>`__ (v2.3)
   6.  `PROVEAN <http://provean.jcvi.org/index.php>`__ (v1.1 Jan 2015)
   7.  `FATHMM\_MKL <http://fathmm.biocompute.org.uk/fathmmMKL.htm>`__
   8.  `DBNSFP\_CONSENSUS\_SVM <https://www.ncbi.nlm.nih.gov/pubmed/25552646>`__
       (Ensembl/consensus prediction, based on support vector machines)
   9.  `DBNSFP\_CONSENSUS\_LR <https://www.ncbi.nlm.nih.gov/pubmed/25552646>`__
       (Ensembl/consensus prediction, logistic regression based)
   10. `SPLICE\_SITE\_EFFECT\_ADA <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       adaptive boosting)
   11. `SPLICE\_SITE\_EFFECT\_RF <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       random forest)
   12. `M-CAP <http://bejerano.stanford.edu/MCAP>`__
   13. `MutPred <http://mutpred.mutdb.org>`__
   14. `GERP <http://mendel.stanford.edu/SidowLab/downloads/gerp/>`__

-  SIFT\_DBNSFP - predicted effect from SIFT (dbNSFP)
-  PROVEAN\_DBNSFP - predicted effect from PROVEAN (dbNSFP)
-  MUTATIONTASTER\_DBNSFP - predicted effect from MUTATIONTASTER
   (dbNSFP)
-  MUTATIONASSESSOR\_DBNSFP - predicted effect from MUTATIONASSESSOR
   (dbNSFP)
-  M\_CAP\_DBNSFP - predicted effect from M-CAP (dbNSFP)
-  MUTPRED\_DBNSFP - score from MUTPRED (dbNSFP)
-  FATHMM\_DBNSFP - predicted effect from FATHMM (dbNSFP)
-  FATHMM\_MKL\_DBNSFP - predicted effect from FATHMM-mkl (dbNSFP)
-  META\_LR\_DBNSFP - predicted effect from ensemble prediction
   (logistic regression - dbNSFP)
-  SPLICE\_SITE\_RF\_DBNSFP - predicted effect of splice site
   disruption, using random forest (dbscSNV)
-  SPLICE\_SITE\_ADA\_DBNSFP - predicted effect of splice site
   disruption, using boosting (dbscSNV)

*Variant frequencies/annotations in germline/somatic databases*
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

-  AFR\_AF\_GNOMAD - African/American germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  AMR\_AF\_GNOMAD - American germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  GLOBAL\_AF\_GNOMAD - Adjusted global germline allele frequency
   (`Genome Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  SAS\_AF\_GNOMAD - South Asian germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  EAS\_AF\_GNOMAD - East Asian germline allele frequency (`Genome
   Aggregation Database release
   21 <http://gnomad.broadinstitute.org/>`__)
-  FIN\_AF\_GNOMAD - Finnish germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  NFE\_AF\_GNOMAD - Non-Finnish European germline allele frequency
   (`Genome Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  OTH\_AF\_GNOMAD - Other germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  ASJ\_AF\_GNOMAD - Ashkenazi Jewish allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  AFR\_AF\_1KG - `1000G Project - phase
   3 <http://www.1000genomes.org>`__ germline allele frequency for
   samples from AFR (African)
-  AMR\_AF\_1KG - `1000G Project - phase
   3 <http://www.1000genomes.org>`__ germline allele frequency for
   samples from AMR (Ad Mixed American)
-  EAS\_AF\_1KG - `1000G Project - phase
   3 <http://www.1000genomes.org>`__ germline allele frequency for
   samples from EAS (East Asian)
-  EUR\_AF\_1KG - `1000G Project - phase
   3 <http://www.1000genomes.org>`__ germline allele frequency for
   samples from EUR (European)
-  SAS\_AF\_1KG - `1000G Project - phase
   3 <http://www.1000genomes.org>`__ germline allele frequency for
   samples from SAS (South Asian)
-  GLOBAL\_AF\_1KG - `1000G Project - phase
   3 <http://www.1000genomes.org>`__ germline allele frequency for all
   1000G project samples (global)
-  DBSNPRSID - `dbSNP <http://www.ncbi.nlm.nih.gov/SNP/>`__ reference
   ID, as provided by VEP
-  COSMIC\_MUTATION\_ID - Mutation identifier in `Catalog of somatic
   mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__
   database, as provided by VEP
-  TCGA\_PANCANCER\_COUNT - Raw variant count across all TCGA tumor
   types
-  TCGA\_FREQUENCY - Frequency of variant across TCGA tumor types.
   Format: tumortype\| percent affected\|affected cases\|total cases
-  ICGC\_PCAWG\_OCCURRENCE - Mutation occurrence in
   `ICGC-PCAWG <http://docs.icgc.org/pcawg/>`__. By project:
   project\_code\|affected\_donors\|tested\_donors\|frequency)
-  ICGC\_PCAWG\_AFFECTED\_DONORS - Number of donors with the current
   mutation in `ICGC-PCAWG <http://docs.icgc.org/pcawg/>`__

*Clinical associations*
'''''''''''''''''''''''

-  CLINVAR\_MSID - `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
   Measure Set/Variant ID
-  CLINVAR\_ALLELE\_ID -
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ allele ID
-  CLINVAR\_PMID - Associated Pubmed IDs for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - germline
   state-of-origin
-  CLINVAR\_HGVSP - Protein variant expression using HGVS nomenclature
-  CLINVAR\_PMID\_SOMATIC - Associated Pubmed IDs for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - somatic
   state-of-origin
-  CLINVAR\_CLNSIG - Clinical significance for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - germline
   state-of-origin
-  CLINVAR\_CLNSIG\_SOMATIC - Clinical significance for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - somatic
   state-of-origin
-  CLINVAR\_MEDGEN\_CUI - Associated
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept identifiers
   (*CUIs*) - germline state-of-origin
-  CLINVAR\_MEDGEN\_CUI\_SOMATIC - Associated
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept identifiers
   (*CUIs*) - somatic state-of-origin
-  CLINVAR\_VARIANT\_ORIGIN - Origin of variant (somatic, germline, de
   novo etc.) for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
-  DOCM\_PMID - Associated Pubmed IDs for variant in `Database of
   Curated Mutations <http://docm.genome.wustl.edu>`__

*Other*
'''''''

-  CHEMBL\_COMPOUND\_ID - antineoplastic drugs targeting the encoded
   protein (from `Drug-Gene Interaction
   Database <http://dgidb.genome.wustl.edu/>`__, drugs are listed as
   `ChEMBL <https://www.ebi.ac.uk/chembl/>`__ compound identifiers)
-  CIVIC\_ID, CIVIC\_ID\_2 - Variant identifiers in the `CIViC
   database <http://civic.genome.wustl.edu>`__, CIVIC\_ID refers to
   markers mapped at variant level, CIVIC\_ID\_2 refers to region
   markers (codon, exon etc.)
-  CBMDB\_ID - Variant identifier in the `Cancer Biomarkers
   database <https://www.cancergenomeinterpreter.org/biomarkers>`__

Tab-separated values (TSV)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotated List of all SNVs/InDels
'''''''''''''''''''''''''''''''''

We provide a tab-separated values file with most important annotations
for SNVs/InDels. The file has the following naming convention:

**sample\_id**.\ **tier\_model**.\ **genome\_assembly**.snvs\_indels.tiers.tsv

The SNVs/InDels are organized into different **tiers** (as defined above
for the HTML report)

The following variables are included in the tiered TSV file:

::

    1. GENOMIC_CHANGE - Identifier for variant at the genome (VCF) level, e.g. 1:g.152382569A>G
          Format: (<chrom>:g.<position><ref_allele>><alt_allele>)
    2. GENOME_VERSION - Assembly version, e.g. GRCh37
    3. VCF_SAMPLE_ID - Sample identifier
    4. VARIANT_CLASS - Variant type, e.g. SNV/insertion/deletion
    5. SYMBOL - Gene symbol
    6. GENE_NAME - Gene description
    7. CCDS - CCDS identifier
    8. ENTREZ_ID - Entrez gene identifier
    9. UNIPROT_ID - UniProt protein identifier
    10. ENSEMBL_TRANSCRIPT_ID - Ensembl transcript identifier
    11. ENSEMBL_GENE_ID - Ensembl gene identifier
    12. REFSEQ_MRNA - RefSeq mRNA identifier
    13. ONCOSCORE - Literature-derived score for cancer gene relevance
    14. ONCOGENE - Gene is predicted as an oncogene according to literature mining (CancerMine)
    15. TUMOR_SUPPRESSOR - Gene is predicted as tumor suppressor according to literature mining(CancerMine)
        candidate according to TSGene
    16. DISGENET_CUI - Associated tumor types from DisGeNET (MedGen concept IDs)
    17. DISGENET_TERMS - Associated tumor types from DisGeNET (MedGen concept terms)
    18. CONSEQUENCE - Variant consequence (as defined above for VCF output:
        Consequence)
    19. PROTEIN_CHANGE - Protein change (HGVSp without reference accession)
    20. PROTEIN_DOMAIN - Protein domain
    21. CDS_CHANGE - composite VEP-based variable for coding change, format:
        Consequence:Feature:cDNA_position:EXON:HGVSp_short
    22. HGVSp
    23. HGVSc
    24. EFFECT_PREDICTIONS - as defined above for VCF
    25. MUTATION_HOTSPOT - mutation hotspot codon in
        cancerhotspots.org. Format: gene_symbol | codon | q-value
    26. MUTATION_HOTSPOT_TRANSCRIPT - hotspot-associated transcripts (Ensembl transcript ID)
    27. MUTATION_HOTSPOT_CANCERTYPE - hotspot-associated cancer types (from cancerhotspots.org)
    28. INTOGEN_DRIVER_MUT - Indicates if existing variant is predicted as
        driver mutation from IntoGen Catalog of Driver Mutations
    29. VEP_ALL_CONSEQUENCE - all VEP consequences
    30. DBSNPRSID - dbSNP reference cluster ID
    31. COSMIC_MUTATION_ID - COSMIC mutation ID
    32. TCGA_PANCANCER_COUNT - Raw variant count across all TCGA tumor types
    33. TCGA_FREQUENCY - Frequency of variant across TCGA tumor types. Format: tumortype|
    percent affected|affected cases|total cases
    34. ICGC_PCAWG_OCCURRENCE - Mutation occurrence in ICGC-PCAWG by project:
    project_code|affected_donors|tested_donors|frequency
    35. CHEMBL_COMPOUND_ID - Compounds (as ChEMBL IDs) that target the encoded protein (from DGIdb)
    36. CHEMBL_COMPOUND_TERMS - Compounds (as drug names) that target the encoded protein (from DGIdb)
    37. CLINVAR - ClinVar association: variant origin and associated traits
    38. CLINVAR_CLNSIG - clinical significance of ClinVar variant
    39. GLOBAL_AF_GNOMAD - global germline allele frequency in gnomAD
    40. GLOBAL_AF_1KG - 1000G Project - phase 3, germline allele frequency
    41. CALL_CONFIDENCE - confidence indicator for somatic variant
    42. DP_TUMOR - sequencing depth at variant site (tumor)
    43. AF_TUMOR - allelic fraction of alternate allele (tumor)
    44. DP_NORMAL - sequencing depth at variant site (normal)
    45. AF_NORMAL - allelic fraction of alternate allele (normal)
    46. TIER
    47. TIER_DESCRIPTION

**NOTE**: The user has the possibility to append the TSV file with data
from other tags in the input VCF of interest (i.e. using the
*custom\_tags* option in the TOML configuration file)

Output files - somatic copy number aberrations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Tab-separated values (TSV)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copy number segments are intersected with the genomic coordinates of all
transcripts from `GENCODE's basic gene
annotation <https://www.gencodegenes.org/releases/current.html>`__. In
addition, PCGR attaches cancer-relevant annotations for the affected
transcripts. The naming convention of the compressed TSV file is as
follows:

**sample\_id**.\ **tier\_model**.\ **genome\_assembly**.cna\_segments.tsv.gz

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
    17. tsgene - tumor suppressor gene status (TSgene database)
    18. tsgene_oncogene - oncogene status (TSgene database)
    19. intogen_drivers - predicted driver gene status (IntoGen Cancer Drivers Database)
    20. chembl_compound_id - antineoplastic drugs targeting the encoded protein
       (from DGIdb, drugs are listed as ChEMBL compound identifiers)
    21. gencode_gene_biotype
    22. gencode_tag
    23. gencode_release
