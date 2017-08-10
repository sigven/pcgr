Input & output
--------------

Input
~~~~~

The PCGR workflow accepts two types of input files:

-  An unannotated, single-sample VCF file (>= v4.2) with called somatic
   variants (SNVs/InDels)
-  A copy number segment file

**IMPORTANT NOTE: GRCh37 is the reference genome build currently
supported by PCGR**

PCGR can be run with either or both of the two input files present.

VCF
^^^

-  We **strongly** recommend that the input VCF is compressed and
   indexed using `bgzip <http://www.htslib.org/doc/tabix.html>`__ and
   `tabix <http://www.htslib.org/doc/tabix.html>`__
-  If the input VCF contains multi-allelic sites, these will be subject
   to `decomposition <http://genome.sph.umich.edu/wiki/Vt#Decompose>`__
-  Variants used for reporting should be designated as 'PASS' in the VCF
   FILTER column

**IMPORTANT NOTE 1**: Considering the VCF output for the `numerous
somatic SNV/InDel callers <https://www.biostars.org/p/19104/>`__ that
have been developed, we have a experienced a general lack of uniformity
and robustness for the representation of somatic variant genotype data
(e.g. variant allelic depths (tumor/normal), genotype quality etc.).
Variant genotype data can be specified as optional arguments top the
PCGR workflow, which in turn will be used for filtering and output in
the tumor report.

**IMPORTANT NOTE 2**: PCGR generates a number of VCF INFO annotation
tags that is appended to the query VCF. We will therefore encourage the
users to submit query VCF files that have not been subject to
annotations by other means, but rather a VCF file that comes directly
from variant calling. If not, there are likely to be INFO tags in the
query VCF file that coincide with those produced by PCGR.

Copy number segments
^^^^^^^^^^^^^^^^^^^^

The tab-separated values file with copy number aberrations **MUST**
contain the following four columns:

-  Chromosome
-  Start
-  End
-  Segment\_Mean

Here, *Chromosome*, *Start*, and *End* denote the chromosomal segment
(GRCh37), and **Segment\_Mean** denotes the log(2) ratio for a
particular segment, which is a common output of somatic copy number
alteration callers. Below shows the initial part of a copy number
segment file that is formatted correctly according to PCGR's
requirements:

::

      Chromosome    Start   End Segment_Mean
      1 3218329 3550598 0.0024
      1 3552451 4593614 0.1995
      1 4593663 6433129 -1.0277

Output - Interactive HTML report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An interactive and tier-structured HTML report that shows the most
relevant findings in the query cancer genome is provided with the
following naming convention:

**sample\_id**.pcgr.html

The **sample\_id** is provided as input by the user, and reflects a
unique identifier of the tumor-normal sample pair to be analyzed.

The report is structured in six main sections, described in more detail
below:

1. **Annotation sources**

   -  Lists underlying tools and annotation sources (versions)

2. **Somatic SNVs/InDels**

   -  *Summary statistics* - indicate number of SNVs/InDels as well as
      number of coding/non-coding variants
   -  *Tier statistics* - indicate number of variant found in each tier
      (see below)
   -  *Global distribution - allelic support* - distribution (histogram)
      of variant allelic support for somatic variants (will only be
      present in the report if specific fields in input VCF is defined
      and specified by the user)
   -  *Global variant browser* - permits filtering of the whole
      SNV/Indel dataset by various criteria (call confidence, variant
      sequencing depth, variant consequence etc.)
   -  Variants are organized into five tiers (interactive datatables)
      according to clinical utility

      -  *Tier 1* - constitutes variants linked to predictive,
         prognostic, diagnostic, and predisposision biomarkers in the
         `CIViC database <http://civic.genome.wustl.edu>`__ and the
         `Cancer Biomarkers
         Database <https://www.cancergenomeinterpreter.org/biomarkers>`__
      -  *Tier 2* - includes other coding variants that are found in
         known cancer mutation hotspots, predicted as cancer driver
         mutations, or curated as disease-causing
      -  *Tier 3* - includes other coding variants found in oncogenes,
         tumor suppressor genes, or cancer census genes
      -  *Tier 4* - includes other coding variants
      -  *Tier 5* - includes non-coding variants

         -  will only present if specified by the user
            ('--list\_noncoding')

      -  **NOTE**: The tier structure is inspired by recommended variant
         prioritization by `Dienstmann et al.,
         2014 <https://www.ncbi.nlm.nih.gov/pubmed/24768039>`__. The
         table below shows the correspondence between the terminology
         for reportable variants used by Dienstmann et al. and the tiers
         in PCGR:

      +-----------------------------------------------------------------------------------------+-------------------+
      | `Dienstmann et al., 2014 <https://www.ncbi.nlm.nih.gov/pubmed/24768039>`__ (Figure 2)   | PCGR              |
      +=========================================================================================+===================+
      | *Actionable*                                                                            | Tier 1            |
      +-----------------------------------------------------------------------------------------+-------------------+
      | *Other relevant variants*                                                               | Tier 2 & Tier 3   |
      +-----------------------------------------------------------------------------------------+-------------------+
      | *Unknown*\ \*                                                                           | Tier 4 (Tier 5)   |
      +-----------------------------------------------------------------------------------------+-------------------+

      \*While Dienstmann et al. suggest that the *Unknown* category
      should be categorized according to pathways, the PCGR employs an
      arrangement of variant results according to genes, using a
      literature-derived score of oncogenic potential (KEGG pathway
      information is also linked).

3. **Somatic CNA analysis**

   -  *Segments - amplifications and homozygous deletions*

      -  Based on user-defined/default log-ratio thresholds of
         gains/losses, the whole CNA dataset can be navigated further
         through filters (e.g. cytoband and type of event (focal or
         broad))

   -  *Proto-oncogenes subject to copy number amplifications*

      -  Datatable listing known proto-oncogenes covered by
         user-defined/default amplifications

   -  *Tumor suppressor genes subject to homozygous deletions*

      -  Datatable listing known tumor suppressor genes covered by
         user-defined/default losses

   -  *Copy number aberrations as biomarkers for prognosis, diagnosis,
      predisposition, and drug response*

      -  Interactive data table where the user can navigate aberrations
         acting as biomarkers across therapeutic contexts, tumor types,
         evidence level etc.

4. **MSI status**

   -  Indicates predicted microsatellite stability from the somatic
      mutation profile and supporting evidence (details of the
      underlying MSI statistical classifier can be found
      `here <http://rpubs.com/sigven/msi>`__)
   -  Note that the MSI classifier was trained on exome samples.
   -  Will only be present in the report if specified by the user
      ('--msi\_predict')

5. **Mutational signatures**

   -  Estimation of relative contribution of `30 known mutational
      signatures <http://cancer.sanger.ac.uk/cosmic/signatures>`__ in
      tumor sample (using
      `deconstructSigs <https://github.com/raerose01/deconstructSigs>`__
      as the underlying framework)
   -  Datatable with signatures and proposed underlying etiologies
   -  Will only be present in the report if specified by the user
      ('--msig\_identify')
   -  `Trimer (i.e. DNA 3-mer)
      normalization <https://github.com/raerose01/deconstructSigs>`__
      can be configured according to sequencing approach used (WES, WXS
      etc.) using the '--msig\_normalization' option

6. **References**

   -  Supporting scientific literature (key report elements)

-  `View an example report for a breast tumor sample
   (TCGA) <http://folk.uio.no/sigven/tumor_sample.BRCA.0.4.2.pcgr.html>`__
-  `View an example report for a colon adenocarcinoma sample
   (TCGA) <http://folk.uio.no/sigven/tumor_sample.COAD.0.4.2.pcgr.html>`__
-  `View an example report for a lung adenocarcinoma sample
   (TCGA) <http://folk.uio.no/sigven/tumor_sample.LUAD.0.4.2.pcgr.html>`__

The HTML reports have been tested using the following browsers:

-  Safari (10.0.3)
-  Mozilla Firefox (52.0.2)
-  Google Chrome (57.0.2987.110)

Output - Somatic SNVs/InDels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Variant call format - VCF
^^^^^^^^^^^^^^^^^^^^^^^^^

A VCF file containing annotated, somatic calls (single nucleotide
variants and insertion/deletions) is generated with the following naming
convention:

**sample\_id**.pcgr.vcf.gz

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
   Allele\|Consequence\|IMPACT\|SYMBOL\|Gene\|Feature\_type\|Feature\|BIOTYPE\|
   EXON\|INTRON\|HGVSc\|HGVSp\|cDNA\_position\|CDS\_position\|Protein\_position\|Amino\_acids\|
   Codons\|Existing\_variation\|ALLELE\_NUM\|DISTANCE\|STRAND\|FLAGS\|PICK\|VARIANT\_CLASS\|
   SYMBOL\_SOURCE\|HGNC\_ID\|CANONICAL\|APPRIS\|CCDS\|ENSP\|SWISSPROT\|TREMBL\|
   UNIPARC\|RefSeq\|DOMAINS\|HGVS\_OFFSET\|CLIN\_SIG\|SOMATIC\|PHENO\|MOTIF\_NAME\|
   MOTIF\_POS\|HIGH\_INF\_POS\|MOTIF\_SCORE\_CHANGE
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
-  CANCER\_CENSUS\_SOMATIC - Gene with known cancer association -
   `Cancer Gene Census,
   WTSI <http://cancer.sanger.ac.uk/cancergenome/projects/census/>`__
-  CANCER\_CENSUS\_GERMLINE - Gene with known cancer association -
   `Cancer Gene Census,
   WTSI <http://cancer.sanger.ac.uk/cancergenome/projects/census/>`__
-  TUMOR\_SUPPRESSOR - Gene is predicted as tumor suppressor candidate
   according to (`TSGene
   v2.0 <http://bioinfo.mc.vanderbilt.edu/TSGene/>`__)
-  ONCOGENE - Gene is curated as an oncogene according to (`TSGene
   v2.0 <http://bioinfo.mc.vanderbilt.edu/TSGene/>`__)
-  ONCOSCORE - Literature-derived score for cancer gene relevance
   `Bioconductor/OncoScore <http://bioconductor.org/packages/release/bioc/html/OncoScore.html>`__,
   range from 0 (low oncogenic potential) to 1 (high oncogenic
   potential)
-  INTOGEN\_DRIVER - Gene is predicted as a cancer driver in the
   `IntoGen Cancer Drivers Database -
   2014.12 <https://www.intogen.org/downloads>`__

*Variant effect and protein-coding information*
'''''''''''''''''''''''''''''''''''''''''''''''

-  CANCER\_MUTATION\_HOTSPOT - mutation hotspot codon in
   `cancerhotspots.org <http://cancerhotspots.org/>`__. Format:
   gene\_symbol \| codon \| q-value
-  UNIPROT\_FEATURE - Overlapping protein annotations from `UniProt
   KB <http://www.uniprot.org>`__
-  INTOGEN\_DRIVER\_MUT - Indicates if existing variant is predicted as
   driver mutation from IntoGen Catalog of Driver Mutations
-  EFFECT\_PREDICTIONS - Predictions of effect of variant on protein
   function and pre-mRNA splicing from `database of non-synonymous
   functional predictions - dbNSFP
   v3.4 <https://sites.google.com/site/jpopgen/dbNSFP>`__. Predicted
   effects are provided by different sources/algorithms (separated by
   '&'):

   1.  `SIFT <http://provean.jcvi.org/index.php>`__ (Jan 2015)
   2.  `PolyPhen2-HDIV <http://genetics.bwh.harvard.edu/pph2/>`__ (v
       2.2.2)
   3.  `PolyPhen2-HVAR <http://genetics.bwh.harvard.edu/pph2/>`__ (v
       2.2.2)
   4.  `LRT <http://www.genetics.wustl.edu/jflab/lrt_query.html>`__
       (2009)
   5.  `MutationTaster <http://www.mutationtaster.org/>`__ (data release
       Nov 2015)
   6.  `MutationAssessor <http://mutationassessor.org/>`__ (release 3)
   7.  `FATHMM <http://fathmm.biocompute.org.uk>`__ (v2.3)
   8.  `PROVEAN <http://provean.jcvi.org/index.php>`__ (v1.1 Jan 2015)
   9.  `FATHMM\_MKL <http://fathmm.biocompute.org.uk/fathmmMKL.htm>`__
   10. `CADD <http://cadd.gs.washington.edu/>`__ (v1.3)
   11. `DBNSFP\_CONSENSUS\_SVM <https://www.ncbi.nlm.nih.gov/pubmed/25552646>`__
       (Ensembl/consensus prediction, based on support vector machines)
   12. `DBNSFP\_CONSENSUS\_LR <https://www.ncbi.nlm.nih.gov/pubmed/25552646>`__
       (Ensembl/consensus prediction, logistic regression based)
   13. `SPLICE\_SITE\_EFFECT\_ADA <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       adaptive boosting)
   14. `SPLICE\_SITE\_EFFECT\_RF <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       random forest)
   15. `M-CAP <http://bejerano.stanford.edu/MCAP>`__
   16. `MutPred <http://mutpred.mutdb.org>`__
   17. `GERP <http://mendel.stanford.edu/SidowLab/downloads/gerp/>`__

*Variant frequencies/annotations in germline/somatic databases*
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

-  AFR\_AF\_EXAC - African/American germline allele frequency (`Exome
   Aggregation Consortium release
   1 <http://exac.broadinstitute.org/>`__)
-  AMR\_AF\_EXAC - American germline allele frequency (`Exome
   Aggregation Consortium release
   1 <http://exac.broadinstitute.org/>`__)
-  GLOBAL\_AF\_EXAC - Adjusted global germline allele frequency (`Exome
   Aggregation Consortium release
   1 <http://exac.broadinstitute.org/>`__)
-  EAS\_AF\_EXAC - East Asian germline allele frequency (`Exome
   Aggregation Consortium release
   1 <http://exac.broadinstitute.org/>`__)
-  FIN\_AF\_EXAC - Finnish germline allele frequency (`Exome Aggregation
   Consortium release 1 <http://exac.broadinstitute.org/>`__)
-  NFE\_AF\_EXAC - Non-Finnish European germline allele frequency
   (`Exome Aggregation Consortium release
   1 <http://exac.broadinstitute.org/>`__)
-  OTH\_AF\_EXAC - Other germline allele frequency (`Exome Aggregation
   Consortium release 1 <http://exac.broadinstitute.org/>`__)
-  SAS\_AF\_EXAC - South Asian germline allele frequency (`Exome
   Aggregation Consortium release
   1 <http://exac.broadinstitute.org/>`__)
-  AFR\_AF\_GNOMAD - African/American germline allele frequency (`Genome
   Aggregation Database release
   1 <http://gnomad.broadinstitute.org/>`__)
-  AMR\_AF\_GNOMAD - American germline allele frequency (`Genome
   Aggregation Database release
   1 <http://gnomad.broadinstitute.org/>`__)
-  GLOBAL\_AF\_GNOMAD - Adjusted global germline allele frequency
   (`Genome Aggregation Database release
   1 <http://gnomad.broadinstitute.org/>`__)
-  SAS\_AF\_GNOMAD - South Asian germline allele frequency (`Genome
   Aggregation Database release
   1 <http://gnomad.broadinstitute.org/>`__)
-  EAS\_AF\_GNOMAD - East Asian germline allele frequency (`Genome
   Aggregation Database release
   11 <http://gnomad.broadinstitute.org/>`__)
-  FIN\_AF\_GNOMAD - Finnish germline allele frequency (`Genome
   Aggregation Database release
   1 <http://gnomad.broadinstitute.org/>`__)
-  NFE\_AF\_GNOMAD - Non-Finnish European germline allele frequency
   (`Genome Aggregation Database release
   1 <http://gnomad.broadinstitute.org/>`__)
-  OTH\_AF\_GNOMAD - Other germline allele frequency (`Genome
   Aggregation Database release
   1 <http://gnomad.broadinstitute.org/>`__)
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
-  DBSNPRSID - `dbSNP <http://www.ncbi.nlm.nih.gov/SNP/>`__ reference ID
-  DBSNPBUILDID - Initial `dbSNP <http://www.ncbi.nlm.nih.gov/SNP/>`__
   build ID for rsID
-  DBSNP\_MAPPINGSTATUS - Status with respect to the genomic mappability
   of the flanking sequence of the rsID
-  DBSNP\_VALIDATION - Categories of evidence that support the variant
   in `dbSNP <http://www.ncbi.nlm.nih.gov/SNP/>`__
-  DBSNP\_SUBMISSIONS - Number of individual submissions to rsID
-  GWAS\_CATALOG\_PMID - Variant is linked to phenotype through the
   `GWAS Catalog <https://www.ebi.ac.uk/gwas/>`__, literature in PMID
   list
-  GWAS\_CATALOG\_TRAIT\_URI - List of trait URIs for GWAS-associated
   variant
-  COSMIC\_MUTATION\_ID - Mutation identifier in `Catalog of somatic
   mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__
   database
-  COSMIC\_CODON\_FRAC\_GW - For different tumor types, number of
   samples mutated at associated codon position (format:
   codon\_number:tumor\_type:fraction\_mutated). Samples subject to
   exome/genome-wide screens only `Catalog of somatic mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_CODON\_COUNT\_GW - For different tumor types, number of
   samples mutated at associated codon position (format:
   codon\_number:tumor\_type:frequency). Samples subject to
   exome/genome-wide screens only `Catalog of somatic mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__
-  COSMIC\_COUNT\_GW - Global frequency of variant in `Catalog of
   somatic mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_SITE\_HISTOLOGY - Primary site/histology distribution across
   tumor types in `Catalog of somatic mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_CANCER\_TYPE\_GW - Frequency of variant across different
   tumor types in `Catalog of somatic mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__
   - samples subject to exome/genome-wide screens only
-  COSMIC\_CANCER\_TYPE\_ALL - Frequency of variant across different
   tumor types in `Catalog of somatic mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__
-  COSMIC\_SAMPLE\_SOURCE - Sample source distribution for variant in
   `Catalog of somatic mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_DRUG\_RESISTANCE - Targeted drugs/therapies subject to
   resistance in tumors that carry the mutation. `Catalog of somatic
   mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_FATHMM\_PRED - Variant effect prediction from COSMIC's FATHMM
   algorithm (COSMIC variants only) `Catalog of somatic mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_VARTYPE - COSMIC variant type `Catalog of somatic mutations
   in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_CONSEQUENCE - COSMIC consequence type `Catalog of somatic
   mutations in
   cancer <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  ICGC\_PROJECTS - Variant frequency count in different `ICGC Project
   IDs <https://dcc.icgc.org/repository/current/Projects>`__

*Clinical associations*
'''''''''''''''''''''''

-  CLINVAR\_MSID - `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
   Measure Set/Variant ID
-  CLINVAR\_PMIDS - Associated Pubmed IDs for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
-  CLINVAR\_SIG - Clinical significance for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
-  CLINVAR\_VARIANT\_ORIGIN - Origin of variant (somatic, germline, de
   novo etc.) for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
-  DOCM\_DISEASE - Associated disease types for variant in `Database of
   Curated Mutations <http://docm.genome.wustl.edu>`__
-  DOCM\_PMID - Associated Pubmed IDs for variant in `Database of
   Curated Mutations <http://docm.genome.wustl.edu>`__

*Other*
'''''''

-  ANTINEOPLASTIC\_DRUG\_INTERACTION - Approved and experimental
   antineoplastic drugs interacting with the mutated gene, as retrieved
   from the `Drug-Gene Interaction
   Database <http://dgidb.genome.wustl.edu/>`__
-  CIVIC\_ID, CIVIC\_ID\_2 - Variant identifiers in the `CIViC
   database <http://civic.genome.wustl.edu>`__
-  CBMDB\_ID - Variant identifier in the `Cancer bioMarkers
   database <https://www.cancergenomeinterpreter.org/biomarkers>`__

Tab-separated values (TSV)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotated List of all SNVs/InDels
'''''''''''''''''''''''''''''''''

We provide a tab-separated values file with most important annotations
for SNVs/InDels. The file has the following naming convention:

**sample\_id**.pcgr.snvs\_indels.tiers.tsv

The SNVs/InDels are organized into different **tiers** (as defined above
for the HTML report)

The following variables are included in the tiered TSV file:

::

    1. GENOMIC_CHANGE - Identifier for genomic variant, e.g. g.chr1:152382569:A>G
    2. GENOME_VERSION - Assembly version, e.g. GRCh37
    3. VCF_SAMPLE_ID - Sample identifier
    4. VARIANT_CLASS - Variant type, e.g. SNV/insertion/deletion
    5. SYMBOL - Gene symbol
    6. GENE_NAME - Gene description
    7. CCDS - CCDS identifier
    8. ENTREZ_ID - Entrez gene identifier
    9. UNIPROT_ID - UniProt protein identifier
    10. ONCOSCORE - Literature-derived score for cancer gene relevance
    11. ONCOGENE - Gene is curated as an oncogene according to TSGene
    12. TUMOR_SUPPRESSOR - Gene is predicted as tumor suppressor
        candidate according to TSGene
    13. INTOGEN_DRIVER - Gene is predicted as a cancer driver in the
        IntoGen Cancer Drivers Database - 2014.12
    14. CANCER_CENSUS_SOMATIC - Gene with known cancer association -
        Cancer Gene Census, WTSI
    15. CANCER_CENSUS_GERMLINE - Gene with known cancer association -
        Cancer Gene Census, WTSI
    16. CONSEQUENCE - Variant consequence (as defined above for VCF output:
        Consequence)
    17. PROTEIN_CHANGE - Protein change (as defined above for VCF output:
        HGVSp_short)
    18. PROTEIN_DOMAIN - Protein domain
    19. CDS_CHANGE - composite variable for coding change, format:
        Consequence:Feature:cDNA_position:EXON:HGVSp_short
    20. EFFECT_PREDICTIONS - as defined above for VCF
    21. CANCER_MUTATION_HOTSPOT - mutation hotspot codon in
        cancerhotspots.org. Format: gene_symbol | codon | q-value
    22. INTOGEN_DRIVER_MUT - Indicates if existing variant is predicted as
        driver mutation from IntoGen Catalog of Driver Mutations
    23. VEP_ALL_CONSEQUENCE - all VEP consequences
    24. DBSNP - dbSNP reference cluster ID
    25. COSMIC - COSMIC mutation ID
    26. COSMIC_SITE_HISTOLOGY - distribution of tumor sites/histology types
        for COSMIC mutation
    27. COSMIC_DRUG_RESISTANCE - variant associated with resistance to a
        particular antineoplastic drug
    28. CLINVAR - variant origin and associated traits associated with variant
    29. CLINVAR_SIG - clinical significance of CLINVAR variant
    30. GLOBAL_AF_EXAC - adjusted global germline allele frequency in ExAC
    31. GLOBAL_AF_1KG - 1000G Project - phase 3, germline allele frequency
        for all 1000G project samples (global)
    32. CALL_CONFIDENCE - confidence indicator for somatic variant
    33. DP_TUMOR - sequencing depth at variant site (tumor)
    34. AF_TUMOR - allelic fraction of alternate allele (tumor)
    35. DP_NORMAL - sequencing depth at variant site (normal)
    36. AF_NORMAL - allelic fraction of alternate allele (normal)
    37. TIER
    38. TIER_DESCRIPTION

Biomarkers among SNVs/InDEls
''''''''''''''''''''''''''''

For tumor samples that have variant hits in **Tier 1** we provide an
additional file with all associated `clinical evidence
items <https://civic.genome.wustl.edu/#/help/evidence/overview>`__. The
file has the following naming convention:

**sample\_id**.pcgr.snvs\_indels.biomarkers.tsv

The format of the biomarker TSV file is as follows:

::

    1. GENOMIC_CHANGE - Identifier for genomic variant, e.g. g.chr1:152382569:A>G
    2. GENOME_VERSION - Assembly version, e.g. GRCh37
    3. VCF_SAMPLE_ID - Sample identifier
    4. SYMBOL - Gene symbol
    5. CONSEQUENCE - Variant consequence
    6. BM_CLINICAL_SIGNIFICANCE - The association with diagnostic/prognostic end point or treatment
    7. BM_EVIDENCE_LEVEL - The type of experiment from which the evidence is curated (validated, clinical, pre-clinical, case study, and inferential)
    8. BM_EVIDENCE_TYPE - Category of clinical action/relevance implicated by event (Predictive, Prognostic, Predisposing and Diagnostic)
    9. BM_EVIDENCE_DIRECTION - An indicator of whether the evidence statement supports or refutes the clinical significance of an event
    10. BM_CANCER_TYPE - Specific disease or disease subtype that is associated with this event and its clinical implication
    11. BM_THERAPEUTIC_CONTEXT - For predictive evidence, indicates the therapy for which sensitivity or resistance is indicated
    12. BM_RATING - A rating on a 5-star scale, portraying the curators trust in the experiments from which the evidence is curated
    13. BM_CITATION - Publication(s) where the event was described/explored/guidelines/trials
    14. TIER
    15. TIER_DESCRIPTION

Mutational signatures
'''''''''''''''''''''

For each tumor sample, we apply the `deconstructSigs
package <https://github.com/raerose01/deconstructSigs>`__ to delineate
the known mutational signatures. The inferred, weighted contributions by
each signature and their underlying, proposed etiologies are given in a
TSV file with the following naming convention:

**sample\_id**.pcgr.mutational\_signatures.tsv

The format of the mutational signatures TSV file is as follows:

::

    1. Signature_ID - ID of signature from COSMIC's 30 reference signatures
    2. Weight - inferred weight of signature in the tumor sample
    3. Cancer_types - cancer types in which the signature has been observed
    4. Proposed_aetiology - proposed underlying etiology
    5. Trimer_normalization_method - method used for trimer count normalization (deconstructSigs)
    6. SampleID - Sample identifier

Output - Somatic copy number aberrations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Tab-separated values (TSV)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copy number segments are intersected with the genomic coordinates of all
transcripts from (`ENSEMBL/GENCODE's basic gene
annotation <https://www.gencodegenes.org/releases/25lift37.html>`__). In
addition, we attach cancer-relevant annotations for the affected
transcripts. The naming convention of the compressed TSV file is as
follows:

**sample\_id**.pcgr.cna\_segments.tsv.gz

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
    15. gene_biotype - type of gene
    16. cancer_census_germline - gene implicated with germline predisposition to various cancer subtypes
    17. cancer_census_somatic - gene for which somatic mutations have been causally implicated in tumor development
    18. tsgene - tumor suppressor gene status (TSgene database)
    19. tsgene_oncogene - oncogene status (TSgene database)
    20. intogen_drivers - predicted driver gene status (IntoGen Cancer Drivers Database)
    21. antineoplastic_drugs_dgidb - validated and experimental antineoplastic drugs interacting with gene
    22. gencode_transcript_type -
    23. gencode_tag -
    24. gencode_v19 - transcript is part of GENCODE V19
