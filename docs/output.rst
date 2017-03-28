Input & output
--------------

Input
~~~~~

The PCGR workflow accepts two types of input files:

-  A single-sample VCF file with somatic variants (SNVs/InDels)
-  A copy number segment file

PCGR can be run with either or both of the two input files present.

**IMPORTANT NOTE**: Only the GRCh37 version of the human genome is
currently supported.

VCF preprocessing
^^^^^^^^^^^^^^^^^

The following requirements **MUST** be met by the input VCF for PCGR to
work properly:

1. Variants in the raw VCF that contain multiple alternative alleles
   (e.g. "multiple ALTs") must be split into variants with a single
   alternative allele. A description on how this can be done with the
   help of `vt <https://github.com/atks/vt>`__ is described within the
   `documentation page for
   vcfanno <http://brentp.github.io/vcfanno/#preprocessing>`__
2. The contents of the VCF must be sorted correctly (i.e. according to
   chromosomal order and chromosomal position). This can be obtained by
   `vcftools <https://vcftools.github.io/perl_module.html#vcf-sort>`__.

   -  We recommend that the input VCF is compressed and indexed using
      `bgzip <http://www.htslib.org/doc/tabix.html>`__ and
      `tabix <http://www.htslib.org/doc/tabix.html>`__
   -  'chr' must be stripped from the chromosome names

Output - Interactive HTML report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO

Output - Somatic SNVs/InDels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Variant call format - VCF
^^^^^^^^^^^^^^^^^^^^^^^^^

Each tumor-normal sample pair is provided with a VCF file containing
annotated, somatic calls (single nucleotide variants and
insertion/deletions) that has the following naming convention:

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
   v3.2 <https://sites.google.com/site/jpopgen/dbNSFP>`__. Predicted
   effects are provided by 14 different sources/algorithms (separated by
   '&'):

   1.  `SIFT <http://provean.jcvi.org/index.php>`__ (Jan 2015)
   2.  `PolyPhen2 <http://genetics.bwh.harvard.edu/pph2/>`__ (v 2.2.2,
       predictions based on
       `HumDiv <http://genetics.bwh.harvard.edu/pph2/dokuwiki/overview#prediction>`__)
   3.  `LRT <http://www.genetics.wustl.edu/jflab/lrt_query.html>`__
       (2009)
   4.  `MutationTaster <http://www.mutationtaster.org/>`__ (data release
       Nov 2015)
   5.  `MutationAssessor <http://mutationassessor.org/>`__ (release 3)
   6.  [FATHMM] (http://fathmm.biocompute.org.uk) (v2.3)
   7.  `PROVEAN <http://provean.jcvi.org/index.php>`__ (v1.1 Jan 2015)
   8.  `FATHMM\_MKL <http://fathmm.biocompute.org.uk/fathmmMKL.htm>`__
   9.  `CADD <http://cadd.gs.washington.edu/>`__ (v1.3)
   10. `DBNSFP\_CONSENSUS\_SVM <https://www.ncbi.nlm.nih.gov/pubmed/25552646>`__
       (Ensembl/consensus prediction, based on support vector machines)
   11. `DBNSFP\_CONSENSUS\_LR <https://www.ncbi.nlm.nih.gov/pubmed/25552646>`__
       (Ensembl/consensus prediction, logistic regression based)
   12. `SPLICE\_SITE\_EFFECT\_ADA <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       adaptive boosting)
   13. `SPLICE\_SITE\_EFFECT\_RF <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       adaptive boosting)

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
   mutations in cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__
   database
-  COSMIC\_CODON\_FRAC\_GW - For different tumor types, number of
   samples mutated at associated codon position (format:
   codon\_number:tumor\_type:fraction\_mutated). Samples subject to
   exome/genome-wide screens only `Catalog of somatic mutations in
   cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_CODON\_COUNT\_GW - For different tumor types, number of
   samples mutated at associated codon position (format:
   codon\_number:tumor\_type:frequency). Samples subject to
   exome/genome-wide screens only `Catalog of somatic mutations in
   cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__
-  COSMIC\_COUNT\_GW - Global frequency of variant in `Catalog of
   somatic mutations in cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_SITE\_HISTOLOGY - Primary site/histology distribution across
   tumor types in `Catalog of somatic mutations in cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_CANCER\_TYPE\_GW - Frequency of variant across different
   tumor types in `Catalog of somatic mutations in cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__ -
   samples subject to exome/genome-wide screens only
-  COSMIC\_CANCER\_TYPE\_ALL - Frequency of variant across different
   tumor types in `Catalog of somatic mutations in cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__
-  COSMIC\_SAMPLE\_SOURCE - Sample source distribution for variant in
   `Catalog of somatic mutations in cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_DRUG\_RESISTANCE - Targeted drugs/therapies subject to
   resistance in tumors that carry the mutation. `Catalog of somatic
   mutations in cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_FATHMM\_PRED - Variant effect prediction from COSMIC's FATHMM
   algorithm (COSMIC variants only) `Catalog of somatic mutations in
   cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_VARTYPE - COSMIC variant type `Catalog of somatic mutations
   in cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
-  COSMIC\_CONSEQUENCE - COSMIC consequence type `Catalog of somatic
   mutations in cancer - COSMIC
   v78 <http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/>`__.
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

The SNVs/InDels are organized into different **tiers** that reflect
relevance for therapeutics/tumorigenesis:

-  **Tier 1** constitute variants recorded as prognostic/diagnostic/drug
   sensitivity biomarkers in the `CIViC
   database <http://civic.genome.wustl.edu>`__ and the `Cancer
   Biomarkers
   Database <https://www.cancergenomeinterpreter.org/biomarkers>`__
-  **Tier 2** includes other coding variants that are found in known
   mutational hotspots, predicted as cancer driver mutations, or curated
   as disease-causing
-  **Tier 3** includes other coding variants found in oncogenes, tumor
   suppressor genes, or cancer census genes
-  **Tier 4** includes other coding variants
-  **Tier 5** includes non-coding variants

**Note**: '*coding variants*' refer to the set of variants with the
following consequences: - missense variant - splice donor/splice
acceptor alteration - stop gained/stop lost - frameshift/non-frameshift
variants

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
    30. GLOBAL_AF_EXAC - adjusted global germline allele frequency in ExAC release 0.3.1
    31. GLOBAL_AF_1KG - 1000G Project - phase 3, germline allele frequency
        for all 1000G project samples (global)
    32. DP_TUMOR - tumor depth
    33. AF_TUMOR - allelic fraction tumor
    34. DP_NORMAL - normal depth
    35. AF_NORMAL - allelic fraction normal
    36. TIER
    37. TIER_DESCRIPTION

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
    10. BM_DISEASE_NAME - Specific disease or disease subtype that is associated with this event and its clinical implication
    11. BM_DRUG_NAMES - For predictive evidence, indicates the therapy for which sensitivity or resistance is indicated
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
    5. SampleID - Sample identifier

Output - Somatic copy number abberations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Tab-separated values (TSV)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copy number segments are intersected with the genomic coordinates of all
transcripts from (`ENSEMBL/GENCODE's basic gene
annotation <https://www.gencodegenes.org/releases/25lift37.html>`__. In
adddition, we attach cancer-relevant annotations for the affected
transcripts. The naming convention of the compressed TSV file is as
follows:

**sample\_id**.pcgr.cnv.tsv.gz

The format of the compressed TSV file is the following:

::

    1. chrom - chromosome  
    2. segment_start - start of copy number segment
    3. segment_end - end of copy number segment
    4. segment_length - length of segment in Mb
    5. LogR - Copy log-ratio
    6. ensembl_gene_id
    7. symbol - gene symbol
    8. ensembl_transcript_id
    9. transcript_start
    10. transcript_end
    11. transcript_overlap_percent - percent of transcript length covered by CN segment
    12. name - gene name description
    13. gene_biotype - type of gene
    14. cancer_census_germline - gene implicated with germline predisposition to various cancer subtypes
    15. cancer_census_somatic - gene for which somatic mutations have been causally implicated in tumor development
    16. tsgene - tumor suppressor gene status (TSgene 2.0 database)
    17. tsgene_oncogene - oncogene status (TSgene 2.0 database)
    18. intogen_drivers - predicted driver gene status (IntoGen Cancer Drivers Database)
    19. antineoplastic_drugs_dgidb - validated and experimental antineoplastic drugs interacting with gene
    20. gencode_v19 - transcript is part of GENCODE V19
