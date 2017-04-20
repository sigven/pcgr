## Input & output

### Input

The PCGR workflow accepts two types of input files:

  * An unannotated, single-sample VCF file (>= v4.2) with called somatic variants (SNVs/InDels)
  * A copy number segment file

PCGR can be run with either or both of the two input files present.

__IMPORTANT NOTE__: Only the GRCh37 version of the human genome is currently supported.

#### VCF

* We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html)
* If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose)

__IMPORTANT NOTE 1__: Considering the VCF output for the [numerous somatic SNV/InDel callers](https://www.biostars.org/p/19104/) that have been developed, we have a experienced a general lack of uniformity and robustness for the representation of somatic variant genotype data (e.g. variant allelic depths (tumor/normal), genotype quality etc.). In the output results provided within the current version of PCGR, we are considering PASSed variants only, and variant genotype data (i.e. as found in the VCF SAMPLE columns) are not handled or parsed. As improved standards for this matter may emerge, we will strive to include this information in the annotated output files.

__IMPORTANT NOTE 2__: PCGR generates a number of VCF INFO annotation tags that is appended to the query VCF. We will therefore encourage the users to submit query VCF files that have not been subject to annotations by other means, but rather a VCF file that comes directly from variant calling. If not, there are likely to be INFO tags in the query VCF file that coincide with those produced by PCGR.

#### Copy number segments

The tab-separated values file with copy number aberrations __MUST__ contain the following four columns:

* Chromosome
* Start
* End
* Segment_Mean

Here, _Chromosome_, _Start_, and _End_ denote the chromosomal segment (GRCh37), and __Segment_Mean__ denotes the log(2) ratio for a particular segment, which is a common output of somatic copy number alteration callers. Below shows the initial part of a copy number segment file that is formatted correctly according to PCGR's requirements:

      Chromosome	Start	End	Segment_Mean
      1 3218329 3550598 0.0024
      1 3552451 4593614 0.1995
      1 4593663 6433129 -1.0277


### Output - Interactive HTML report

An interactive and tier-structured HTML report that shows the most relevant findings in the query cancer genome is provided with the following naming convention:

__sample_id__.pcgr.html

The __sample_id__ is provided as input by the user, and reflects a unique identifier of the tumor-normal sample pair to be analyzed.

* [View an example report for a breast tumor sample (TCGA)](http://folk.uio.no/sigven/tumor_sample.BRCA.pcgr.html)
* [View an example report for a colorectal tumor sample (TCGA)](http://folk.uio.no/sigven/tumor_sample.COAD.pcgr.html)

The HTML reports have been tested using the following browsers:

* Safari (10.0.3)
* Mozilla Firefox (52.0.2)
* Google Chrome (57.0.2987.110)

### Output - Somatic SNVs/InDels

#### Variant call format - VCF

A VCF file containing annotated, somatic calls (single nucleotide variants and insertion/deletions) is generated with the following naming convention:

__sample_id__.pcgr.vcf.gz

Here, the __sample_id__ is provided as input by the user, and reflects a unique identifier of the tumor-normal sample pair to be analyzed. Following common standards, the annotated VCF file is compressed with [bgzip](http://www.htslib.org/doc/tabix.html) and indexed with [tabix](http://www.htslib.org/doc/tabix.html). Below follows a description of all annotations/tags present in the VCF INFO column after processing with the PCGR annotation pipeline:

##### _VEP consequence annotations_
  - CSQ - Complete consequence annotations from VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature\_type|Feature|BIOTYPE| EXON|INTRON|HGVSc|HGVSp|cDNA\_position|CDS\_position|Protein\_position|Amino\_acids| Codons|Existing\_variation|ALLELE\_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT\_CLASS| SYMBOL\_SOURCE|HGNC\_ID|CANONICAL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL| UNIPARC|RefSeq|DOMAINS|HGVS\_OFFSET|CLIN\_SIG|SOMATIC|PHENO|MOTIF_NAME| MOTIF\_POS|HIGH\_INF\_POS|MOTIF\_SCORE\_CHANGE
  - Consequence	 - Impact modifier for the consequence type (picked by VEP's --flag\_pick\_allele option)
  - Gene - Ensembl stable ID of affected gene (picked by VEP's --flag\_pick\_allele option)
  - Feature_type - Type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature (picked by VEP's --flag\_pick\_allele option)
  - Feature - Ensembl stable ID of feature (picked by VEP's --flag\_pick\_allele option)
  - cDNA_position - Relative position of base pair in cDNA sequence (picked by VEP's --flag\_pick\_allele option)
  - CDS_position - Relative position of base pair in coding sequence (picked by VEP's --flag\_pick\_allele option)
  - CDS\_CHANGE - Coding, transcript-specific sequence annotation (picked by VEP's --flag\_pick\_allele option)
  - Protein_position - Relative position of amino acid in protein (picked by VEP's --flag\_pick\_allele option)
  - Amino_acids - Only given if the variant affects the protein-coding sequence (picked by VEP's --flag\_pick\_allele option)
  - Codons - The alternative codons with the variant base in upper case (picked by VEP's --flag\_pick\_allele option)
  - IMPACT - Impact modifier for the consequence type (picked by VEP's --flag\_pick\_allele option)
  - VARIANT_CLASS - Sequence Ontology variant class (picked by VEP's --flag\_pick\_allele option)
  - SYMBOL - Gene symbol (picked by VEP's --flag\_pick\_allele option)
  - SYMBOL_SOURCE - The source of the gene symbol (picked by VEP's --flag\_pick\_allele option)
  - STRAND - The DNA strand (1 or -1) on which the transcript/feature lies (picked by VEP's --flag\_pick\_allele option)
  - ENSP - The Ensembl protein identifier of the affected transcript (picked by VEP's --flag\_pick\_allele option)
  - FLAGS - Transcript quality flags: cds\_start\_NF: CDS 5', incomplete cds\_end\_NF: CDS 3' incomplete (picked by VEP's --flag\_pick\_allele option)
  - SWISSPROT - Best match UniProtKB/Swiss-Prot accession of protein product (picked by VEP's --flag\_pick\_allele option)
  - TREMBL - Best match UniProtKB/TrEMBL accession of protein product (picked by VEP's --flag\_pick\_allele option)
  - UNIPARC - Best match UniParc accession of protein product (picked by VEP's --flag\_pick\_allele option)
  - HGVSc - The HGVS coding sequence name (picked by VEP's --flag\_pick\_allele option)
  - HGVSp - The HGVS protein sequence name (picked by VEP's --flag\_pick\_allele option)
  - HGVSp_short - The HGVS protein sequence name, short version (picked by VEP's --flag\_pick\_allele option)
  - HGVS_OFFSET - Indicates by how many bases the HGVS notations for this variant have been shifted (picked by VEP's --flag\_pick\_allele option)
  - MOTIF_NAME - The source and identifier of a transcription factor binding profile aligned at this position (picked by VEP's --flag\_pick\_allele option)
  - MOTIF_POS - The relative position of the variation in the aligned TFBP (picked by VEP's --flag\_pick\_allele option)
  - HIGH\_INF\_POS - A flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP) (picked by VEP's --flag\_pick\_allele option)
  - MOTIF\_SCORE\_CHANGE - The difference in motif score of the reference and variant sequences for the TFBP (picked by VEP's --flag\_pick\_allele option)
  - CELL_TYPE - List of cell types and classifications for regulatory feature (picked by VEP's --flag\_pick\_allele option)
  - CANONICAL - A flag indicating if the transcript is denoted as the canonical transcript for this gene (picked by VEP's --flag\_pick\_allele option)
  - CCDS - The CCDS identifier for this transcript, where applicable (picked by VEP's --flag\_pick\_allele option)
  - INTRON - The intron number (out of total number) (picked by VEP's --flag\_pick\_allele option)
  - EXON - The exon number (out of total number) (picked by VEP's --flag\_pick\_allele option)
  - DOMAINS - The source and identifier of any overlapping protein domains (picked by VEP's --flag\_pick\_allele option)
  - DISTANCE - Shortest distance from variant to transcript (picked by VEP's --flag\_pick\_allele option)
  - BIOTYPE - Biotype of transcript or regulatory feature (picked by VEP's --flag\_pick\_allele option)
  - TSL - Transcript support level (picked by VEP's --flag\_pick\_allele option)>
  - PUBMED - PubMed ID(s) of publications that cite existing variant - VEP
  - PHENO - Indicates if existing variant is associated with a phenotype, disease or trait - VEP
  - GENE_PHENO - Indicates if overlapped gene is associated with a phenotype, disease or trait - VEP
  - ALLELE_NUM - Allele number from input; 0 is reference, 1 is first alternate etc - VEP
  - REFSEQ_MATCH - The RefSeq transcript match status; contains a number of flags indicating whether this RefSeq transcript matches the underlying reference sequence and/or an Ensembl transcript (picked by VEP's --flag\_pick\_allele option)
  - PICK - Indicates if this block of consequence data was picked by VEP's --flag\_pick\_allele option
  - VEP\_ALL\_CONSEQUENCE - All transcript consequences (Consequence:SYMBOL:Feature_type:Feature:BIOTYPE) - VEP

##### _Gene information_
  - ENTREZ_ID - [Entrez](http://www.ncbi.nlm.nih.gov/gene) gene identifier
  - APPRIS - Principal isoform flags according to the [APPRIS principal isoform database](http://appris.bioinfo.cnio.es/#/downloads)
  - UNIPROT_ID - [UniProt](http://www.uniprot.org) identifier
  - CANCER\_CENSUS_SOMATIC - Gene with known cancer association - [Cancer Gene Census, WTSI](http://cancer.sanger.ac.uk/cancergenome/projects/census/)
  - CANCER\_CENSUS_GERMLINE - Gene with known cancer association - [Cancer Gene Census, WTSI](http://cancer.sanger.ac.uk/cancergenome/projects/census/)
  - TUMOR_SUPPRESSOR - Gene is predicted as tumor suppressor candidate according to ([TSGene v2.0](http://bioinfo.mc.vanderbilt.edu/TSGene/))
  - ONCOGENE - Gene is curated as an oncogene according to ([TSGene v2.0](http://bioinfo.mc.vanderbilt.edu/TSGene/))
  - ONCOSCORE - Literature-derived score for cancer gene relevance [Bioconductor/OncoScore](http://bioconductor.org/packages/release/bioc/html/OncoScore.html), range from 0 (low oncogenic potential) to 1 (high oncogenic potential)
  - INTOGEN_DRIVER - Gene is predicted as a cancer driver in the [IntoGen Cancer Drivers Database - 2014.12](https://www.intogen.org/downloads)


##### _Variant effect and protein-coding information_
  - CANCER\_MUTATION\_HOTSPOT - mutation hotspot codon in [cancerhotspots.org](http://cancerhotspots.org/). Format: gene_symbol | codon | q-value
  - UNIPROT\_FEATURE - Overlapping protein annotations from [UniProt KB](http://www.uniprot.org)
  - INTOGEN\_DRIVER\_MUT - Indicates if existing variant is predicted as driver mutation from IntoGen Catalog of Driver Mutations
  - EFFECT\_PREDICTIONS - Predictions of effect of variant on protein function and pre-mRNA splicing from [database of non-synonymous functional predictions - dbNSFP v3.4](https://sites.google.com/site/jpopgen/dbNSFP). Predicted effects are provided by different sources/algorithms (separated by '&'):

    1. [SIFT](http://provean.jcvi.org/index.php) (Jan 2015)
    2. [PolyPhen2-HDIV](http://genetics.bwh.harvard.edu/pph2/) (v 2.2.2)
    3. [PolyPhen2-HVAR](http://genetics.bwh.harvard.edu/pph2/) (v 2.2.2)
    4. [LRT](http://www.genetics.wustl.edu/jflab/lrt_query.html) (2009)
    5. [MutationTaster](http://www.mutationtaster.org/) (data release Nov 2015)
    6. [MutationAssessor](http://mutationassessor.org/) (release 3)
    7. [FATHMM] (http://fathmm.biocompute.org.uk) (v2.3)
    8. [PROVEAN](http://provean.jcvi.org/index.php) (v1.1 Jan 2015)
    9. [FATHMM_MKL](http://fathmm.biocompute.org.uk/fathmmMKL.htm)
    10. [CADD](http://cadd.gs.washington.edu/) (v1.3)
    11. [DBNSFP\_CONSENSUS\_SVM](https://www.ncbi.nlm.nih.gov/pubmed/25552646) (Ensembl/consensus prediction, based on support vector machines)
    12. [DBNSFP\_CONSENSUS\_LR](https://www.ncbi.nlm.nih.gov/pubmed/25552646) (Ensembl/consensus prediction, logistic regression based)
    13. [SPLICE\_SITE\_EFFECT_ADA](http://nar.oxfordjournals.org/content/42/22/13534) (Ensembl/consensus prediction of splice-altering SNVs, based on adaptive boosting)
    14. [SPLICE\_SITE\_EFFECT_RF](http://nar.oxfordjournals.org/content/42/22/13534) (Ensembl/consensus prediction of splice-altering SNVs, based on adaptive boosting)
    15. [M-CAP](http://bejerano.stanford.edu/MCAP)
    16. [REVEL](https://www.ncbi.nlm.nih.gov/pubmed/27666373)
    17. [MutPred](http://mutpred.mutdb.org)
    18. [GERP](http://mendel.stanford.edu/SidowLab/downloads/gerp/)


##### _Variant frequencies/annotations in germline/somatic databases_
  - AFR\_AF\_EXAC - African/American germline allele frequency ([Exome Aggregation Consortium release 1](http://exac.broadinstitute.org/))
  - AMR\_AF\_EXAC - American germline allele frequency ([Exome Aggregation Consortium release 1](http://exac.broadinstitute.org/))
  - GLOBAL\_AF\_EXAC - Adjusted global germline allele frequency ([Exome Aggregation Consortium release 1](http://exac.broadinstitute.org/))
  - EAS\_AF\_EXAC - East Asian germline allele frequency ([Exome Aggregation Consortium release 1](http://exac.broadinstitute.org/))
  - FIN\_AF\_EXAC - Finnish germline allele frequency ([Exome Aggregation Consortium release 1](http://exac.broadinstitute.org/))
  - NFE\_AF\_EXAC - Non-Finnish European germline allele frequency ([Exome Aggregation Consortium release 1](http://exac.broadinstitute.org/))
  - OTH\_AF\_EXAC - Other germline allele frequency ([Exome Aggregation Consortium release 1](http://exac.broadinstitute.org/))
  - SAS\_AF\_EXAC - South Asian germline allele frequency ([Exome Aggregation Consortium release 1](http://exac.broadinstitute.org/))
  - AFR\_AF\_GNOMAD - African/American germline allele frequency ([Genome Aggregation Database release 1](http://gnomad.broadinstitute.org/))
  - AMR\_AF\_GNOMAD - American germline allele frequency ([Genome Aggregation Database release 1](http://gnomad.broadinstitute.org/))
  - GLOBAL\_AF\_GNOMAD - Adjusted global germline allele frequency ([Genome Aggregation Database release 1](http://gnomad.broadinstitute.org/))
  - SAS\_AF\_GNOMAD - South Asian germline allele frequency ([Genome Aggregation Database release 1](http://gnomad.broadinstitute.org/))
  - EAS\_AF\_GNOMAD - East Asian germline allele frequency ([Genome Aggregation Database release 11](http://gnomad.broadinstitute.org/))
  - FIN\_AF\_GNOMAD - Finnish germline allele frequency ([Genome Aggregation Database release 1](http://gnomad.broadinstitute.org/))
  - NFE\_AF\_GNOMAD - Non-Finnish European germline allele frequency ([Genome Aggregation Database release 1](http://gnomad.broadinstitute.org/))
  - OTH\_AF\_GNOMAD - Other germline allele frequency ([Genome Aggregation Database release 1](http://gnomad.broadinstitute.org/))
  - AFR\_AF\_1KG - [1000G Project - phase 3](http://www.1000genomes.org) germline allele frequency for samples from AFR (African)
  - AMR\_AF\_1KG - [1000G Project - phase 3](http://www.1000genomes.org) germline allele frequency for samples from AMR (Ad Mixed American)
  - EAS\_AF\_1KG - [1000G Project - phase 3](http://www.1000genomes.org) germline allele frequency for samples from EAS (East Asian)
  - EUR\_AF\_1KG - [1000G Project - phase 3](http://www.1000genomes.org) germline allele frequency for samples from EUR (European)
  - SAS\_AF\_1KG - [1000G Project - phase 3](http://www.1000genomes.org) germline allele frequency for samples from SAS (South Asian)
  - GLOBAL\_AF\_1KG - [1000G Project - phase 3](http://www.1000genomes.org) germline allele frequency for all 1000G project samples (global)
  - DBSNPRSID - [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) reference ID
  - DBSNPBUILDID - Initial [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) build ID for rsID
  - DBSNP_MAPPINGSTATUS - Status with respect to the genomic mappability of the flanking sequence of the rsID
  - DBSNP_VALIDATION - Categories of evidence that support the variant in [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/)
  - DBSNP_SUBMISSIONS - Number of individual submissions to rsID
  - GWAS\_CATALOG_PMID - Variant is linked to phenotype through the [GWAS Catalog](https://www.ebi.ac.uk/gwas/), literature in PMID list
  - GWAS\_CATALOG\_TRAIT_URI - List of trait URIs for GWAS-associated variant
  - COSMIC\_MUTATION_ID - Mutation identifier in [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/) database
  - COSMIC\_CODON\_FRAC\_GW - For different tumor types, number of samples mutated at associated codon position (format: codon\_number:tumor\_type:fraction_mutated). Samples subject to exome/genome-wide screens only [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/).
  - COSMIC\_CODON\_COUNT\_GW - For different tumor types, number of samples mutated at associated codon position (format: codon\_number:tumor\_type:frequency). Samples subject to exome/genome-wide screens only [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/)
  - COSMIC\_COUNT\_GW - Global frequency of variant in [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/).
  - COSMIC\_SITE\_HISTOLOGY - Primary site/histology distribution across tumor types in [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/).
  - COSMIC\_CANCER\_TYPE\_GW - Frequency of variant across different tumor types in [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/) - samples subject to exome/genome-wide screens only
  - COSMIC\_CANCER\_TYPE\_ALL - Frequency of variant across different tumor types in [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/)
  - COSMIC\_SAMPLE\_SOURCE - Sample source distribution for variant in [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/).
  - COSMIC\_DRUG\_RESISTANCE - Targeted drugs/therapies subject to resistance in tumors that carry the mutation. [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/).
  - COSMIC\_FATHMM\_PRED - Variant effect prediction from COSMIC's FATHMM algorithm (COSMIC variants only) [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/).
  - COSMIC_VARTYPE - COSMIC variant type [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/).
  - COSMIC_CONSEQUENCE - COSMIC consequence type [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/).
  - ICGC_PROJECTS - Variant frequency count in different [ICGC Project IDs](https://dcc.icgc.org/repository/current/Projects)


##### _Clinical associations_
  - CLINVAR_MSID - [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) Measure Set/Variant ID
  - CLINVAR_PMIDS - Associated Pubmed IDs for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar)
  - CLINVAR_SIG - Clinical significance for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar)
  - CLINVAR\_VARIANT\_ORIGIN - Origin of variant (somatic, germline, de novo etc.) for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar)
  - DOCM_DISEASE - Associated disease types for variant in [Database of Curated Mutations](http://docm.genome.wustl.edu)
  - DOCM_PMID - Associated Pubmed IDs for variant in [Database of Curated Mutations](http://docm.genome.wustl.edu)

##### _Other_
  - ANTINEOPLASTIC\_DRUG\_INTERACTION - Approved and experimental antineoplastic drugs interacting with the mutated gene, as retrieved from the [Drug-Gene Interaction Database](http://dgidb.genome.wustl.edu/)
  - CIVIC\_ID, CIVIC\_ID_2 - Variant identifiers in the [CIViC database](http://civic.genome.wustl.edu)
  - CBMDB_ID - Variant identifier in the [Cancer bioMarkers database](https://www.cancergenomeinterpreter.org/biomarkers)

#### Tab-separated values (TSV)

##### Annotated List of all SNVs/InDels
We provide a tab-separated values file with most important annotations for SNVs/InDels. The file has the following naming convention:

__sample_id__.pcgr.snvs\_indels.tiers.tsv

The SNVs/InDels are organized into different __tiers__ that reflect relevance for therapeutics/tumorigenesis:

 - __Tier 1__ constitute variants recorded as prognostic/diagnostic/drug sensitivity biomarkers in the [CIViC database](http://civic.genome.wustl.edu) and the [Cancer Biomarkers Database](https://www.cancergenomeinterpreter.org/biomarkers)
 - __Tier 2__ includes other coding variants that are found in known mutational hotspots, predicted as cancer driver mutations, or curated as disease-causing
 - __Tier 3__ includes other coding variants found in oncogenes, tumor suppressor genes, or cancer census genes
 - __Tier 4__ includes other coding variants
 - __Tier 5__ includes non-coding variants

 __Note__: '_coding variants_' refer to the set of variants with the following consequences:
   - missense variant
   - splice donor/splice acceptor alteration
   - stop gained/stop lost
   - frameshift/non-frameshift variants

The following variables are included in the tiered TSV file:

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
    32. TIER
    33. TIER_DESCRIPTION

##### Biomarkers among SNVs/InDEls

For tumor samples that have variant hits in __Tier 1__ we provide an additional file with all associated [clinical evidence items](https://civic.genome.wustl.edu/#/help/evidence/overview). The file has the following naming convention:

__sample_id__.pcgr.snvs\_indels.biomarkers.tsv

The format of the biomarker TSV file is as follows:

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

##### Mutational signatures

For each tumor sample, we apply the [deconstructSigs package](https://github.com/raerose01/deconstructSigs) to delineate the known mutational signatures. The inferred, weighted contributions by each signature and their underlying, proposed etiologies are given in a TSV file with the following naming convention:

__sample_id__.pcgr.mutational\_signatures.tsv

The format of the mutational signatures TSV file is as follows:

	1. Signature_ID - ID of signature from COSMIC's 30 reference signatures
	2. Weight - inferred weight of signature in the tumor sample
	3. Cancer_types - cancer types in which the signature has been observed
	4. Proposed_aetiology - proposed underlying etiology
	5. SampleID - Sample identifier


### Output - Somatic copy number abberations

#### 1. Tab-separated values (TSV)

 Copy number segments are intersected with the genomic coordinates of all transcripts from ([ENSEMBL/GENCODE's basic gene annotation](https://www.gencodegenes.org/releases/25lift37.html). In adddition, we attach cancer-relevant annotations for the affected transcripts. The naming convention of the compressed TSV file is as follows:

__sample_id__.pcgr.cna_segments.tsv.gz

The format of the compressed TSV file is the following:

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
    16. tsgene - tumor suppressor gene status (TSgene database)
    17. tsgene_oncogene - oncogene status (TSgene database)
    18. intogen_drivers - predicted driver gene status (IntoGen Cancer Drivers Database)
    19. antineoplastic_drugs_dgidb - validated and experimental antineoplastic drugs interacting with gene
    20. gencode_transcript_type -
    21. gencode_tag -
    22. gencode_v19 - transcript is part of GENCODE V19
