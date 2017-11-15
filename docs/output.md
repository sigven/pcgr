## Input & output

### Input

The PCGR workflow accepts two types of input files:

  * An unannotated, single-sample VCF file (>= v4.2) with called somatic variants (SNVs/InDels)
  * A copy number segment file

  __IMPORTANT NOTE: GRCh37 is the reference genome build currently supported by PCGR__

PCGR can be run with either or both of the two input files present.

#### VCF

* We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html)
* If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose)
* Variants used for reporting should be designated as 'PASS' in the VCF FILTER column

__IMPORTANT NOTE 1__: Considering the VCF output for the [numerous somatic SNV/InDel callers](https://www.biostars.org/p/19104/) that have been developed, we have a experienced a general lack of uniformity and robustness for the representation of somatic variant genotype data (e.g. variant allelic depths (tumor/normal), genotype quality etc.). Variant genotype data found as INFO tags in the input VCF can be specified as optional arguments to the PCGR workflow, which in turn can be used for interactive filtering in the tumor report.

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

The report is structured in six main sections, described in more detail below:

  1. __Annotation sources__
      * Lists underlying tools and annotation sources (versions)
  2. __Somatic SNVs/InDels__
      * _Summary statistics_ - indicate number of SNVs/InDels as well as number of coding/non-coding variants
	   	- Note that coding refers to protein-altering and variants at canonical splice-sites (donor, acceptor)
	 * _Mutational burden (TMB)_ - given a coding target region size specified by the user, an estimate of the mutational burden is provided
      * _Tier statistics_ - indicate number of variant found in each tier (see below)
      * _Global distribution - allelic support_ - distribution (histogram) of variant allelic support for somatic variants (will only be present in the report if specific fields in input VCF is defined and specified by the user)
      * _Global variant browser_ - permits filtering of the whole SNV/Indel dataset by various criteria (call confidence, variant sequencing depth, variant consequence etc.)
      * Variants are organized into five tiers (interactive datatables) according to clinical utility
        - _Tier 1_ - constitutes variants linked to predictive, prognostic, diagnostic, and predisposision biomarkers in the [CIViC database](http://civic.genome.wustl.edu) and the [Cancer Biomarkers Database](https://www.cancergenomeinterpreter.org/biomarkers)
        - _Tier 2_ - includes other coding variants that are found in known cancer mutation hotspots, predicted as cancer driver mutations, or curated as disease-causing
        - _Tier 3_ - includes other coding variants found in oncogenes or tumor suppressor genes
        - _Tier 4_ - includes other coding variants
        - _Tier 5_ - includes non-coding variants
            - will only present if specified by the user ('--list_noncoding')
        - __NOTE__: The tier structure is inspired by recommended variant prioritization by [Dienstmann et al., 2014](https://www.ncbi.nlm.nih.gov/pubmed/24768039). The table below shows the correspondence between the terminology for reportable variants used by Dienstmann et al. and the tiers in PCGR:
        <br><br>

        [Dienstmann et al., 2014](https://www.ncbi.nlm.nih.gov/pubmed/24768039) (Figure 2) | PCGR
         --- | ---
         _Actionable_ | Tier 1
         _Other relevant variants_ | Tier 2 & Tier 3
         _Unknown_* | Tier 4 (Tier 5)

         \*While Dienstmann et al. suggest that the _Unknown_ category should be categorized according to pathways, the PCGR employs an arrangement of variant results according to genes, using a literature-derived score of oncogenic potential (KEGG pathway information is also linked).


  3. __Somatic CNA analysis__
      * _Segments - amplifications and homozygous deletions_
         -  Based on user-defined/default log-ratio thresholds of gains/losses, the whole CNA dataset can be navigated further through filters:
	      * cytoband
		 * type of CNA event - *focal* (less than 25% of chromosome arm affected) or *broad*
		 * log ratio
      * _Proto-oncogenes subject to copy number amplifications_
         - Datatable listing known proto-oncogenes covered by user-defined/default amplifications and potential targeted therapies
      * _Tumor suppressor genes subject to homozygous deletions_
         - Datatable listing known tumor suppressor genes covered by user-defined/default losses and potential targeted therapies
      * _Copy number aberrations as biomarkers for prognosis, diagnosis, predisposition, and drug response_
         - Interactive data table where the user can navigate aberrations acting as biomarkers across therapeutic contexts, tumor types, evidence level etc.
  4. __MSI status__
      * Indicates predicted microsatellite stability from the somatic mutation profile and supporting evidence (details of the underlying MSI statistical classifier can be found [here](http://rpubs.com/sigven/msi))
      * Note that the MSI classifier was trained on exome samples.
      * Will only be present in the report if specified by the user in the configuration file ('msi = true')
  5. __Mutational signatures__
      * Estimation of relative contribution of [30 known mutational signatures](http://cancer.sanger.ac.uk/cosmic/signatures) in tumor sample (using [deconstructSigs](https://github.com/raerose01/deconstructSigs) as the underlying framework)
      * Datatable with signatures and proposed underlying etiologies
      * Will only be present in the report if specified by the user in the configuration file ('mutsignatures = true')
      * [Trimer (i.e. DNA 3-mer) normalization](https://github.com/raerose01/deconstructSigs) can be configured according to sequencing approach used (WES, WXS etc.) using the 'mutsignatures_normalization' option, as can the minimum number of mutations required for analysis (option 'mutsignatures_mutation_limit') and the maximum number of mutational signatures in the search space (option 'mutsignatures_signature_limit')
  6. __References__
      * Supporting scientific literature (key report elements)


* [View an example report for a breast tumor sample (TCGA)](http://folk.uio.no/sigven/tumor_sample.BRCA.0.5.0.pcgr.html)
* [View an example report for a colon adenocarcinoma sample (TCGA)](http://folk.uio.no/sigven/tumor_sample.COAD.0.5.0.pcgr.html)


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
  - CSQ - Complete consequence annotations from VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|
  EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|
  Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|
  FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|APPRIS|CCDS|
  ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|DOMAINS|HGVS_OFFSET|HGVSg|AF|
  AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF|gnomAD_AFR_AF|
  gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|
  gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|MOTIF_NAME|
  MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE
  - Consequence - Impact modifier for the consequence type (picked by VEP's --flag\_pick\_allele option)
  - Gene - Ensembl stable ID of affected gene (picked by VEP's --flag\_pick\_allele option)
  - Feature_type - Type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature (picked by VEP's --flag\_pick\_allele option)
  - Feature - Ensembl stable ID of feature (picked by VEP's --flag\_pick\_allele option)
  - cDNA_position - Relative position of base pair in cDNA sequence (picked by VEP's --flag\_pick\_allele option)
  - CDS_position - Relative position of base pair in coding sequence (picked by VEP's --flag\_pick\_allele option)
  - CDS\_CHANGE - Coding, transcript-specific sequence annotation (picked by VEP's --flag\_pick\_allele option)
  - Amino_acid_start - Protein position indicating absolute start of amino acid altered (fetched from Protein_position)
  - Amino_acid_end -  Protein position indicating absolute end of amino acid altered (fetched from Protein_position)
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
  - DISGENET_CUI - Tumor types associated with gene, as found in DisGeNET. Tumor types are listed as [MedGen](https://www.ncbi.nlm.nih.gov/medgen/) concept IDs
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
    2. [LRT](http://www.genetics.wustl.edu/jflab/lrt_query.html) (2009)
    3. [MutationTaster](http://www.mutationtaster.org/) (data release Nov 2015)
    4. [MutationAssessor](http://mutationassessor.org/) (release 3)
    5. [FATHMM](http://fathmm.biocompute.org.uk) (v2.3)
    6. [PROVEAN](http://provean.jcvi.org/index.php) (v1.1 Jan 2015)
    7. [FATHMM_MKL](http://fathmm.biocompute.org.uk/fathmmMKL.htm)
    8. [DBNSFP\_CONSENSUS\_SVM](https://www.ncbi.nlm.nih.gov/pubmed/25552646) (Ensembl/consensus prediction, based on support vector machines)
    9. [DBNSFP\_CONSENSUS\_LR](https://www.ncbi.nlm.nih.gov/pubmed/25552646) (Ensembl/consensus prediction, logistic regression based)
    10. [SPLICE\_SITE\_EFFECT_ADA](http://nar.oxfordjournals.org/content/42/22/13534) (Ensembl/consensus prediction of splice-altering SNVs, based on adaptive boosting)
    11. [SPLICE\_SITE\_EFFECT_RF](http://nar.oxfordjournals.org/content/42/22/13534) (Ensembl/consensus prediction of splice-altering SNVs, based on random forest)
    12. [M-CAP](http://bejerano.stanford.edu/MCAP)
    13. [MutPred](http://mutpred.mutdb.org)
    14. [GERP](http://mendel.stanford.edu/SidowLab/downloads/gerp/)


##### _Variant frequencies/annotations in germline/somatic databases_
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
  - COSMIC\_MUTATION\_ID - Mutation identifier in [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/) database, as provided by VEP
  - *COSMIC\_CODON_FRAC\_GW - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC\_CODON_COUNT\_GW - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC\_SITE\_HISTOLOGY - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC\_CANCER\_TYPE\_GW - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC\_CANCER\_TYPE\_ALL - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC\_SITE\_HISTOLOGY - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC\_CANCER_TYPE\_GW - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC_SAMPLE_SOURCE - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC\_DRUG\_RESISTANCE - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC_FATHMM_PRED - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC_VARTYPE - Deprecated in 0.5.1 due to COSMIC licensing restrictions*
  - *COSMIC_CONSEQUENCE - Deprecated in 0.5.1 due to COSMIC licensing restrictions*


<!--- COSMIC\_MUTATION_ID - Mutation identifier in [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/) database
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
  - COSMIC_CONSEQUENCE - COSMIC consequence type [Catalog of somatic mutations in cancer](http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/)
  - ICGC_PROJECTS - Variant frequency count in different [ICGC Project IDs](https://dcc.icgc.org/repository/current/Projects)
--->

##### _Clinical associations_
  - CLINVAR_MSID - [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) Measure Set/Variant ID
  - CLINVAR_PMIDS - Associated Pubmed IDs for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar)
  - CLINVAR_SIG - Clinical significance for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar)
  - CLINVAR\_VARIANT\_ORIGIN - Origin of variant (somatic, germline, de novo etc.) for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar)
  - DOCM_DISEASE - Associated disease types for variant in [Database of Curated Mutations](http://docm.genome.wustl.edu)
  - DOCM_PMID - Associated Pubmed IDs for variant in [Database of Curated Mutations](http://docm.genome.wustl.edu)

##### _Other_
  - CHEMBL_COMPOUND_ID - antineoplastic drugs targeting the encoded protein (from [Drug-Gene Interaction Database](http://dgidb.genome.wustl.edu/), drugs are listed as [ChEMBL](https://www.ebi.ac.uk/chembl/) compound identifiers)
  - CIVIC\_ID, CIVIC\_ID_2 - Variant identifiers in the [CIViC database](http://civic.genome.wustl.edu), CIVIC_ID refers to markers mapped at variant level, CIVIC_ID_2 refers to region markers (codon, exon etc.)
  - CBMDB_ID - Variant identifier in the [Cancer Biomarkers database](https://www.cancergenomeinterpreter.org/biomarkers)

#### Tab-separated values (TSV)

##### Annotated List of all SNVs/InDels
We provide a tab-separated values file with most important annotations for SNVs/InDels. The file has the following naming convention:

__sample_id__.pcgr.snvs\_indels.tiers.tsv

The SNVs/InDels are organized into different __tiers__ (as defined above for the HTML report)

The following variables are included in the tiered TSV file:

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
    10. ONCOSCORE - Literature-derived score for cancer gene relevance
    11. ONCOGENE - Gene is curated as an oncogene according to TSGene
    12. TUMOR_SUPPRESSOR - Gene is predicted as tumor suppressor
        candidate according to TSGene
    13. DISGENET_CUI - Associated tumor types from DisGeNET (MedGen concept IDs)
    14. DISGENET_TERMS - Associated tumor types from DisGeNET (MedGen concept terms)
    15. CONSEQUENCE - Variant consequence (as defined above for VCF output:
        Consequence)
    16. PROTEIN_CHANGE - Protein change (HGVSp without reference accession)
    17. PROTEIN_DOMAIN - Protein domain
    18. CDS_CHANGE - composite VEP-based variable for coding change, format:
        Consequence:Feature:cDNA_position:EXON:HGVSp_short
    19. HGVSp
    20. HGVSc
    21. EFFECT_PREDICTIONS - as defined above for VCF
    22. CANCER_MUTATION_HOTSPOT - mutation hotspot codon in
        cancerhotspots.org. Format: gene_symbol | codon | q-value
    23. INTOGEN_DRIVER_MUT - Indicates if existing variant is predicted as
        driver mutation from IntoGen Catalog of Driver Mutations
    24. VEP_ALL_CONSEQUENCE - all VEP consequences
    25. DBSNP - dbSNP reference cluster ID
    26. DBSNPBUILDID - initial dbSNP build ID for rsID
    27. COSMIC_MUTATION_ID - COSMIC mutation ID
    28. TCGA_PANCANCER_COUNT - Global count of variant across TCGA cohorts
    29. CHEMBL_COMPOUND_ID - Compounds (as ChEMBL IDs) that target the encoded protein (from DGIdb)
    30. CHEMBL_COMPOUND_TERMS - Compounds (as drug names) that target the encoded protein (from DGIdb)
    31. CLINVAR - variant origin and associated traits associated with variant
    32. CLINVAR_SIG - clinical significance of CLINVAR variant
    33. GLOBAL_AF_GNOMAD - global germline allele frequency in gnomAD
    34. GLOBAL_AF_1KG - 1000G Project - phase 3, germline allele frequency
        for all 1000G project samples (global)
    35. CALL_CONFIDENCE - confidence indicator for somatic variant
    36. DP_TUMOR - sequencing depth at variant site (tumor)
    37. AF_TUMOR - allelic fraction of alternate allele (tumor)
    38. DP_NORMAL - sequencing depth at variant site (normal)
    39. AF_NORMAL - allelic fraction of alternate allele (normal)
    40. TIER
    41. TIER_DESCRIPTION

##### Biomarkers among SNVs/InDEls

For tumor samples that have variant hits in __Tier 1__ we provide an additional file with all associated [clinical evidence items](https://civic.genome.wustl.edu/#/help/evidence/overview). The file has the following naming convention:

__sample_id__.pcgr.snvs\_indels.biomarkers.tsv

The format of the biomarker TSV file is as follows:

    1. GENOMIC_CHANGE - Identifier for variant at the genome (VCF) level, e.g. 1:g.152382569A>G
    	  Format: (<chrom>:g.<position><ref_allele>><alt_allele>)
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

##### Mutational signatures

For each tumor sample, we apply the [deconstructSigs package](https://github.com/raerose01/deconstructSigs) to delineate the known mutational signatures. The inferred, weighted contributions by each signature and their underlying, proposed etiologies are given in a TSV file with the following naming convention:

__sample_id__.pcgr.mutational\_signatures.tsv

The format of the mutational signatures TSV file is as follows:

  	1. Signature_ID - ID of signature from COSMIC's 30 reference signatures
  	2. Weight - inferred weight of signature in the tumor sample
  	3. Cancer_types - cancer types in which the signature has been observed
  	4. Proposed_aetiology - proposed underlying etiology
     5. Trimer_normalization_method - method used for trimer count normalization (deconstructSigs)
  	6. SampleID - Sample identifier


### Output - Somatic copy number aberrations

#### 1. Tab-separated values (TSV)

 Copy number segments are intersected with the genomic coordinates of all transcripts from ([ENSEMBL/GENCODE's basic gene annotation](https://www.gencodegenes.org/releases/27lift37.html)). In addition, we attach cancer-relevant annotations for the affected transcripts. The naming convention of the compressed TSV file is as follows:

__sample_id__.pcgr.cna_segments.tsv.gz

The format of the compressed TSV file is the following:

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
