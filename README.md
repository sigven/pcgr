## Personal Cancer Genome Reporter (PCGR)- variant interpretation report for precision oncology

### Overview

The Personal Cancer Genome Reporter (PCGR) is a stand-alone software package intended for analysis and clinical interpretation of individual cancer genomes. It interprets both somatic SNVs/InDels and copy number aberrations. The software extends basic gene and variant annotations from the [Ensembl’s Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) with oncology-relevant, up-to-date annotations retrieved flexibly through [vcfanno](https://github.com/brentp/vcfanno), and produces HTML reports that can be navigated by clinical oncologists (Figure 1).

![PCGR overview](PCGR_workflow.png)

### Documentation

[![Documentation Status](https://readthedocs.org/projects/pcgr/badge/?version=latest)](http://pcgr.readthedocs.io/en/latest/?badge=latest)


### Annotation resources included in PCGR (v1.2)

* [VEP v85](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor release 85 (GENCODE v19 as the gene reference dataset)
* [COSMIC v80](http://cancer.sanger.ac.uk/cosmic/) - Catalogue of somatic mutations in cancer (February 2017)
* [dBNSFP v3.4](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (March 2017)
* [ICGC v23](https://dcc.icgc.org/) - Somatic mutations discovered in all ICGC (International Cancer Genomics Consortium) tumor cohorts (Dec 2016)
* [ExAC r1](http://exac.broadinstitute.org/) - Germline variant frequencies exome-wide (February 2017)
* [gnomAD r1](http://gnomad.broadinstitute.org/) - Germline variant frequencies exome-wide (March 2017)
* [dbSNP b147](http://www.ncbi.nlm.nih.gov/SNP/) - Database of short genetic variants (April 2016)
* [1000 Genomes Project - phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) - Germline variant frequencies genome-wide (May 2013)
* [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of clinically related variants (March 2017)
* [DoCM](http://docm.genome.wustl.edu) - Database of curated mutations (v3.2, April 2016)
* [CIViC](http://civic.genome.wustl.edu) - Clinical interpretations of variants in cancer (March 17th 2017)
* [CBMDB](http://www.cancergenomeinterpreter.org/biomarkers) - Cancer Biomarkers database (Feb 8th 2017)
* [IntOGen catalog of driver mutations](https://www.intogen.org/downloads) - (May 2016)
* [Cancer Hotspots](http://cancerhotspots.org) - Resource for statistically significant mutations in cancer (2016)
* [UniProt/SwissProt KnowledgeBase 2016_09](http://www.uniprot.org) - Resource on protein sequence and functional information (March 2017)
* [Pfam v30](http://pfam.xfam.org) - Database of protein families and domains (March 2017)
* [DGIdb](http://dgidb.genome.wustl.edu) - Database of interactions betweeen antineoplastic drugs and human proteins (v2.22, February 2016)
* [TSGene v2.0](http://bioinfo.mc.vanderbilt.edu/TSGene/) - Tumor suppressor/oncogene database (November 2015)

### Getting started

#### STEP 1: Installation of Docker

1. [Install the Docker engine](https://docs.docker.com/engine/installation/) on your preferred platform
2. Test that Docker is running
2. Adjust the computing resources dedicated to the Docker, i.e.:
   - Memory: minimum 5GB
   - CPUs: minimum 4

#### STEP 2: Installation of PCGR

1. Download and unpack the [latest release](https://github.com/sigven/pcgr/releases/tag/v1.2)
2. Download and unpack the data bundle (approx. 17Gb) in the PCGR directory
   * Download [data bundle](https://drive.google.com/open?id=0B8aYD2TJ472mb1dqZlpJM2w4aE0) from Google Drive to `~/pcgr-X.X` (replace _X.X_ with the version number)
   * Decompress and untar the bundle, e.g. `gzip -dc pcgr.databundle.vX.X.GRCh37.tgz | tar xvf -`

    A _data/_ folder within the _pcgr-X.X_ software folder should now have been produced
3. Pull the PCGR Docker image from DockerHub:
   * `docker pull sigven/pcgr:latest` (PCGR annotation engine)

#### STEP 3: Input preprocessing

The PCGR workflow accepts two types of input files:

  * A single-sample VCF file with somatic variants (SNVs/InDels)
  * A copy number segment file

PCGR can be run with either or both of the two input files present.

The following requirements __MUST__ be met by the input VCF for PCGR to work properly:

1. Variants in the raw VCF that contain multiple alternative alleles (e.g. "multiple ALTs") must be split into variants with a single alternative allele. A description on how this can be done with the help of [vt](https://github.com/atks/vt) is described within the [documentation page for vcfanno](http://brentp.github.io/vcfanno/#preprocessing)
2. The contents of the VCF must be sorted correctly (i.e. according to chromosomal order and chromosomal position). This can be obtained by [vcftools](https://vcftools.github.io/perl_module.html#vcf-sort).
   * We recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html)
   * 'chr' must be stripped from the chromosome names

#### STEP 4: Run example

Generate report for the example tumor genome (VCF file with SNVs/InDels only and CNV segment file) present in the _examples/_ folder with the following command (_X.X_ should be replaced with the PCGR version number):

`python run_pcgr.py --input_vcf tumor_sample.COAD.vcf.gz --input_cnv_segments tumor_sample.COAD.cna.tsv ~/pcgr-X.X ~/pcgr-X.X/examples tumor_sample.COAD`

This command will run the Docker-based PCGR workflow and produce the following output files in the _examples_ folder:

  1. __tumor_sample.COAD.pcgr.html__ - An interactive HTML report for clinical interpretation
  2. __tumor_sample.COAD.pcgr.vcf.gz__ - VCF file with rich set of annotations for precision oncology
  2. __tumor_sample.COAD.pcgr.snvs_indels.tiers.tsv__ - Tab-separated values file with variants organized according to tiers of functional relevance
  3. __tumor_sample.COAD.pcgr.mutational_signatures.tsv__ - Tab-separated values file with estimated contributions by known mutational signatures and associated underlying etiologies
  4. __tumor_sample.COAD.pcgr.snvs_indels.biomarkers.tsv__ - Tab-separated values file with clinical evidence items associated with biomarkers for diagnosis, prognosis or drug sensitivity/resistance
