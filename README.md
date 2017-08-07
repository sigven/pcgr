## Personal Cancer Genome Reporter (PCGR)- variant interpretation report for precision oncology

### Overview

The Personal Cancer Genome Reporter (PCGR) is a stand-alone software package for functional annotation and translation of individual cancer genomes for precision oncology. It interprets both somatic SNVs/InDels and copy number aberrations. The software extends basic gene and variant annotations from the [Ensembl’s Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) with oncology-relevant, up-to-date annotations retrieved flexibly through [vcfanno](https://github.com/brentp/vcfanno), and produces interactive HTML reports intended for clinical interpretation (Figure 1).

![PCGR overview](PCGR_workflow.png)

### Example reports
* [Report for a breast tumor sample (TCGA)](http://folk.uio.no/sigven/tumor_sample.BRCA.0.4.1.pcgr.html)
* [Report for a colon adenocarcinoma sample (TCGA)](http://folk.uio.no/sigven/tumor_sample.COAD.0.4.1.pcgr.html)
* [Report for a lung adenocarcinoma sample (TCGA)](http://folk.uio.no/sigven/tumor_sample.LUAD.0.4.1.pcgr.html)

### PCGR documentation

[![Documentation Status](https://readthedocs.org/projects/pcgr/badge/?version=latest)](http://pcgr.readthedocs.io/en/latest/?badge=latest)

If you use PCGR, please cite our preprint paper:

Sigve Nakken, Ghislain Fournous, Daniel Vodák, Lars Birger Aaasheim, and Eivind Hovig. __Personal Cancer Genome Reporter: Variant Interpretation Report For Precision Oncology__ (2017). bioRxiv. doi:[10.1101/122366](https://doi.org/10.1101/122366)

### Annotation resources included in PCGR (v0.4.1)

* [VEP v85](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor release 85 (GENCODE v19 as the gene reference dataset)
* [COSMIC v81](http://cancer.sanger.ac.uk/cosmic/) - Catalogue of somatic mutations in cancer (May 2017)
* [dBNSFP v3.4](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (March 2017)
* [ICGC v23](https://dcc.icgc.org/) - Somatic mutations discovered in all ICGC (International Cancer Genomics Consortium) tumor cohorts (Dec 2016)
* [ExAC r1](http://exac.broadinstitute.org/) - Germline variant frequencies exome-wide (February 2017)
* [gnomAD r1](http://gnomad.broadinstitute.org/) - Germline variant frequencies exome-wide (March 2017)
* [dbSNP b147](http://www.ncbi.nlm.nih.gov/SNP/) - Database of short genetic variants (April 2016)
* [1000 Genomes Project - phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) - Germline variant frequencies genome-wide (May 2013)
* [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of clinically related variants (March 2017)
* [DoCM](http://docm.genome.wustl.edu) - Database of curated mutations (v3.2, April 2016)
* [CIViC](http://civic.genome.wustl.edu) - Clinical interpretations of variants in cancer (August 3rd 2017)
* [CBMDB](http://www.cancergenomeinterpreter.org/biomarkers) - Cancer Biomarkers database (Feb 8th 2017)
* [IntOGen catalog of driver mutations](https://www.intogen.org/downloads) - (May 2016)
* [Cancer Hotspots](http://cancerhotspots.org) - Resource for statistically significant mutations in cancer (2016)
* [UniProt/SwissProt KnowledgeBase 2016_09](http://www.uniprot.org) - Resource on protein sequence and functional information (June 2017)
* [Pfam v30](http://pfam.xfam.org) - Database of protein families and domains (March 2017)
* [DGIdb](http://dgidb.genome.wustl.edu) - Database of interactions betweeen antineoplastic drugs and human proteins (v2.22, February 2016)
* [TSGene v2.0](http://bioinfo.mc.vanderbilt.edu/TSGene/) - Tumor suppressor/oncogene database (November 2015)

### Getting started

#### STEP 0: Python

A local installation of Python (it has been tested with [version 2.7.13](https://www.python.org/downloads/)) is required to run PCGR. Check that Python is installed by typing `python --version` in a terminal window.

#### STEP 1: Installation of Docker

1. [Install the Docker engine](https://docs.docker.com/engine/installation/) on your preferred platform
2. Test that Docker is running, e.g. by typing `docker ps` or `docker images` in the terminal window
3. Adjust the computing resources dedicated to the Docker, i.e.:
   - Memory: minimum 5GB
   - CPUs: minimum 4
   - [How to - Windows](https://docs.docker.com/docker-for-windows/#advanced) / [How to - Mac OS X](https://docs.docker.com/docker-for-mac/#advanced)  

#### STEP 2: Download PCGR

__August 7th 2017__: New release (0.4.1)

1. Download and unpack the [latest software release (0.4.1)](https://github.com/sigven/pcgr/releases/latest)
2. Download and unpack the data bundle (approx. 17Gb) in the PCGR directory
   * Download [the latest data bundle](https://drive.google.com/file/d/0B8aYD2TJ472mNnpLOFNXdFV3bVE/) from Google Drive to `~/pcgr-X.X` (replace _X.X_ with the version number, e.g `~/pcgr-0.4.1`)
   * Unpack the data bundle, e.g. through the following Unix command: `gzip -dc pcgr.databundle.GRCh37.YYYYMMDD.tgz | tar xvf -`
   * __IMPORTANT NOTE__: The data bundle contains elements (e.g. [COSMIC](http://cancer.sanger.ac.uk/cancergenome/assets/COSMIC_academic_license_march2015.pdf), [IntoGen driver mutations](https://www.intogen.org/downloads), [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP)) which are freely available for non-commercial use only

    A _data/_ folder within the _pcgr-X.X_ software folder should now have been produced
3. Pull the [PCGR Docker image (0.4.1)](https://hub.docker.com/r/sigven/pcgr/) from DockerHub (2.99Gb):
   * `docker pull sigven/pcgr:0.4.1` (PCGR annotation engine)

#### STEP 3: Input preprocessing

The PCGR workflow accepts two types of input files:

  * An unannotated, single-sample VCF file (>= v4.2) with called somatic variants (SNVs/InDels)
  * A copy number segment file

  __NOTE: GRCh37 is currently supported as the reference genome build__

PCGR can be run with either or both of the two input files present.

* We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html)
* If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose)
* Variants used for reporting should be designated as 'PASS' in the VCF FILTER column



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


#### STEP 4: Run example

A tumor sample report is generated by calling the Python script __pcgr.py__, which takes the following arguments and options:

    usage: pcgr.py [-h] [--input_vcf INPUT_VCF] [--input_cna INPUT_CNA]
               [--logR_gain LOGR_GAIN] [--logR_homdel LOGR_HOMDEL]
               [--msig_identify] [--msig_n MSIG_N]
               [--msig_normalization {default,exome,genome,exome2genome}]
               [--msi_predict] [--list_noncoding]
               [--n_vcfanno_proc N_VCFANNO_PROC] [--n_vep_forks N_VEP_FORKS]
               [--tumor_dp_tag TUMOR_DP_TAG] [--tumor_af_tag TUMOR_AF_TAG]
               [--normal_dp_tag NORMAL_DP_TAG] [--normal_af_tag NORMAL_AF_TAG]
               [--call_conf_tag CALL_CONF_TAG] [--force_overwrite] [--version]
               pcgr_dir output_dir sample_id

    Personal Cancer Genome Reporter (PCGR) workflow for clinical interpretation of
    somatic nucleotide variants and copy number aberration segments

    positional arguments:
    pcgr_dir              PCGR base directory with accompanying data directory,
                        e.g. ~/pcgr-0.4.1
    output_dir            Output directory
    sample_id             Tumor sample/cancer genome identifier - prefix for
                        output files

    optional arguments:
    -h, --help            show this help message and exit
    --input_vcf INPUT_VCF
                        VCF input file with somatic query variants
                        (SNVs/InDels). Note: GRCh37 is currently the only
                        reference genome build supported (default: None)
    --input_cna INPUT_CNA
                        Somatic copy number alteration segments (tab-separated
                        values) (default: None)
    --logR_gain LOGR_GAIN
                        Log(2) ratio treshold for reporting copy number gains
                        in genome report (default: 0.8)
    --logR_homdel LOGR_HOMDEL
                        Log(2) ratio treshold for reporting homozygous
                        deletions in genome report (default: -0.8)
    --msig_identify       Identify relative contribution of 30 known mutational
                        signatures (COSMIC) through the deconstructSigs
                        framework (default: False)
    --msig_n MSIG_N       deconstructSigs option: number of mutational
                        signatures (to limit the search to (default: 6)
    --msig_normalization {default,exome,genome,exome2genome}
                        deconstructSigs option: type of trimer count
                        normalization for inference of known mutational
                        signatures, see explanation at
                        https://github.com/raerose01/deconstructSigs (default:
                        default)
    --msi_predict         Perform statistical prediction of tumor microsatellite
                        instability from distribution of somatic indels/snvs
                        in input VCF (default: False)
    --list_noncoding      List non-coding variants in HTML report (default:
                        False)
    --n_vcfanno_proc N_VCFANNO_PROC
                        Number of processes used during vcfanno annotation
                        (default: 4)
    --n_vep_forks N_VEP_FORKS
                        Number of forks (--forks option in VEP) used during
                        VEP annotation (default: 4)
    --tumor_dp_tag TUMOR_DP_TAG
                        Tag in input VCF (INFO column) that denotes total read
                        depth at variant site in tumor sample (default: _na)
    --tumor_af_tag TUMOR_AF_TAG
                        Tag in input VCF (INFO column) that denotes fraction
                        of alternate allele reads in tumor sample (default:
                        _na)
    --normal_dp_tag NORMAL_DP_TAG
                        Tag in input VCF (INFO column) that denotes total read
                        depth at variant site in control/normal sample
                        (default: _na)
    --normal_af_tag NORMAL_AF_TAG
                        Tag in input VCF (INFO column) that denotes fraction
                        of alternate allele reads control/normal sample
                        (default: _na)
    --call_conf_tag CALL_CONF_TAG
                        Tag in input VCF (INFO column) that denotes level of
                        confidence in variant call (default: _na)
    --force_overwrite     By default, the script will fail with an error if any
                        output file already exists. You can force the
                        overwrite of existing result files by using this flag
                        (default: False)
    --version             show program's version number and exit



The _examples_ folder contain sample files from TCGA. A report for a colorectal tumor case can be generated through the following command:

`python pcgr.py --input_vcf ~/pcgr-0.4.1/examples/tumor_sample.COAD.vcf.gz --input_cna ~/pcgr-0.4.1/examples/tumor_sample.COAD.cna.tsv`
`--msig_identify --msi_predict --tumor_af_tag TVAF --tumor_dp_tag TDP --call_conf_tag TAL --list_noncoding`
`~/pcgr-0.4.1 ~/pcgr-0.4.1/examples tumor_sample.COAD`

This command will run the Docker-based PCGR workflow and produce the following output files in the _examples_ folder:

1. __tumor_sample.COAD.pcgr.html__ - An interactive HTML report for clinical interpretation
2. __tumor_sample.COAD.pcgr.vcf.gz__ - VCF file with rich set of annotations for precision oncology
3. __tumor_sample.COAD.pcgr.maf__ - A basic [MAF](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification) file for use as input in downstream analyses with other tools (e.g. [2020plus](https://github.com/KarchinLab/2020plus), MutSigCV)
4. __tumor_sample.COAD.pcgr.snvs_indels.tiers.tsv__ - Tab-separated values file with variants organized according to tiers of functional relevance
5. __tumor_sample.COAD.pcgr.mutational_signatures.tsv__ - Tab-separated values file with estimated contributions by known mutational signatures and associated underlying etiologies
6. __tumor_sample.COAD.pcgr.snvs_indels.biomarkers.tsv__ - Tab-separated values file with clinical evidence items associated with biomarkers for diagnosis, prognosis or drug sensitivity/resistance
7. __tumor_sample.COAD.pcgr.cna_segments.tsv.gz__ - Tab-separated values file with annotations of gene transcripts that overlap with somatic copy number aberrations

## Contact

sigven@ifi.uio.no
