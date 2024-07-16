## Personal Cancer Genome Reporter (PCGR) <a href="https://sigven.github.io/pcgr/"><img src="pcgrr/man/figures/logo.png" align="right" height="118" width="100"/></a>

[![Conda install ver](https://anaconda.org/pcgr/pcgr/badges/version.svg)](https://anaconda.org/pcgr/pcgr)
[![Conda install lrd](https://anaconda.org/pcgr/pcgr/badges/latest_release_date.svg)](https://anaconda.org/pcgr/pcgr)
[<img src="https://img.shields.io/badge/dockerhub-images-important.svg?logo=LOGO">](https://hub.docker.com/r/sigven/pcgr)

### Overview

The Personal Cancer Genome Reporter (PCGR) is a stand-alone software package for functional annotation and translation of individual tumor genomes for precision cancer medicine. It interprets primarily somatic SNVs/InDels and copy number aberrations, and has additional support for interpretation of bulk RNA-seq expression data. The software classifies variants both with respect to _oncogenicity_, and _actionability_. Interactive HTML output reports allow the user to interrogate the clinical impact of the molecular findings in an individual tumor.

- Variant classification
  - according to *oncogenicity*: evaluating the oncogenic potential of somatic DNA aberrations ([VICC/CGC/ClinGen guidelines](https://pubmed.ncbi.nlm.nih.gov/35101336/))
  - according to *actionability*: mapping the therapeutic, diagnostic, and prognostic implications of somatic DNA aberrations ([AMP/ASCO/CAP guidelines](https://pubmed.ncbi.nlm.nih.gov/27993330/))
- Tumor mutational burden (TMB) estimation
- Mutational signature analysis
- Microsatellite instability (MSI) classification
- RNA-seq support - gene expression outlier detection, sample similarity analysis, and immune contexture profiling

PCGR supports both of the most recent human genome assemblies (GRCh37/GRCh38), and accepts variant calls from both tumor-control and tumor-only sequencing assays. Much of the functionality is intended for whole-exome/whole-genome sequencing assays, but you can also apply PCGR to output from targeted sequencing panels. If you are interested in the interrogation of germline variants and their relation to cancer predisposition, we recommend trying the accompanying tool [Cancer Predisposition Sequencing Reporter (CPSR)](https://github.com/sigven/cpsr).

Example screenshots from the [quarto](https://quarto.org)-based cancer genome report by PCGR:

![PCGR screenshot 1](pcgrr/pkgdown/assets/img/sc2.png)
![PCGR screenshot 2](pcgrr/pkgdown/assets/img/sc1.png)
![PCGR screenshot 3](pcgrr/pkgdown/assets/img/sc3.png)

PCGR originates from the [Norwegian Cancer Genomics Consortium (NCGC)](http://cancergenomics.no), at the [Institute for Cancer Research, Oslo University Hospital, Norway](http://radium.no).

### Top News

- *July 16th 2024*: **2.0.2 release** 
  - patch to ensure correct reference to actionability guidelines
  - [CHANGELOG](http://sigven.github.io/pcgr/articles/CHANGELOG.html)

- *July 7th 2024*: **2.0.1 release** 
  - patch with bug fix for mitochondrial input variants ([pr245](https://github.com/sigven/pcgr/pull/245))
  - [CHANGELOG](http://sigven.github.io/pcgr/articles/CHANGELOG.html)

- *June 2024*: **2.0.0 release**
  - Details in [CHANGELOG](http://sigven.github.io/pcgr/articles/CHANGELOG.html)
  - Massive reference data bundle upgrade, new report layout, oncogenicity classification++
  - Support for Singularity/Apptainer
  - Major data/software updates:
    - Ensembl VEP `v112`
    - ClinVar (2024-06)
    - CIViC (2024-06-21)
    - GENCODE `v46/v19` (GRCh38/GRCh37)
    - CancerMine `v50` (2023-03)
    - UniProt KB `v2024_03`

- *February 2023*: **1.3.0 release**
  - Details in [CHANGELOG](http://sigven.github.io/pcgr/articles/CHANGELOG.html)
  - prioritize protein-coding BIOTYPE csq ([pr201](https://github.com/sigven/pcgr/pull/201))
  - expose `--pcgrr_conda` option to flexibly activate pcgrr env via a non-default pcgrr name
  - `cpsr_validate_input.py`: refactor for efficient custom gene egrep

- *November 2022*: **1.2.0 release**
  -    Keep only autosomal, X, Y, M/MT chromosomes
  -    Import bcftools as dependency

### Example reports

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12752833.svg)](https://doi.org/10.5281/zenodo.12752833)

### Why use PCGR?

The great complexity of acquired mutations in individual tumor genomes poses a severe challenge for clinical interpretation. PCGR aims to be a comprehensive reporting platform that can

- systematically interrogate tumor-specific variants in the context of known therapeutic, diagnostic, and prognostic biomarkers
- highlight genomic aberrations with likely oncogenic potential
- provide a structured and concise summary of the most relevant findings
- present the results in a format accessible to clinical experts

PCGR integrates a [comprehensive set of knowledge resources](https://sigven.github.io/pcgr/articles/annotation_resources.html) related to tumor biology and therapeutic biomarkers, both at the gene, and at the level of individual variants. The software generates a comprehensive molecular interpretation report that supports the translation of individual cancer genomes towards molecularly guided treatment strategies.

### Getting started

- [Installation instructions](https://sigven.github.io/pcgr/articles/installation.html)
- [Run through an example](https://sigven.github.io/pcgr/articles/running.html#example-run)
- Learn more about:

    1) Details regarding [PCGR input files](https://sigven.github.io/pcgr/articles/input.html), and how they should be formatted
    2) Configuration of [key settings](https://sigven.github.io/pcgr/articles/running.html)
    3) The types and contents of [PCGR output files](https://sigven.github.io/pcgr/articles/output.html)
    4) [Variant classifications implemented in PCGR](https://sigven.github.io/pcgr/articles/variant_classification.html)
    5) [Primary tumor sites used in PCGR](https://sigven.github.io/pcgr/articles/primary_tumor_sites.html)
    6) The list of [gene and variant annotation resources](https://sigven.github.io/pcgr/articles/annotation_resources.html) used in PCGR annotation

- [Frequently asked questions (FAQ)](https://sigven.github.io/pcgr/articles/faq.html)

### Citation

If you use PCGR, please cite our publication:

Sigve Nakken, Ghislain Fournous, Daniel Vodák, Lars Birger Aaasheim, Ola Myklebost, and Eivind Hovig. **Personal Cancer Genome Reporter: variant interpretation report for precision oncology** (2017). *Bioinformatics*. 34(10):1778--1780. [doi.org/10.1093/bioinformatics/btx817](https://doi.org/10.1093/bioinformatics/btx817)

## Contact

sigven AT ifi.uio.no
