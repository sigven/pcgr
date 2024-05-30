---
editor_options: 
  markdown: 
    wrap: 72
---


# Personal Cancer Genome Reporter <a href="https://sigven.github.io/pcgr/"><img src="man/figures/logo.png" align="right" height="106" width="90"/></a>

<br>

The Personal Cancer Genome Reporter (PCGR) is a stand-alone software package for functional annotation and translation of individual tumor genomes for precision cancer medicine. It interprets primarily somatic SNVs/InDels and copy number aberrations, and also have support for interpretation of bulk RNA-seq expression data. The software [classifies variants](articles/variant_classification.html) both with respect to _oncogenicity_, and _actionability_. Interactive HTML output reports allow the user to interrogate the clinical impact of the molecular findings in an individual tumor.

Example views from the dashboard HTML output:

![](img/pcgr_dashboard_views.png)

PCGR originates from the [Norwegian Cancer Genomics Consortium (NCGC)](http://cancergenomics.no), at the [Institute for Cancer Research, Oslo University Hospital, Norway](http://radium.no).

## Example reports

-   [Cervical cancer sample (tumor-control)](http://insilico.hpc.uio.no/pcgr/example_reports/latest/cervix_tumor_control.grch37.flexdb.html)
-   [Stomach cancer sample (tumor-control)](http://insilico.hpc.uio.no/pcgr/example_reports/latest/esophagus_stomach_tumor_control.grch37.flexdb.html)
-   [Breast cancer sample (tumor-only)](http://insilico.hpc.uio.no/pcgr/example_reports/latest/breast_tumor_only.grch37.flexdb.html)

(to view the rmarkdown-based reports, simply remove *.flexdb.* in the file names for the flexdashboard reports)

## Why use PCGR?

The great complexity of acquired mutations in individual tumor genomes poses a severe challenge for clinical interpretation. PCGR aims to be a comprehensive reporting platform that can

-   systematically interrogate tumor-specific variants in the context of known therapeutic and prognostic biomarkers
-   highlight genomic aberrations with likely oncogenic potential
-   provide a structured and concise summary of the most relevant findings
-   present the results in a format accessible to clinical experts

PCGR integrates a [comprehensive set of knowledge resources](articles/annotation_resources.html) related to tumor biology and therapeutic biomarkers, both at the gene and variant level. The software generates a comprehensive molecular interpretation report that supports the translation of individual cancer genomes towards molecularly guided treatment strategies.

## Getting started

-   [Installation instructions](articles/installation.html)
-   [Run through an example](articles/running.html#example-run)
-   Learn more about
    1)  Details regarding [PCGR input files](articles/input.html), and how they should be formatted
    2)  Configuration of [key settings](articles/running.html#key-settings)
    3)  The types and contents of [PCGR output files](articles/output.html)
    4)  [Variant classifications implemented in PCGR](articles/variant_classification.html)
    5)  [Primary tumor sites used in PCGR](articles/primary_tumor_sites.html)
    6)  The list of [gene and variant annotation resources](articles/annotation_resources.html) available in PCGR  
    <br>
-   [Frequenty asked questions (FAQ)](articles/faq.html)

## Citation

If you use PCGR, please cite our publication:

Sigve Nakken, Ghislain Fournous, Daniel Vodák, Lars Birger Aaasheim, Ola Myklebost, and Eivind Hovig. **Personal Cancer Genome Reporter: variant interpretation report for precision oncology** (2017). *Bioinformatics*. 34(10):1778--1780. [doi.org/10.1093/bioinformatics/btx817](https://doi.org/10.1093/bioinformatics/btx817)

## Contact

[sigven\@ifi.uio.no](mailto:sigven@ifi.uio.no){.email}
