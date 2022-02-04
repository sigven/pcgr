## Personal Cancer Genome Reporter (PCGR) - variant interpretation for precision cancer medicine

### Overview

The Personal Cancer Genome Reporter (PCGR) is a stand-alone software package for translation of individual tumor genomes for precision cancer medicine.

PCGR interprets primarily **somatic SNVs/InDels and copy number aberrations**. The software extends basic gene and variant annotations from the [Ensembl's Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) with oncology-relevant, up-to-date annotations retrieved flexibly through [vcfanno](https://github.com/brentp/vcfanno), and produces interactive HTML reports intended for clinical interpretation. PCGR can perform multiple types of analyses, including:

- Somatic variant classification (ACMG/AMP)
	- mapping the therapeutic and prognostic implications of somatic DNA aberrations
- Tumor mutational burden (TMB) estimation
- Tumor-only analysis (variant filtering)
- Mutational signature analysis
- Kataegis detection
- Microsatellite instability (MSI) classification

If you want to generate a personal report with respect to cancer predisposition for germline variants, we recommend trying the accompanying tool [Cancer Predisposition Sequencing Reporter (CPSR)](https://github.com/sigven/cpsr).

![PCGR overview](pcgrr/pkgdown/assets/img/pcgr_dashboard_views.png)

### News

-   *February 2022*: **x.x.x release**

    -   Complete restructure of Python and R components. Installation now relies on two separate [conda](https://docs.conda.io/en/latest/) packages, `pcgr` (Python component) and `pcgrr` (R component). Direct Docker support remains, with the Dockerfile simplified to rely exclusively on the installation of the above Conda packages. Significant contributon by the great [@pdiakumis](https://github.com/pdiakumis)
    -   VCF validation step removed. Feedback from users suggested that Ensembl's `vcf-validator` was often too stringent so its use has been deprecated. The `--no_vcf_validate` option remains for backwards compatibility.
    -   New documentation site (<https://sigven.github.io/pcgr>)
    -   Data bundle updates (CIViC, ClinVar, Open Targets Platform, CancerMine, UniProt KB, Pfam)
    -   [CHANGELOG](http://sigven.github.io/pcgr/articles/CHANGELOG.html)

-   *June 30th 2021*: **0.9.2 release**

    -   Data bundle updates (CIViC, ClinVar, CancerMine, UniProt KB, PFAM)
    -   Software upgrades: VEP (104), R v4.1/BioConductor 3.13
    -   **NEW**: TOML configuration removed - all options to PCGR are now command-line based
    -   **NEW**: Feed PCGR with a [CPSR report](https://github.com/sigven/cpsr) to view key germline findings in the tumor report
    -   [CHANGELOG](http://pcgr.readthedocs.io/en/latest/CHANGELOG.html)
    -   Planned for next release: Support for analysis of RNA fusions

-   *November 30th 2020*: **0.9.1 release**

    -   Data bundle updates (CIViC, ClinVar, CancerMine, UniProt KB)
    -   [CHANGELOG](http://pcgr.readthedocs.io/en/latest/CHANGELOG.html)

### Example reports

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5045309.svg)](https://doi.org/10.5281/zenodo.5045309)

### Getting started

-   [Installation instructions](https://sigven.github.io/pcgr/articles/installation.html)
-   [Run through an example](https://sigven.github.io/pcgr/articles/running.html#example-run)
-   Learn more about
	*  Details regarding [PCGR input files](https://sigven.github.io/pcgr/articles/input.html), and how they should be formatted
	*  Configuration of [key settings](https://sigven.github.io/pcgr/articles/running.html#key-settings)
	* The types and contents of [PCGR output files](https://sigven.github.io/pcgr/articles/output.html)
	* The [variant tier system](https://sigven.github.io/pcgr/articles/variant_classification.html) implemented in PCGR
	* The list of [gene and variant annotation resources](https://sigven.github.io/pcgr/articles/virtual_panels.html) used in PCGR annotation

-   [Frequenty asked questions (FAQ)](https://sigven.github.io/pcgr/articles/faq.html)

### Citation

If you use PCGR, please cite the publication:

Sigve Nakken, Ghislain Fournous, Daniel Vodák, Lars Birger Aaasheim, Ola Myklebost, and Eivind Hovig. **Personal Cancer Genome Reporter: variant interpretation report for precision oncology** (2017). *Bioinformatics*. 34(10):1778–1780. [doi:10.1093/bioinformatics/btx817](https://doi.org/10.1093/bioinformatics/btx817)

## Contact

sigven AT ifi.uio.no
