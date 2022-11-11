## Personal Cancer Genome Reporter (PCGR) <a href="https://sigven.github.io/pcgr/"><img src="pcgrr/man/figures/logo.png" align="right" height="118" width="100"/></a>

### Overview

The Personal Cancer Genome Reporter (PCGR) is a stand-alone software package for translation of individual tumor genomes for precision cancer medicine.

PCGR interprets primarily **somatic SNVs/InDels and copy number aberrations**. The software extends basic gene and variant annotations from the [Ensembl's Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) with oncology-relevant, up-to-date annotations retrieved flexibly through [vcfanno](https://github.com/brentp/vcfanno), and produces interactive HTML reports intended for clinical interpretation. PCGR can perform multiple types of analyses, including:

-   Somatic variant classification (ACMG/AMP)
    -   mapping the therapeutic and prognostic implications of somatic DNA aberrations
-   Tumor mutational burden (TMB) estimation
-   Tumor-only analysis (variant filtering)
-   Mutational signature analysis
-   Kataegis detection
-   Microsatellite instability (MSI) classification

If you want to interrogate germline variants and their relation to cancer predisposition, we recommend trying the accompanying tool [Cancer Predisposition Sequencing Reporter (CPSR)](https://github.com/sigven/cpsr).

![PCGR overview](pcgrr/pkgdown/assets/img/pcgr_dashboard_views.png)

### News

-   *November 2022*: **1.2.0 release**
    -    Keep only autosomal, X, Y, M/MT chromosomes
    -    Import bcftools as dependency


-   *October 2022*: **1.1.0 release**

    -   Remove Docker command wrappers
    -   Deprecate `--no_docker` and `--docker_uid` CLI options
    -   Merged PRs [pr192](https://github.com/sigven/pcgr/pull/192), [pr193](https://github.com/sigven/pcgr/pull/193), [pr194](https://github.com/sigven/pcgr/pull/194), [pr196](https://github.com/sigven/pcgr/pull/196).
    -   See [CHANGELOG](http://sigven.github.io/pcgr/articles/CHANGELOG.html) for a few more changes.

-   *May 2022*: **1.0.3 release**

    -   Merged [PR #191](https://github.com/sigven/pcgr/pull/191)

-   *March 2022*: **1.0.2 release**

    -   Fixed [CPSR issue #44](https://github.com/sigven/cpsr/issues/44)

-   *March 2022*: **1.0.1 release**

    -   Fixed bug for huge input sets that cause JSON output crash
        -   huge input variant sets (WGS) are now reduced prior to reporting with R, i.e. exclusion of intronic and intergenic variants, as well as upstream/downstream gene variants ([#178](https://github.com/sigven/pcgr/issues/178)).
    -   Fixed bug for cases where mutational signature analysis reports \> 18 different aetiologies after fitting ([#187](https://github.com/sigven/pcgr/issues/187)).
    -   [CHANGELOG](http://sigven.github.io/pcgr/articles/CHANGELOG.html)

-   *February 2022*: **1.0.0 release**

    -   Complete restructure of Python and R components. Installation now relies on two separate [conda](https://docs.conda.io/en/latest/) packages, `pcgr` (Python component) and `pcgrr` (R component). Direct Docker support remains, with the Dockerfile simplified to rely exclusively on the installation of the above Conda packages. Significant contributon by the great [\@pdiakumis](https://github.com/pdiakumis)
    -   VCF validation step removed. Feedback from users suggested that Ensembl's `vcf-validator` was often too stringent so its use has been deprecated. The `--no_vcf_validate` option remains for backwards compatibility.
    -   New documentation site (<https://sigven.github.io/pcgr>)
    -   Data bundle updates (CIViC, ClinVar, Open Targets Platform, CancerMine, UniProt KB, Pfam)
    -   [CHANGELOG](http://sigven.github.io/pcgr/articles/CHANGELOG.html)

-   *June 30th 2021*: **0.9.2 release**

    -   Data bundle updates (CIViC, ClinVar, CancerMine, UniProt KB, PFAM)
    -   Software upgrades: VEP (104), R v4.1/BioConductor 3.13
    -   **NEW**: TOML configuration removed - all options to PCGR are now command-line based
    -   **NEW**: Feed PCGR with a [CPSR report](https://github.com/sigven/cpsr) to view key germline findings in the tumor report
    -   [CHANGELOG](http://sigven.github.io/pcgr/articles/CHANGELOG.html)
    -   Planned for next release: Support for analysis of RNA fusions

-   *November 30th 2020*: **0.9.1 release**

    -   Data bundle updates (CIViC, ClinVar, CancerMine, UniProt KB)
    -   [CHANGELOG](http://sigven.github.io/pcgr/articles/CHANGELOG.html)

### Example reports

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6275299.svg)](https://doi.org/10.5281/zenodo.6275299)

### Getting started

-   [Installation instructions](https://sigven.github.io/pcgr/articles/installation.html)
-   [Run through an example](https://sigven.github.io/pcgr/articles/running.html#example-run)
-   Learn more about
    -   Details regarding [PCGR input files](https://sigven.github.io/pcgr/articles/input.html), and how they should be formatted
    -   Configuration of [key settings](https://sigven.github.io/pcgr/articles/running.html#key-settings)
    -   The types and contents of [PCGR output files](https://sigven.github.io/pcgr/articles/output.html)
    -   The [variant tier system](https://sigven.github.io/pcgr/articles/variant_classification.html) implemented in PCGR
    -   The list of [gene and variant annotation resources](https://sigven.github.io/pcgr/articles/virtual_panels.html) used in PCGR annotation
-   [Frequenty asked questions (FAQ)](https://sigven.github.io/pcgr/articles/faq.html)

### Citation

If you use PCGR, please cite the publication:

Sigve Nakken, Ghislain Fournous, Daniel Vodák, Lars Birger Aaasheim, Ola Myklebost, and Eivind Hovig. **Personal Cancer Genome Reporter: variant interpretation report for precision oncology** (2017). *Bioinformatics*. 34(10):1778--1780. [doi:10.1093/bioinformatics/btx817](https://doi.org/10.1093/bioinformatics/btx817)

## Contact

sigven AT ifi.uio.no
