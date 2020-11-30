About
-----

What is the Personal Cancer Genome Reporter (PCGR)?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Personal Cancer Genome Reporter (PCGR) is a stand-alone software
package for functional annotation and translation of individual cancer
genomes for precision cancer medicine. It interprets both somatic
SNVs/InDels and copy number aberrations. The software extends basic gene
and variant annotations from the `Ensembl’s Variant Effect Predictor
(VEP) <http://www.ensembl.org/info/docs/tools/vep/index.html>`__ with
oncology-relevant, up-to-date annotations retrieved flexibly through
`vcfanno <https://github.com/brentp/vcfanno>`__. Variants are classified
into `tiers of clinical significance <tier_systems.md>`__, and
interactive HTML output reports permits exploration of the final
results.

Example views from the dashboard HTML output:

|image1|

The Personal Cancer Genome Reporter has been developed by scientists
affiliated with the `Norwegian Cancer Genomics
Consortium <http://cancergenomics.no>`__, at the `Institute for Cancer
Research/Oslo University Hospital <http://radium.no>`__.

Example reports
~~~~~~~~~~~~~~~

-  `Cervical cancer sample
   (tumor-only) <http://insilico.hpc.uio.no/pcgr/example_reports/0.9.1/TCGA-FU-A3HZ-01A_TO.pcgr_acmg.grch37.flexdb.html>`__
-  `Lung cancer sample
   (tumor-control) <http://insilico.hpc.uio.no/pcgr/example_reports/0.9.1/TCGA-95-7039-01A.pcgr_acmg.grch37.flexdb.html>`__
-  `Breast cancer sample
   (tumor-control) <http://insilico.hpc.uio.no/pcgr/example_reports/0.9.1/TCGA-EW-A1J5-01A.pcgr_acmg.grch37.flexdb.html>`__
-  `Brain cancer sample
   (tumor-control) <http://insilico.hpc.uio.no/pcgr/example_reports/0.9.1/TCGA-14-0866-01B.pcgr_acmg.grch37.flexdb.html>`__

(to view the rmarkdown-based reports, simply remove *.flexdb.* in the
file names for the flexdashboard reports)

Why use PCGR?
~~~~~~~~~~~~~

The great complexity of acquired mutations in individual tumor genomes
poses a severe challenge for clinical interpretation. There is a general
scarcity of tools that can *i)* systematically interrogate cancer
genomes in the context of diagnostic, prognostic, and therapeutic
biomarkers, *ii)* prioritize and highlight the most important findings,
and *iii)* present the results in a format accessible to clinical
experts. PCGR integrates a comprehensive set of knowledge resources
related to tumor biology and therapeutic biomarkers, both at the gene
and variant level. The application generates a tiered report that will
aid the interpretation of individual cancer genomes in a clinical
setting.

If you use PCGR, please cite our recent publication:

Sigve Nakken, Ghislain Fournous, Daniel Vodák, Lars Birger Aaasheim, Ola
Myklebost, and Eivind Hovig. **Personal Cancer Genome Reporter: variant
interpretation report for precision oncology** (2017). *Bioinformatics*.
34(10):1778–1780.
doi:`10.1093/bioinformatics/btx817 <https://doi.org/10.1093/bioinformatics/btx817>`__

Docker-based technology
~~~~~~~~~~~~~~~~~~~~~~~

The PCGR workflow is developed using the `Docker
technology <https://www.docker.com/what-docker>`__. The software is thus
packaged into isolated containers, in which the installation of all
software libraries/tools and required dependencies have been taken care
of. In addition to the bundled software, in the form of a Docker image,
the workflow only needs to be attached with an `annotation data bundle
for precision oncology <annotation_resources.html>`__.

|image2|

Contact
~~~~~~~

sigven@ifi.uio.no

.. |image1| image:: pcgr_dashboard_views.png
.. |image2| image:: docker-logo50.png
