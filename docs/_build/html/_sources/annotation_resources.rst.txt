PCGR annotation resources
-------------------------

Basic variant consequence annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `VEP <http://www.ensembl.org/info/docs/tools/vep/index.html>`__ -
   Variant Effect Predictor release 101 (`GENCODE
   v35 <https://www.gencodegenes.org/human/>`__ as gene reference
   database (v19 for grch37))

*Insilico* predictions of effect of coding variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `dBNSFP <https://sites.google.com/site/jpopgen/dbNSFP>`__ - database
   of non-synonymous functional predictions (v4.1, June 2020)

Variant frequency databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `gnomAD <http://exac.broadinstitute.org/>`__ - germline variant
   frequencies exome-wide (r2.1, October 2018)
-  `dbSNP <http://www.ncbi.nlm.nih.gov/SNP/>`__ - database of short
   genetic variants (build 153)
-  `Cancer Hotspots <http://cancerhotspots.org>`__ - a resource for
   statistically significant mutations in cancer (v2, 2017)
-  `TCGA <https://portal.gdc.cancer.gov/>`__ - somatic mutations
   discovered across 33 tumor type cohorts (release 27.0, October 2020)
-  `ICGC-PCAWG <http://docs.icgc.org/pcawg/>`__ - ICGC Pancancer
   Analysis of Whole Genomes - (release 28, March 17th, 2019)

Variant databases of clinical utility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar/>`__ - database of
   clinically related variants (November 2020)
-  `DoCM <http://docm.genome.wustl.edu>`__ - database of curated
   mutations (v3.2, April 2016)
-  `CIViC <http://civic.genome.wustl.edu>`__ - clinical interpretations
   of variants in cancer (November 18th 2020)
-  `CGI <http://www.cancergenomeinterpreter.org/biomarkers>`__ - Cancer
   Genome Interpreter Cancer Biomarkers Database (CGI) (January 17th
   2018)
-  `ChEMBL <https://www.ebi.ac.uk/chembl/>`__ - database of drugs,
   drug-like small molecules and their targets (ChEMBL_27, May 2020)

Protein domains/functional features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `UniProt/SwissProt KnowledgeBase <http://www.uniprot.org>`__ -
   resource on protein sequence and functional information (2020_05,
   October 2020)
-  `Pfam <http://pfam.xfam.org>`__ - database of protein families and
   domains (v33, May 2020)

Knowledge resources on gene and protein targets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `CancerMine <https://zenodo.org/record/3525385#.XcHblUVKiL4>`__ -
   Literature-mined database of tumor suppressor genes/proto-oncogenes
   (v30, November 2020)
-  `Open Targets Platform <https://www.targetvalidation.org/>`__ -
   Database on disease-target associations, molecularly targeted drugs
   and tractability aggregated from multiple sources (literature,
   pathways, mutations) (2020_09)
-  `TCGA driver genes <https://www.ncbi.nlm.nih.gov/pubmed/29625053>`__
   - predicted cancer driver genes based on application of multiple
   driver gene prediction tools on TCGA pan-cancer cohort

Pathway databases
~~~~~~~~~~~~~~~~~

-  `KEGG PATHWAY Database <http://www.genome.jp/kegg/pathway.htm>`__ -
   November 20th 2020
-  `Oncogenic Signaling Pathways -
   TCGA <https://www.ncbi.nlm.nih.gov/pubmed/29625050>`__ - Sanchez-Vega
   et al., *Cell*, 2018

Notes on variant annotation datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Genome mapping
^^^^^^^^^^^^^^

A requirement for PCGR variant annotation datasets is that variants have
been mapped unambiguously to the reference human genome. For most
datasets this requirement is not an issue (i.e. dbSNP, ClinVar etc.). A
fraction of variants in the annotation datasets related to clinical
interpretation, CIViC and CBMDB, has however not been mapped to the
genome. Whenever possible, we have utilized
`TransVar <http://bioinformatics.mdanderson.org/transvarweb/>`__ to
identify the actual genomic variants (e.g. *g.chr7:140453136A>T*) that
correspond to variants reported at the amino acid level or with other
HGVS nomenclature (e.g. *p.V600E*).

For variants that have been mapped to a specific build (GRCh37/GRCh38),
we have utilized the `crossmap <http://crossmap.sourceforge.net/>`__
package to lift the datasets to the other build.

Data quality
^^^^^^^^^^^^

**Clinical biomarkers**

Clinical biomarkers included in PCGR are limited to the following:

-  Evidence items for specific markers in CIViC must be *accepted*
   (*submitted* evidence items are not considered)
-  Markers reported at the variant level (e.g. **BRAF p.V600E**)
-  Markers reported at the codon level (e.g. **KRAS p.G12**)
-  Markers reported at the exon level (e.g. **KIT exon 11 mutation**)
-  Within the `Cancer bioMarkers database
   (CGI) <https://www.cancergenomeinterpreter.org/biomarkers>`__, only
   markers collected from FDA/NCCN guidelines, scientific literature,
   and clinical trials are included (markers collected from conference
   abstracts etc. are not included)
-  Copy number gains/losses

See also comment on a `closed GitHib
issue <https://github.com/sigven/pcgr/issues/37#issuecomment-391966286>`__

**IMPORTANT NOTE**: The variant consequence reported by CIViC may
deviate from what is reported by PCGR. PCGR picks the variant
consequence according to VEP’s *pick* option (depending on a ranked list
of criteria that can be configured by the user), and this particular
transcript consequence may differ from what has been reported in the
literature.

**Molecularly targeted drugs**

-  For targeted drugs extracted from `Open Targets
   Platform <https://www.targetvalidation.org>`__, we distinguish
   between drugs in late clinical development (phase 3-4), versus those
   in early clinical development (phase 1-2).

**Gene-disease associations**

-  Cancer phenotype associations retrieved from the `Open Targets
   Platform <https://www.targetvalidation.org/>`__ are largely based on
   the `association
   score <https://docs.targetvalidation.org/getting-started/scoring>`__
   developed by the Open Targets Platform, with a couple of extra
   post-processing steps:

   -  Phenotype associations in Open Targets Platform are assembled from
      `a variety of different data
      sources <https://docs.targetvalidation.org/data-sources/data-sources>`__.
      Target-disease associations included in PCGR must be supported by
      **at least two distinct sources**
   -  The weakest associations, here defined as those with an
      association score < 0.4 (scale from 0 to 1), are ommitted
   -  As is done within the Open Targets Platform, association scores
      (for genes) are represented with varying shades of blue: the
      darker the blue, the stronger the association. Variant hits in
      tier 3/4 and the noncoding section are arranged according to this
      association score. If several disease subtypes are associated with
      a gene, the maximum association score is chosen.

**Tumor suppressor genes/proto-oncogenes**

-  Status as oncogenes and/or tumor suppressors genes are done according
   to the following scheme in PCGR:

   -  Five or more publications in the biomedical literature that
      suggests an oncogenic/tumor suppressor role for a given gene (as
      collected from the `CancerMine text-mining
      resource <http://bionlp.bcgsc.ca/cancermine/>`__), **OR**
   -  At least two publications from CancerMine that suggests an
      oncogenic/tumor suppressor role for a given gene **AND** an
      existing record for the same gene as a tumor suppressor/oncogene
      in the `Network of Cancer Genes (NCG) <http://ncg.kcl.ac.uk/>`__
   -  Status as oncogene is ignored if a given gene has three times as
      much (literature evidence) support for a role as a tumor
      suppressor gene (and vice versa)
   -  Oncogenes/tumor suppressor candidates from CancerMine/NCG that are
      found in the `curated list of false positive cancer drivers
      compiled by Bailey et al. (Cell,
      2018) <https://www.ncbi.nlm.nih.gov/pubmed/30096302>`__ have been
      excluded

**TCGA somatic calls**

-  TCGA employs four different variant callers for detection of somatic
   variants (SNVs/InDels): *mutect2, varscan2, somaticsniper and muse*.
   In the TCGA dataset bundled with PCGR, somatic SNVs are restricted to
   those that are detected by at least two independent callers
   (i.e. calls found by a single algorithm are considered low-confident
   and disregarded)
