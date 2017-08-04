Annotation resources
--------------------

Basic variant consequence annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `VEP v85 <http://www.ensembl.org/info/docs/tools/vep/index.html>`__ -
   Variant Effect Predictor release 85 (`GENCODE
   v19 <https://www.gencodegenes.org/releases/19.html>`__ as gene
   reference database)

*Insilico* predictions of effect of coding variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `dBNSFP v3.4 <https://sites.google.com/site/jpopgen/dbNSFP>`__ -
   database of non-synonymous functional predictions (March 2017)
-  `IntOGen catalogs of driver
   mutations/genes <https://www.intogen.org/downloads>`__ - (May 2016)

Variant frequency databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `COSMIC v81 <http://cancer.sanger.ac.uk/cosmic/>`__ - catalogue of
   somatic mutations in cancer (May 2017)
-  `ICGC v23 <https://dcc.icgc.org/>`__ - Somatic mutations discovered
   in all ICGC (International Cancer Genomics Consortium) tumor cohorts
   (Dec 2016)
-  `ExAC r1 <http://exac.broadinstitute.org/>`__ - germline variant
   frequencies exome-wide (February 2017)
-  `gnomAD r1 <http://exac.broadinstitute.org/>`__ - germline variant
   frequencies exome-wide (March 2017)
-  `dbSNP b147 <http://www.ncbi.nlm.nih.gov/SNP/>`__ - database of short
   genetic variants (April 2016)
-  `1000Genomes
   phase3 <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>`__
   - germline variant frequencies genome-wide (May 2013)
-  `Cancer Hotspots <http://cancerhotspots.org>`__ - a resource for
   statistically significant mutations in cancer (2016)

Variant databases of clinical utility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar/>`__ - database of
   clinically related variants (March 2017)
-  `DoCM <http://docm.genome.wustl.edu>`__ - database of curated
   mutations (v3.2, April 2016)
-  `CIViC <http://civic.genome.wustl.edu>`__ - clinical interpretations
   of variants in cancer (August 3rd 2017)
-  `CBMDB <http://www.cancergenomeinterpreter.org/biomarkers>`__ -
   Cancer BioMarkers database (Februay 8th 2017)
-  `DGIdb <http://dgidb.genome.wustl.edu>`__ - database of interactions
   betweeen antineoplastic drugs and human proteins (v2.22, February
   2016)

Protein domains/functional features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `UniProt/SwissProt KnowledgeBase 2016\_09 <http://www.uniprot.org>`__
   - resource on protein sequence and functional information (June 2017)
-  `Pfam v31 <http://pfam.xfam.org>`__ - database of protein families
   and domains (March 2017)

Cancer gene knowledge bases
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `TSGene v2.0 <http://bioinfo.mc.vanderbilt.edu/TSGene/>`__ - tumor
   suppressor/oncogene database (November 2015)
-  `Cancer Gene Cencus <http://cancer.sanger.ac.uk/cosmic/>`__ - (v81,
   May 2017)

Notes on variant annotation datasets
------------------------------------

Genome mapping
~~~~~~~~~~~~~~

A requirement for all variant annotation datasets used in PCGR is that
they have been mapped unambiguously to the human genome (GRCh37). For
most datasets this is already the case (i.e. dbSNP, COSMIC, ClinVar
etc.). A significant proportion of variants in the annotation datasets
related to clinical interpretation, CIViC and CBMDB, is however not
mapped to the genome. Whenever possible, we have utilized
`TransVar <http://bioinformatics.mdanderson.org/transvarweb/>`__ to
identify the actual genomic variants (e.g. *g.chr7:140453136A>T*) that
correspond to variants reported with other HGVS nomenclature (e.g.
*p.V600E*).

Other data quality concerns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Clinical biomarkers**

Clinical biomarkers included in PCGR are limited to the following:

-  Markers reported at the variant level (e.g. **BRAF p.V600E**)
-  Markers reported at the codon level (e.g. **KRAS p.G12**)
-  Markers reported at the exon level (e.g. **KIT exon 11 mutation**)
-  Within the `Cancer bioMarkers database
   (CBMDB) <https://www.cancergenomeinterpreter.org/biomarkers>`__, only
   markers collected from FDA/NCCN guidelines, scientific literature,
   and clinical trials are included (markers collected from conference
   abstracts etc. are not included)

**COSMIC variants**

The COSMIC dataset that is part of the PCGR annotation bundle is the
subset of variants that satisfy the following criteria:

-  **Mutation somatic status** is either '*confirmed\_somatic*' or
   '*reported\_in\_another\_cancer\_sample\_as\_somatic*'.
-  **Site/histology** must be known and the sample must come from a
   malignant tumor (i.e. not polyps/adenomas, which are also found in
   COSMIC)
