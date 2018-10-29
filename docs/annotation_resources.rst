Annotation resources
--------------------

Basic variant consequence annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `VEP v94 <http://www.ensembl.org/info/docs/tools/vep/index.html>`__ -
   Variant Effect Predictor release 93 (`GENCODE
   v28 <https://www.gencodegenes.org/releases/28.html>`__ as gene
   reference database (v19 for grch37))

*Insilico* predictions of effect of coding variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `dBNSFP v3.5 <https://sites.google.com/site/jpopgen/dbNSFP>`__ -
   database of non-synonymous functional predictions (August 2017)
-  `IntOGen catalogs of driver
   mutations/genes <https://www.intogen.org/downloads>`__ - (May 2016)

Variant frequency databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `gnomAD r2 <http://exac.broadinstitute.org/>`__ - germline variant
   frequencies exome-wide (October 2017)
-  `dbSNP b151 <http://www.ncbi.nlm.nih.gov/SNP/>`__ - database of short
   genetic variants (build 150 for grch37)
-  `1000Genomes
   phase3 <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>`__
   - germline variant frequencies genome-wide (May 2013)
-  `Cancer Hotspots <http://cancerhotspots.org>`__ - a resource for
   statistically significant mutations in cancer (v2, 2017)
-  `TCGA release 13.0 <https://portal.gdc.cancer.gov/>`__ - somatic
   mutations discovered across 33 tumor type cohorts (The Cancer Genome
   Atlas)
-  `ICGC-PCAWG <http://docs.icgc.org/pcawg/>`__ - ICGC Pancancer
   Analysis of Whole Genomes - release 27, April 30th, 2018

Variant databases of clinical utility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar/>`__ - database of
   clinically related variants (October 2018)
-  `DoCM <http://docm.genome.wustl.edu>`__ - database of curated
   mutations (v3.2, April 2016)
-  `CIViC <http://civic.genome.wustl.edu>`__ - clinical interpretations
   of variants in cancer (October 17th 2018)
-  `CBMDB <http://www.cancergenomeinterpreter.org/biomarkers>`__ -
   Cancer BioMarkers database (January 17th 2018)
-  `DGIdb <http://dgidb.genome.wustl.edu>`__ - database of targeted
   antineoplastic drugs (v3.0.2, January 2018)

Protein domains/functional features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `UniProt/SwissProt KnowledgeBase 2018_09 <http://www.uniprot.org>`__
   - resource on protein sequence and functional information (October
   2018)
-  `Pfam v32 <http://pfam.xfam.org>`__ - database of protein families
   and domains (September 2018)

Cancer gene knowledge bases
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `CancerMine <https://bioinfo.uth.edu/TSGene/>`__ - Literature-mined
   database of tumor suppressor genes/proto-oncogenes (release5, October
   2018)
-  `DisGeNET v5.0 <http://www.disgenet.org>`__ - curated associations
   between human genes and different tumor types
-  `TCGA driver genes <https://www.ncbi.nlm.nih.gov/pubmed/29625053>`__
   - predicted cancer driver genes based on application of multiple
   driver gene prediction tools on TCGA pan-cancer cohort

Pathway databases
~~~~~~~~~~~~~~~~~

-  `KEGG PATHWAY Database <http://www.genome.jp/kegg/pathway.htm>`__ -
   August 21st 2018

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
identify the actual genomic variants (e.g. *g.chr7:140453136A>T*) that
correspond to variants reported at the amino acid level or with other
HGVS nomenclature (e.g. *p.V600E*).

For variants that have been mapped to a specific build (GRCh37/GRCh38),
we have utilized the `crossmap <http://crossmap.sourceforge.net/>`__
package to lift the datasets to the other build.

Data quality
^^^^^^^^^^^^

**Clinical biomarkers**

Clinical biomarkers included in PCGR are limited to the following:

-  Markers in CIViC must be *accepted* (*submitted* evidence items are
   not considered)
-  Markers reported at the variant level (e.g. **BRAF p.V600E**)
-  Markers reported at the codon level (e.g. **KRAS p.G12**)
-  Markers reported at the exon level (e.g. **KIT exon 11 mutation**)
-  Within the `Cancer bioMarkers database
   (CBMDB) <https://www.cancergenomeinterpreter.org/biomarkers>`__, only
   markers collected from FDA/NCCN guidelines, scientific literature,
   and clinical trials are included (markers collected from conference
   abstracts etc. are not included)
-  Copy number gains/losses

See also comment on a `closed GitHib
issue <https://github.com/sigven/pcgr/issues/37#issuecomment-391966286>`__

**Antineoplastic drugs**

-  For drugs extracted from `DGIdb <http://dgidb.genome.wustl.edu>`__,
   we only include antineoplastic drugs subject to direct interaction
   with a target (i.e. as recorded in ChEMBL)

**Gene-disease associations**

-  For gene-disease associations extracted from DisGeNET 5.0, we require
   a `score <http://www.disgenet.org/web/DisGeNET/menu/dbinfo#score>`__
   greater than 0.2 and that the association is suppported by at least
   one PMID (PubMed article). Associations involving non-cancer type of
   diseases are not included.

**TCGA somatic calls**

-  TCGA employs four different variant callers for detection of somatic
   variants (SNVs/InDels): *mutect2, varscan2, somaticsniper and muse*.
   In the TCGA dataset bundled with PCGR, somatic SNVs are restricted to
   those that are detected by at least two independent callers
   (i.e. calls found by a single algorithm are considered low-confident
   and disregarded)
