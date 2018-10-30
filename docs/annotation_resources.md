## Annotation resources

### Basic variant consequence annotation
  * [VEP v94](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor release 94 ([GENCODE v28](https://www.gencodegenes.org/releases/28.html) as gene reference database (v19 for grch37))

###  *Insilico* predictions of effect of coding variants
  * [dBNSFP v3.5](https://sites.google.com/site/jpopgen/dbNSFP) - database of non-synonymous functional predictions (August 2017)
  * [IntOGen catalogs of driver mutations/genes](https://www.intogen.org/downloads) - (May 2016)

###  Variant frequency databases
  * [gnomAD r2](http://exac.broadinstitute.org/) - germline variant frequencies exome-wide (October 2017)
  * [dbSNP b151](http://www.ncbi.nlm.nih.gov/SNP/) - database of short genetic variants (build 150 for grch37)
  * [1000Genomes phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) - germline variant frequencies genome-wide (May 2013)
  * [Cancer Hotspots](http://cancerhotspots.org) - a resource for statistically significant mutations in cancer (v2, 2017)
  * [TCGA release 13.0](https://portal.gdc.cancer.gov/) - somatic mutations discovered across 33 tumor type cohorts (The Cancer Genome Atlas)
  * [ICGC-PCAWG](http://docs.icgc.org/pcawg/) - ICGC Pancancer Analysis of Whole Genomes - release 27, April 30th, 2018

### Variant databases of clinical utility
  * [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - database of clinically related variants (October 2018)
  * [DoCM](http://docm.genome.wustl.edu) - database of curated mutations (v3.2, April 2016)
  * [CIViC](http://civic.genome.wustl.edu) - clinical interpretations of variants in cancer (October 17th 2018)
  * [CBMDB](http://www.cancergenomeinterpreter.org/biomarkers) - Cancer BioMarkers database (January 17th 2018)
  * [DGIdb](http://dgidb.genome.wustl.edu) - database of targeted antineoplastic drugs (v3.0.2, January 2018)

### Protein domains/functional features
  * [UniProt/SwissProt KnowledgeBase 2018_09](http://www.uniprot.org) - resource on protein sequence and functional information (October 2018)
  * [Pfam v32](http://pfam.xfam.org) - database of protein families and domains (September 2018)

### Cancer gene knowledge bases
  * [CancerMine](https://zenodo.org/record/1336650#.W9do9WJKiL4) - Literature-mined database of tumor suppressor genes/proto-oncogenes (release5, October 2018)
  * [DisGeNET v5.0](http://www.disgenet.org) - curated associations between human genes and different tumor types
  * [TCGA driver genes](https://www.ncbi.nlm.nih.gov/pubmed/29625053) - predicted cancer driver genes based on application of multiple driver gene prediction tools on TCGA pan-cancer cohort

### Pathway databases
  * [KEGG PATHWAY Database](http://www.genome.jp/kegg/pathway.htm) - August 21st 2018

### Notes on variant annotation datasets

#### Genome mapping

A requirement for PCGR variant annotation datasets is that variants have been mapped unambiguously to the reference human genome. For most datasets this requirement is not an issue (i.e. dbSNP, ClinVar etc.). A fraction of variants in the annotation datasets related to clinical interpretation, CIViC and CBMDB, has however not been mapped to the genome. Whenever possible, we have utilized [TransVar](http://bioinformatics.mdanderson.org/transvarweb/) to identify the actual genomic variants (e.g. _g.chr7:140453136A>T_) that correspond to variants reported at the amino acid level or with other HGVS nomenclature (e.g. _p.V600E_).

For variants that have been mapped to a specific build (GRCh37/GRCh38), we have utilized the [crossmap](http://crossmap.sourceforge.net/) package to lift the datasets to the other build.

#### Data quality

__Clinical biomarkers__

Clinical biomarkers included in PCGR are limited to the following:

* Markers in CIViC must be *accepted* (*submitted* evidence items are not considered)
* Markers reported at the variant level (e.g. __BRAF p.V600E__)
* Markers reported at the codon level (e.g. __KRAS p.G12__)
* Markers reported at the exon level (e.g. __KIT exon 11 mutation__)
* Within the [Cancer bioMarkers database (CBMDB)](https://www.cancergenomeinterpreter.org/biomarkers), only markers collected from FDA/NCCN guidelines, scientific literature, and clinical trials are included (markers collected from conference abstracts etc. are not included)
* Copy number gains/losses

See also comment on a [closed GitHib issue](https://github.com/sigven/pcgr/issues/37#issuecomment-391966286)

__Antineoplastic drugs__

- For drugs extracted from [DGIdb](http://dgidb.genome.wustl.edu), we only include antineoplastic drugs subject to direct interaction with a target (i.e. as recorded in ChEMBL)

__Gene-disease associations__

- For gene-disease associations extracted from DisGeNET 5.0, we require a [score](http://www.disgenet.org/web/DisGeNET/menu/dbinfo#score) greater than 0.2 and that the association is suppported by at least one PMID (PubMed article). Associations involving non-cancer type of diseases are not included.

__TCGA somatic calls__

- TCGA employs four different variant callers for detection of somatic variants (SNVs/InDels): _mutect2, varscan2, somaticsniper and muse_. In the TCGA dataset bundled with PCGR, somatic SNVs are restricted to those that are detected by at least two independent callers (i.e. calls found by a single algorithm are considered low-confident and disregarded)
