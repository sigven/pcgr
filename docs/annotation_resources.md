## PCGR annotation resources

### Basic variant consequence annotation
  * [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor release 96 ([GENCODE v30](https://www.gencodegenes.org/human/) as gene reference database (v19 for grch37))

###  *Insilico* predictions of effect of coding variants
  * [dBNSFP](https://sites.google.com/site/jpopgen/dbNSFP) - database of non-synonymous functional predictions (v4.0, May 2019)

###  Variant frequency databases
  * [gnomAD](http://exac.broadinstitute.org/) - germline variant frequencies exome-wide (r2.1, October 2018)
  * [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) - database of short genetic variants (b151)
  * [Cancer Hotspots](http://cancerhotspots.org) - a resource for statistically significant mutations in cancer (v2, 2017)
  * [TCGA](https://portal.gdc.cancer.gov/) - somatic mutations discovered across 33 tumor type cohorts (release 16.0, March 2019)
  * [ICGC-PCAWG](http://docs.icgc.org/pcawg/) - ICGC Pancancer Analysis of Whole Genomes - (release 28, March 17th, 2019)

### Variant databases of clinical utility
  * [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - database of clinically related variants (May 2019)
  * [DoCM](http://docm.genome.wustl.edu) - database of curated mutations (v3.2, April 2016)
  * [CIViC](http://civic.genome.wustl.edu) - clinical interpretations of variants in cancer (May 18th 2019)
  * [CBMDB](http://www.cancergenomeinterpreter.org/biomarkers) - Cancer BioMarkers database (January 17th 2018)
  * [DGIdb](http://dgidb.genome.wustl.edu) - database of targeted antineoplastic drugs (v3.0.2, January 2018)
  * [ChEMBL](https://www.ebi.ac.uk/chembl/) - database of drugs, drug-like small molecules and their targets (ChEMBL_25, March 2019)

### Protein domains/functional features
  * [UniProt/SwissProt KnowledgeBase](http://www.uniprot.org) - resource on protein sequence and functional information (2019_04, May 2019)
  * [Pfam](http://pfam.xfam.org) - database of protein families and domains (v32, September 2018)

### Knowledge resources on gene and protein targets
  * [CancerMine](https://zenodo.org/record/2587719#.XJNfS0RKiL4) - Literature-mined database of tumor suppressor genes/proto-oncogenes (v12, May 2019)
  * [Open Targets Platform](https://www.targetvalidation.org/) - Database on disease-target associations and target tractability aggregated from multiple sources (literature, pathways, mutations) (2019_04)
  * [DisGeNET](http://www.disgenet.org) - curated associations between human genes and different tumor types (v6.0, January 2019)
  * [TCGA driver genes](https://www.ncbi.nlm.nih.gov/pubmed/29625053) - predicted cancer driver genes based on application of multiple driver gene prediction tools on TCGA pan-cancer cohort

### Pathway databases
  * [KEGG PATHWAY Database](http://www.genome.jp/kegg/pathway.htm) - March 1st 2019
  * [Oncogenic Signaling Pathways - TCGA](https://www.ncbi.nlm.nih.gov/pubmed/29625050) - Sanchez-Vega et al., *Cell*, 2018

### Notes on variant annotation datasets

#### Genome mapping

A requirement for PCGR variant annotation datasets is that variants have been mapped unambiguously to the reference human genome. For most datasets this requirement is not an issue (i.e. dbSNP, ClinVar etc.). A fraction of variants in the annotation datasets related to clinical interpretation, CIViC and CBMDB, has however not been mapped to the genome. Whenever possible, we have utilized [TransVar](http://bioinformatics.mdanderson.org/transvarweb/) to identify the actual genomic variants (e.g. _g.chr7:140453136A>T_) that correspond to variants reported at the amino acid level or with other HGVS nomenclature (e.g. _p.V600E_).

For variants that have been mapped to a specific build (GRCh37/GRCh38), we have utilized the [crossmap](http://crossmap.sourceforge.net/) package to lift the datasets to the other build.

#### Data quality

__Clinical biomarkers__

Clinical biomarkers included in PCGR are limited to the following:

* Evidence items for specific markers in CIViC must be *accepted* (*submitted* evidence items are not considered)
* Markers reported at the variant level (e.g. __BRAF p.V600E__)
* Markers reported at the codon level (e.g. __KRAS p.G12__)
* Markers reported at the exon level (e.g. __KIT exon 11 mutation__)
* Within the [Cancer bioMarkers database (CBMDB)](https://www.cancergenomeinterpreter.org/biomarkers), only markers collected from FDA/NCCN guidelines, scientific literature, and clinical trials are included (markers collected from conference abstracts etc. are not included)
* Copy number gains/losses

See also comment on a [closed GitHib issue](https://github.com/sigven/pcgr/issues/37#issuecomment-391966286)

__Antineoplastic drugs__

- For drugs extracted from [DGIdb](http://dgidb.genome.wustl.edu), we only include antineoplastic drugs subject to direct interaction with a target (i.e. as recorded in ChEMBL)

__Gene-disease associations__

- For gene-disease associations extracted from DisGeNET, we require a [score](http://www.disgenet.org/web/DisGeNET/menu/dbinfo#score) greater than 0.2 and that the association is suppported by at least one PMID (PubMed article). Associations involving non-cancer type of diseases are not included.
- Cancer phenotype associations retrieved from the [Open Targets platform](https://www.targetvalidation.org/) are largely based on the [association score](https://docs.targetvalidation.org/getting-started/scoring) developed by the Open Targets platform, with a couple of extra post-processing steps:
	- Phenotype associations in OpenTargets are assembled from [20 different data sources](https://docs.targetvalidation.org/data-sources/data-sources). Target-disease associations included in PCGR must be supported by **at least two distinct sources**
	- The weakest associations, here defined as those with an association score < 0.4 (scale from 0 to 1), are ommitted
	- As is done within the Open Targets Platform, association scores (for genes) are represented with varying shades of blue: the darker the blue, the stronger the association. Variant hits in tier 3/4 and the noncoding section are arranged according to this association score. If several disease subtypes are associated with a gene, the maximum association score is chosen.

__Tumor suppressor genes/proto-oncogenes__

- For liteature-derived predictions of tumor suppressor genes/proto-oncogenes from *CancerMine*, we require a *minimum of four PubMed hits*.

__TCGA somatic calls__

- TCGA employs four different variant callers for detection of somatic variants (SNVs/InDels): _mutect2, varscan2, somaticsniper and muse_. In the TCGA dataset bundled with PCGR, somatic SNVs are restricted to those that are detected by at least two independent callers (i.e. calls found by a single algorithm are considered low-confident and disregarded)
