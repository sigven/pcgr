## Annotation resources

### Basic variant consequence annotation
  * [VEP v85](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor release 85 ([GENCODE v19](https://www.gencodegenes.org/releases/19.html) as gene reference database)

###  *Insilico* predictions of effect of coding variants
  * [dBNSFP v3.4](https://sites.google.com/site/jpopgen/dbNSFP) - database of non-synonymous functional predictions (March 2017)
  * [IntOGen catalogs of driver mutations/genes](https://www.intogen.org/downloads) - (May 2016)

###  Variant frequency databases
  * [COSMIC v80](http://cancer.sanger.ac.uk/cosmic/) - catalogue of somatic mutations in cancer (February 2017)
  * [ICGC v23](https://dcc.icgc.org/) - Somatic mutations discovered in all ICGC (International Cancer Genomics Consortium) tumor cohorts (Dec 2016)
  * [ExAC r1](http://exac.broadinstitute.org/) - germline variant frequencies exome-wide (February 2017)
  * [gnomAD r1](http://exac.broadinstitute.org/) - germline variant frequencies exome-wide (March 2017)
  * [dbSNP b147](http://www.ncbi.nlm.nih.gov/SNP/) - database of short genetic variants (April 2016)
  * [1000Genomes phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) - germline variant frequencies genome-wide (May 2013)
  * [Cancer Hotspots](http://cancerhotspots.org) - a resource for statistically significant mutations in cancer (2016)

### Variant databases of clinical utility
  * [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - database of clinically related variants (March 2017)
  * [DoCM](http://docm.genome.wustl.edu) - database of curated mutations (v3.2, April 2016)
  * [CIViC](http://civic.genome.wustl.edu) - clinical interpretations of variants in cancer (April 5th 2017)
  * [CBMDB](http://www.cancergenomeinterpreter.org/biomarkers) - Cancer BioMarkers database (Februay 8th 2017)
  * [DGIdb](http://dgidb.genome.wustl.edu) - database of interactions betweeen antineoplastic drugs and human proteins (v2.22, February 2016)

### Protein domains/functional features
  * [UniProt/SwissProt KnowledgeBase 2016_09](http://www.uniprot.org) - resource on protein sequence and functional information (March 2017)
  * [Pfam v31](http://pfam.xfam.org) - database of protein families and domains (March 2017)

### Cancer gene knowledge bases
  * [TSGene v2.0](http://bioinfo.mc.vanderbilt.edu/TSGene/) - tumor suppressor/oncogene database (November 2015)
  * [Cancer Gene Cencus](http://cancer.sanger.ac.uk/cosmic/) - (February 2017)

## Notes on variant annotation datasets

### Genome mapping

A requirement for all variant annotation datasets used in PCGR is that they have been mapped unambiguously to the human genome (GRCh37). For most datasets this is already the case (i.e. dbSNP, COSMIC, ClinVar etc.). A significant proportion of variants in the annotation datasets related to clinical interpretation, CIViC and CBMDB, is however not mapped to the genome. Whenever possible, we have utilized [TransVar](http://bioinformatics.mdanderson.org/transvarweb/) to identify the actual genomic variants (e.g. _g.chr7:140453136A>T_) that correspond to variants reported with other HGVS nomenclature (e.g. _p.V600E_).

### Other data quality concerns

__Clinical biomarkers__

Clinical biomarkers included in PCGR are limited to the following:

* Markers reported at the variant level (e.g. __BRAF p.V600E__)
* Markers reported at the codon level (e.g. __KRAS p.G12__)
* Markers reported at the exon level (e.g. __KIT exon 11 mutation__)
* Within the [Cancer bioMarkers database (CBMDB)](https://www.cancergenomeinterpreter.org/biomarkers), only markers collected from FDA/NCCN guidelines, scientific literature, and clinical trials are included (markers collected from conference abstracts etc. are not included)

__COSMIC variants__

The COSMIC dataset that is part of the PCGR annotation bundle is the subset of variants that satisfy the following criteria:

* __Mutation somatic status__ is either '_confirmed_somatic_' or '_reported_in_another_cancer_sample_as_somatic_'.
* __Site/histology__ must be known and the sample must come from a malignant tumor (i.e. not polyps/adenomas, which are also found in COSMIC)
