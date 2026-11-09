Homologous recombination deficiency (HRD) can be reported through two paths, 
shown side by side when both are available:

* __External HRD assay/tool__ (`--input_hrd`): a status and/or score from an external 
assay or tool, with the type of score (e.g. genomic scar, mutational signature-based) 
stated in the input file. A status given in the file is used as is, unless a threshold 
on the score is provided (`--input_hrd_threshold`, required if the file has no status), 
from which PCGR derives the status. PCGR sets no default threshold, as these are specific to the assay/tool. 
* __PCGR genomic scar score__ (`--estimate_hrd`): the sum of three scar measures (HRD-LOH, LST 
and TAI), computed on autosomes from allele-specific copy number segments (`--input_cna`) 
following the [scarHRD](https://github.com/sztup/scarHRD) algorithm 
([Sztupinszki et al., 2018]({pubmed_url}29978035)). It requires genome-wide segments 
(i.e. _not_ targeted assays) that were called with a matched normal sample. This applies even 
when PCGR runs with `--tumor_only`, which concerns only the SNV/InDel calls. Without a matched 
normal, germline runs of homozygosity can be mistaken for somatic LOH and inflate the score. The score is intended for research use only, and PCGR 
does not apply a threshold to it - an HRD status is thus reported only for external results. 
