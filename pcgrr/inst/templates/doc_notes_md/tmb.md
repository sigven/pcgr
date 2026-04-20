Tumor mutational load or mutational burden is a measure of the number of mutations 
within a tumor genome, defined as the total number of mutations per coding 
area of a tumour genome. TMB may serve as a proxy for determining the number 
of neoantigens per tumor, which in turn may have implications for response 
to immunotherapy. PCGR estimates TMB according to three different schemes 
(one to be chosen as displayed in the report):

  1) __TMB_coding_silent__: the same approach as was outlined in a large-scale 
  study of TMB ([Chalmers et al., 2017]({pubmed_url}28420421)), 
  i.e. counting all somatic base substitutions and indels in the protein-coding 
  regions of the sequencing assay, including those at canonical splice sites, 
  and including synonymous alterations. 
  2) __TMB_coding_non_silent__: counting all somatic base substitutions and indels 
  in the protein-coding regions of the sequencing assay, including those at canonical 
  splice sites, but excluding synonymous alterations.
  3) __TMB_missense_only__: missense (non-synonymous) variants variants only, 
  i.e. as employed by [Fernandez et al., 2019]({pubmed_url}31475242)

Numbers obtained with 1), 2) or 3) are next divided by the coding target 
size of the sequencing assay, which is an explicit input parameter to PCGR. 
We encourage users to provide accurate estimates of the target size of the sequencing assay. If the users 
utilize VAF/DP filtering for variants included in the TMB calculation, 
the same cutoffs/thresholds (DP) should ideally be applied when estimating 
the effective assay target size. 
