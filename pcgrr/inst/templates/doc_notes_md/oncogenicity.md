Somatic aberrations (SNV/InDels) are evaluated for oncogenicity through an 
implementation of standard operating procedures proposed by ClinGen/CGC/VICC 
[@Horak2022-uh]. Here, various properties of the variants and genes affected are 
assigned specific scores according to several criteria, both negative and 
positive, pending on whether the properties support an oncogenic or benign 
variant type. These scores are in turn aggregated towards an overall oncogenicity score. 

Note that all properties/criteria provided in the SOP's are _not_ readily 
implemented in PCGR, specifically the ones requiring manual curation or 
expert review (i.e. experimental oncogenic variant evidence, requiring 
support from _in vitro_ or _in vivo_ functional studies (criteria _OS2_)). 
Also, the current source for existing oncogenicity classifications (ClinVar), 
which is needed for implementation of criteria _OS1_, _OM3_, is also 
considerably limited. Considering the current limitations, some oncogenic 
variants are likely to be missed and classified with uncertain significance 
(VUS) by PCGR. We highlight conventional ClinVar classifications (i.e. with 
respect to pathogenicity) alongside the current oncogenicity classifications, 
this may add value to the interpretation for uncertain cases. Furthermore, 
taking into the account the nature of the current implementation, we have 
adopted slightly different score thresholds for variant classifications to 
those proposed originally by [@Horak2022-uh]. We are working to further 
improve the oncogenicity classification in PCGR, and welcome feedback on this matter.

Note also that for somatic copy number aberrations, we showcase potential
oncogenic events as **proto-oncogenes subject to amplifications** (where 
level of amplification is configurable by the user), as well 
as **tumor suppressor genes subject to homozygous deletions**.
