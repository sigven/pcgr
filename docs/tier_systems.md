
## Tier models

### *pcgr_acmg*

This tier model attempts to adopt concensus recommendations by ACMG, as outlined in [Li et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/):

 - *Tier 1: Variants of strong clinical significance* - constitutes variants linked to predictive, prognostic, or diagnostic biomarkers in the [CIViC database](http://civic.genome.wustl.edu) and the [Cancer Biomarkers Database](https://www.cancergenomeinterpreter.org/biomarkers) that are
	 - A: found within the same tumor type/class as specified by the user, AND
	 - B: of strong clinical evidence (i.e. part of guidelines, validated or discovered in late clinical trials)
 - *Tier 2: Variants of potential clinical significance* - constitutes other variants linked to predictive, prognostic, or diagnostic biomarkers in the [CIViC database](http://civic.genome.wustl.edu) and the [Cancer Biomarkers Database](https://www.cancergenomeinterpreter.org/biomarkers) that are either
	 - A: of strong clinical evidence in other tumor types/classes than the one specified by the user, OR
	 - B: of weak clinical evidence (early trials, case reports etc.) in the same tumor type/class as specified by the user
 - *Tier 3: Variants of uncertain clinical significance* - includes other coding variants found in oncogenes or tumor suppressor genes
 - *Tier 4* - includes other coding variants

 For **copy number aberrations**, aberrations linked to Tier 1 & 2 are displayed, within following sections in the HTML report:

  * *Copy number aberrations as biomarkers: Aberrations of strong clinical significance* (Tier 1)
  * *Copy number aberrations as biomarkers: Aberrations of potential clinical significance* (Tier 2)
