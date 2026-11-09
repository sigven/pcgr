Microsatellite instability (MSI) is the result of impaired DNA mismatch repair and 
constitutes a cellular phenotype of clinical significance in many cancer types, 
most prominently colorectal cancers, stomach cancers, endometrial cancers, 
and ovarian cancers ([Cortes-Ciriano et al., 2017]({pubmed_url}28585546)). 
MSI status can be reported through two paths, shown side by side when both are available:

* __External MSI caller__ (`--input_msi`): a result from a dedicated MSI caller 
(e.g. MSIsensor-pro). A status reported by the caller is used as is. For callers 
that report a score only (MSIsensor-pro), PCGR derives the status from a threshold 
on the score (default: 20% unstable sites, see `--input_msi_threshold`).
* __PCGR MSI classifier__ (`--estimate_msi`): a statistical classifier built from 
somatic mutation profiles that separates _MSI.H_ (MSI-high) from _MSS_ (MS stable) 
tumors. It was trained using __N = 575__ exome-sequenced TCGA tumor samples with known 
MSI status (i.e. assayed from mononucleotide markers), and obtained a 
[positive predictive value]({wikipedia_url}Positive_and_negative_predictive_values#Positive_predictive_value) 
of 97.9% and a [negative predictive value]({wikipedia_url}Positive_and_negative_predictive_values#Negative_predictive_value) 
of 99.3% on an independent test set of __N = 245__ samples (details 
<a href="https://sigven-pcgr-msi-classifier.share.connect.posit.cloud/" 
target="_blank">here</a>). Given the nature of the training dataset, it can only 
be applied for samples originating from WGS/WES tumor-control 
sequencing assays (i.e. _not_ for tumor-only settings). 

An external MSI-H or MSS call takes precedence over the classifier for the 
assessment of MSI-H biomarker evidence, and any discordance between the two is flagged. 
Otherwise, biomarker evidence is shown only for high-confidence classifier predictions of MSI-H. 
