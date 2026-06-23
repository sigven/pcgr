Microsatellite instability (MSI) is the result of impaired DNA mismatch repair and 
constitutes a cellular phenotype of clinical significance in many cancer types, 
most prominently colorectal cancers, stomach cancers, endometrial cancers, 
and ovarian cancers ([Cortes-Ciriano et al., 2017]({pubmed_url}28585546)). 
We have built a statistical MSI classifier from somatic mutation profiles that 
separates _MSI.H_ (MSI-high) from _MSS_ (MS stable) tumors. The MSI classifier 
was trained using __N = 575__ exome-sequenced TCGA tumor samples with known 
MSI status (i.e. assayed from mononucleotide markers), and obtained a 
[positive predictive value]({wikipedia_url}Positive_and_negative_predictive_values#Positive_predictive_value) 
of 97.9% and a [negative predictive value]({wikipedia_url}Positive_and_negative_predictive_values#Negative_predictive_value) 
of 99.3% on an independent test set of __N = 245__ samples. Details of the MSI 
classification approach can be found <a href="https://rpubs.com/sigven/msi_classifier_2025" 
target="_blank">here</a>.

Note that, given the nature of the training dataset, the MSI classifier can only 
be applied for samples originating from WGS/WES tumor-control 
sequencing assays (i.e. _not_ for tumor-only settings). 
