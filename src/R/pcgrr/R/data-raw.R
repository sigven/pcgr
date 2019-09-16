acmg_evidence_codes <- read.table(file="data-raw/acmg_evidence.tsv",header=T,stringsAsFactors = F,comment.char="",na.strings=c("NA"),sep="\t")

significance_colors <- c("#9E0142", "#D53E4F", "#000000", "#78C679", "#077009")
significance_levels <- c('Pathogenic','Likely_Pathogenic','VUS','Likely_Benign','Benign')


usethis::use_data(acmg_evidence_codes, overwrite = T)
usethis::use_data(significance_colors, overwrite = T)
usethis::use_data(significance_levels, overwrite = T)

heterozygous_states <- c()
ref_allele_index <- 0
while(ref_allele_index < 20){
  alt_allele_index <- ref_allele_index + 1
  while(alt_allele_index <= 20){
    phased_gt_1 <- paste0(ref_allele_index,'|',alt_allele_index)
    phased_gt_2 <- paste0(alt_allele_index,'|',ref_allele_index)
    unphased_gt_1 <- paste0(ref_allele_index,'/',alt_allele_index)
    unphased_gt_2 <- paste0(alt_allele_index,'/',ref_allele_index)
    heterozygous_states <- c(heterozygous_states,phased_gt_1,phased_gt_2,unphased_gt_1,unphased_gt_2)
    alt_allele_index <- alt_allele_index + 1
  }
  ref_allele_index <- ref_allele_index + 1
}
homozygous_states <- c()
hom_allele_index <- 1
while(hom_allele_index <= 10){
  phased_gt <- paste0(hom_allele_index,'|',hom_allele_index)
  unphased_gt <- paste0(hom_allele_index,'/',hom_allele_index)
  homozygous_states <- c(homozygous_states,phased_gt,unphased_gt)
  hom_allele_index <- hom_allele_index + 1
}

usethis::use_data(heterozygous_states,overwrite = T)
usethis::use_data(homozygous_states,overwrite = T)

variant_db_url <- data.frame('name' = 'DBSNP',group_by_var = 'VAR_ID',url_prefix = 'http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',
                       link_key_var = "DBSNPRSID",link_display_var = "DBSNPRSID", stringsAsFactors = F)
variant_db_url <- dplyr::bind_rows(variant_db_url,
                             data.frame('name' = 'CLINVAR',group_by_var = 'VAR_ID',url_prefix = 'http://www.ncbi.nlm.nih.gov/clinvar/variation/',
                                        link_key_var = "CLINVAR_MSID",link_display_var = "CLINVAR_TRAITS_ALL", stringsAsFactors = F),
                             data.frame('name' = 'GENE_NAME',group_by_var = 'VAR_ID',url_prefix = 'http://www.ncbi.nlm.nih.gov/gene/',
                                        link_key_var = "ENTREZ_ID",link_display_var = "GENENAME", stringsAsFactors = F),
                             data.frame('name' = 'PROTEIN_DOMAIN',group_by_var = 'VAR_ID',url_prefix = 'http://pfam.xfam.org/family/',
                                        link_key_var = "PFAM_DOMAIN",link_display_var = "PFAM_DOMAIN_NAME", stringsAsFactors = F),
                             data.frame('name' = 'COSMIC',group_by_var = 'VAR_ID',url_prefix = 'http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=',
                                        link_key_var = "COSMIC_MUTATION_ID",link_display_var = "COSMIC_MUTATION_ID", stringsAsFactors = F)
)

usethis::use_data(variant_db_url,overwrite = T)

effect_prediction_algos <- read.table(file="data-raw/effect_prediction_algorithms.tsv",header=T,sep="\t",quote="",stringsAsFactors = F)
usethis::use_data(effect_prediction_algos,overwrite = T)
