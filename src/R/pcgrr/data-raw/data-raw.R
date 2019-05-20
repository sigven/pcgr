acmg_evidence_codes <- read.table(file="data-raw/acmg_evidence.tsv",header=T,stringsAsFactors = F,comment.char="",na.strings=c("NA"),sep="\t")

significance_colors <- c("#9E0142", "#D53E4F", "#000000", "#78C679", "#077009")
significance_levels <- c('Pathogenic','Likely_Pathogenic','VUS','Likely_Benign','Benign')


usethis::use_data(acmg_evidence_codes, overwrite = T)
usethis::use_data(significance_colors, overwrite = T)
usethis::use_data(significance_levels, overwrite = T)
