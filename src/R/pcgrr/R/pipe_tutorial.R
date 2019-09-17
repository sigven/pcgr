# test_fun <- function(data, pi = T){
#   stopifnot(is.data.frame(data) & "VAR_ID" %in% colnames(data))
#   data %>% dplyr::mutate(GENOMIC_CHANGE = VAR_ID)
# }
#
# df <- data.frame('VAR_ID' = 234, 'CHROM' = '24', 'POS' = 1, stringsAsFactors = F)
# df <- rbind(df, data.frame('VAR_ID' = 234, 'CHROM' = '51', 'POS' = 2, stringsAsFactors = F))
# df <- rbind(df, data.frame('VAR_ID' = 234, 'CHROM' = '44', 'POS' = 2, stringsAsFactors = F))
# df <- rbind(df, data.frame('VAR_ID' = 234, 'CHROM' = '53', 'POS' = 2, stringsAsFactors = F))
# df <- rbind(df, data.frame('VAR_ID' = 234, 'CHROM' = '33', 'POS' = 1, stringsAsFactors = F))
# df <- rbind(df, data.frame('VAR_ID' = 234, 'CHROM' = 'Xy', 'POS' = 1, stringsAsFactors = F))
#
# df %>% test_fun()
#
#
# library(BSgenome.Hsapiens.UCSC.hg19)
# source('/Users/sigven/research/tcga/data-raw/data_cleaning_utils.R')
#
# genome_seq <- BSgenome.Hsapiens.UCSC.hg19
# seqinfo <- GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), genome = 'hg19')
#
# tmp_pathogenic <- xlsx::read.xlsx("/Users/sigven/research/docker/cpsr/tcga_pancan_pathogenic_germline.xlsx",1,startRow = 1,endRow = 854, header = T)
# tmp_pathogenic <- tmp_pathogenic %>% dplyr::select(HUGO_Symbol, HGVSp_short, Variant_Classification, Feature, Positive_Evidence, Negative_Evidence, Chromosome, Reference, Alternate, VARIANT_CLASS, Start, Stop, CharGer_Classification, ClinVar_Pathogenicity)
# tmp_pri_VUS <- xlsx::read.xlsx("/Users/sigven/research/docker/cpsr/tcga_pancan_pathogenic_germline.xlsx",2,startRow = 1,endRow = 541, header = T)
# tmp_pri_VUS <- tmp_pri_VUS %>% dplyr::select(HUGO_Symbol, HGVSp_short, Variant_Classification, Feature, Positive_Evidence, Negative_Evidence, Chromosome, Reference, Alternate, VARIANT_CLASS, Start, Stop, CharGer_Classification, ClinVar_Pathogenicity)

# tcga_variants <- dplyr::bind_rows(tmp_pathogenic, tmp_pri_VUS) %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(CharGer_Classification = stringr::str_replace_all(CharGer_Classification," ","_")) %>%
#   dplyr::mutate(ClinVar_Pathogenicity = stringr::str_replace_all(ClinVar_Pathogenicity," ","_")) %>%
#   dplyr::rename(Reference_Allele = Reference, Tumor_Seq_Allele2 = Alternate, Start_Position = Start, End_Position = Stop, Variant_Type = VARIANT_CLASS) %>%
#   dplyr::rename(CharGer_Consequence = Variant_Classification, CharGer_HGVSp = HGVSp_short, CharGer_Positive_Evidence = Positive_Evidence, CharGer_Negative_Evidence = Negative_Evidence, CharGer_Feature = Feature, CharGer_ClinVar = ClinVar_Pathogenicity) %>%
#   dplyr::mutate(Start_Position = as.integer(Start_Position)) %>%
#   dplyr::mutate(CharGer_HGVSp = stringr::str_replace_all(CharGer_HGVSp,"\\*","X")) %>%
#   dplyr::mutate(End_Position = as.integer(End_Position)) %>%
#   dplyr::mutate(Chromosome = paste0('chr',Chromosome))
#
# tcga_variants2 <- get_proper_maf_alleles(tcga_variants, genome_seq = genome_seq, seqinfo = seqinfo)
# tcga_variants2$QUAL <- '.'
# tcga_variants2$FILTER <- 'PASS'
# tcga_variants2$ID <- '.'
# tcga_variants2$INFO <- paste0('CharGer_Consequence=',tcga_variants2$CharGer_Consequence,';CharGer_HGVSp=',tcga_variants2$CharGer_HGVSp,';CharGer_Symbol=',tcga_variants2$HUGO_Symbol,';CharGer_Feature=',tcga_variants2$CharGer_Feature,';CharGer_Negative_Evidence=',tcga_variants2$CharGer_Negative_Evidence,';CharGer_Positive_Evidence=',tcga_variants2$CharGer_Positive_Evidence,';CharGer_Classification=',tcga_variants2$CharGer_Classification,';CharGer_ClinVar=', tcga_variants2$CharGer_ClinVar)
#
# tcga_variants_vcf <- tcga_variants2 %>%
#   pcgrr::order_variants(chrom_var = "CHROM", pos_var = "POS") %>%
#   dplyr::filter(!(POS == 32972574 & CharGer_Feature == 'ENST00000380152')) %>%
#   dplyr::filter(!(POS == 89831233 & CharGer_Feature == 'ENST00000305699')) %>%
#   dplyr::select(-c(HUGO_Symbol,CharGer_HGVSp, CharGer_Consequence, CharGer_Feature, CharGer_Positive_Evidence, CharGer_Negative_Evidence, Chromosome,Reference_Allele,Tumor_Seq_Allele2,Variant_Type,Start_Position,End_Position,CharGer_Classification,CharGer_ClinVar,GENOMIC_CHANGE)) %>%
#   dplyr::select(CHROM,POS,ID,REF,ALT,QUAL,FILTER, INFO) %>%
#   dplyr::distinct()
#
# header_lines <- c("##fileformat=VCFv4.2","##INFO=<ID=CharGer_Consequence,Number=.,Type=String,Description=\"CharGer Consequence\">","##INFO=<ID=CharGer_Feature,Number=.,Type=String,Description=\"CharGer Feature\">","##INFO=<ID=CharGer_HGVSp,Number=.,Type=String,Description=\"CharGer HGVSp\">","##INFO=<ID=CharGer_Negative_Evidence,Number=.,Type=String,Description=\"CharGer Negative Evidence\">","##INFO=<ID=CharGer_Positive_Evidence,Number=.,Type=String,Description=\"CharGer Positive Evidence\">","##INFO=<ID=CharGer_Symbol,Number=.,Type=String,Description=\"CharGer Symbol\">","##INFO=<ID=CharGer_Classification,Number=.,Type=String,Description=\"CharGer Classification\">","##INFO=<ID=CharGer_ClinVar,Number=.,Type=String,Description=\"CharGer ClinVar Classification\">","#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
# write.table(header_lines,file="tcga_pancancer_germline_pathogenic_vus.grch37.vcf",sep="\n",row.names = F,col.names = F,quote=F)
# write.table(tcga_variants_vcf,file="tcga_pancancer_germline_pathogenic_vus.grch37.vcf",quote=F,sep="\t",col.names = F,row.names = F,append = T)
#
