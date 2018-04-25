

# suppressWarnings(suppressPackageStartupMessages(library(pcgrr)))
# suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
# suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19)))
# suppressWarnings(suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38)))
# suppressWarnings(suppressPackageStartupMessages(library(deconstructSigs)))
# suppressWarnings(suppressPackageStartupMessages(library(randomForest)))
# suppressWarnings(suppressPackageStartupMessages(library(caret)))
# suppressWarnings(suppressPackageStartupMessages(library(RcppTOML)))
# # # #
# # # project_directory <- '/Users/sigven/research/docker/pcgr_predispose'
#pcgr_config_file <- '/Users/sigven/research/docker/pcgr/examples/pcgr_conf.thymic_penile.toml'
#pcgr_config <-  RcppTOML::parseTOML(pcgr_config_file, fromFile = T)
# # # sample_name <- 'GCF0167-0001-N01'
# # # pcgr_version <- '0.1.0'
# # # query_vcf2tsv <- '/Users/sigven/research/docker/pcgr_predispose/GCF0167-0001-N01.pcgr_predispose.pass.tsv.gz'
# # # tsv_gz_file <- query_vcf2tsv
# # # load(file='/Users/sigven/research/docker/pcgr/data/grch37/rda/pcgr_data.rda')
# # # genome_assembly <- 'grch37'
# # # biomarker_mapping_stringency <- 1
# # # callset <- 'somatic_calls'
# # #
# # #
# # ### TEST - TIER MODEL PCGR - COAD
#pcgr_config_file <- '/Users/sigven/research/docker/pcgr_predispose/pcgr_predispose.toml'
#pcgr_config <-  RcppTOML::parseTOML(pcgr_config_file, fromFile = T)
# sample_name <- 'cov50-ensembl-annotated-decomposed'
# pcgr_version <- '0.6.0'
#query_vcf2tsv <- '/Users/sigven/research/tumor_seq_projects/mesothelioma_germline/GCF0167-0001-N01.pcgr_predispose.pass.tsv.gz'
# tsv_gz_file <- query_vcf2tsv
# # # query_cnv <- '/Users/sigven/research/tumor_seq_projects/mesothelioma_somatic/cnv/facets/GCF0167-0001-N01_GCF0167-0001-T01.cna.tsv'
# # # cna_file <- query_cnv
# # # # # # # # cna_segments_tsv <- query_cnv
# # # # # load(file='/Users/sigven/research/docker/pcgr/data/grch37/rda/pcgr_data.rda')
# # # # genome_assembly <- 'grch38'
# genome_seq <- BSgenome.Hsapiens.UCSC.hg19
# biomarker_mapping_stringency <- 1
# callset <- 'somatic_calls'
# # # # # #
# pcgr_config_file <- '/Users/sigven/research/docker/pcgr/examples/pcgr_conf.no_validation.toml'
# pcgr_config <-  RcppTOML::parseTOML(pcgr_config_file, fromFile = T)
# # # # # # query_cnv <- '/Users/sigven/research/tumor_seq_projects/mesothelioma_somatic/cnv/facets/GCF0167-0001-N01_GCF0167-0001-T01.cna.tsv'
# # # # # # cna_file <- query_cnv
# # # # # # cna_segments_tsv <- query_cnv
# load(file='/Users/sigven/research/docker/pcgr/data/grch37/rda/pcgr_data.rda')
# genome_assembly <- 'grch37'
# #genome_seq <- BSgenome.Hsapiens.UCSC.hg19
# # # pcgr_version <- '0.6.0'
# # transcript_overlap_pct <- pcgr_config$cna$cna_overlap_pct
# #
# #
# #
#
#
#
#
#
#

#' A function that converts the INFO tags in a VRanges object into a basic types
#'
#' @param vr VRanges object
#' @return vr Vranges object with INFO tags more simply formatted
#'

# postprocess_vranges_info <- function(vr){
#   ## Convert IntegerLists and CharacterLists to basic character lists
#   vcf_annotations_df <- NULL
#   for(tag in colnames(GenomicRanges::mcols(vr))){
#     mcol_class <- class(GenomicRanges::mcols(vr)[,c(tag)])[1]
#
#     if(mcol_class != "character" & mcol_class != "integer" & mcol_class != "logical" & mcol_class != "numeric"){
#       annotation_track <- NULL
#       #cat("TAG: ",tag, ', type:',mcol_class,'\n')
#       if(mcol_class == "CompressedCharacterList"){
#         annotation_track <- data.frame(val = as.character(Biostrings::unstrsplit(GenomicRanges::mcols(vr)[[tag]], sep=',')))
#       }else{
#         annotation_track <- data.frame(val = as.character(sapply(GenomicRanges::mcols(vr)[[tag]], paste, collapse=",")))
#       }
#       if(is.null(vcf_annotations_df)){
#         vcf_annotations_df <- data.frame(annotation_track$val)
#         names(vcf_annotations_df) <- c(tag)
#         vcf_annotations_df[,c(tag)] <- as.character(vcf_annotations_df[,c(tag)])
#       }
#       else{
#         vcf_annotations_df[,c(tag)] <- as.character(annotation_track$val)
#       }
#       ## add NA to empty values
#       if(nrow(as.data.frame(vcf_annotations_df[nchar(vcf_annotations_df[,c(tag)]) == 0,])) != 0){
#         if(dim(vcf_annotations_df)[2] == 1){
#           vcf_annotations_df[nchar(vcf_annotations_df[,c(tag)]) == 0,] <- NA
#         }
#         else{
#           vcf_annotations_df[nchar(vcf_annotations_df[,c(tag)]) == 0,][,c(tag)] <- NA
#         }
#       }
#     }
#     else{
#       #cat("TAG: ",tag, ', type:',mcol_class,'\n')
#       if(is.null(vcf_annotations_df)){
#         vcf_annotations_df <- data.frame(GenomicRanges::mcols(vr)[,c(tag)])
#         names(vcf_annotations_df) <- c(tag)
#       }
#       else{
#         vcf_annotations_df[,c(tag)] <- GenomicRanges::mcols(vr)[,c(tag)]
#       }
#     }
#   }
#
#   ## add variant_id and sample_id
#   position <- GenomicRanges::start(GenomicRanges::ranges(vr))
#   vcf_annotations_df['VAR_ID'] <- paste(as.character(GenomeInfoDb::seqnames(vr)),position,VariantAnnotation::ref(vr),VariantAnnotation::alt(vr),sep="_")
#   #vcf_annotations_df['VCF_SAMPLE_ID'] <- as.character(VariantAnnotation::sampleNames(vr))
#   GenomicRanges::mcols(vr) <- S4Vectors::DataFrame(vcf_annotations_df)
#
#   return(vr)
#
# }

#' Function that adds PFAM name descriptions to PFAM identifiers
#'
#' @param vcf_data_df
#'
#' @return vcf_data_df_pfam
#'
# add_pfam_domain_links <- function(vcf_data_df, pcgr_data){
#
#   rlogging::message("Extending annotation descriptions related to PFAM protein domains")
#   if("DOMAINS" %in% colnames(vcf_data_df) & "VAR_ID" %in% colnames(vcf_data_df)){
#     pfam_df <- dplyr::select(vcf_data_df,DOMAINS,VAR_ID) %>% dplyr::filter(!is.na(DOMAINS))
#     if(nrow(pfam_df) == 0){
#       vcf_data_df$PROTEIN_DOMAIN <- NA
#       return(vcf_data_df)
#     }
#     pfam_df <- pfam_df %>% dplyr::distinct() %>% tidyr::separate_rows(DOMAINS,sep="&") %>% dplyr::filter(stringr::str_detect(DOMAINS,"Pfam_domain"))
#     pfam_df$DOMAINS <- stringr::str_replace(pfam_df$DOMAINS,"Pfam_domain:","")
#     pfam_df <- dplyr::left_join(pfam_df,pcgr_data$pfam_domains,by=c("DOMAINS" = "pfam_id")) %>% dplyr::select(VAR_ID,url)
#     pfam_df <- dplyr::rename(pfam_df, PD = url)
#     pfam_ret <- as.data.frame(dplyr::group_by(pfam_df, VAR_ID) %>% dplyr::summarise(PROTEIN_DOMAIN = paste(PD, collapse=", ")))
#
#     if(nrow(pfam_ret) > 0){
#       vcf_data_df <- dplyr::left_join(vcf_data_df,pfam_ret,by=c("VAR_ID" = "VAR_ID"))
#     }
#     else{
#       vcf_data_df$PROTEIN_DOMAIN <- NA
#     }
#   }
#   else{
#     vcf_data_df$PROTEIN_DOMAIN <- NA
#   }
#
#   return(vcf_data_df)
# }



# df <- data.frame(
#   x = rep(seq(2, 15, 6.5), 2),
#   y = c(rep(2,3), rep(6.5, 3)),
#   h = rep(4, 6),
#   w = rep(6, 6),
#   info = c("78%\nmeaningless plots",
#            "+10K\nhours wasted",
#            "8/10\nzombies prefer brains",
#            "ALL\ndogs go to heaven",
#            "6\ninfoboxes",
#            "< 0.5\ntarget pvalue"),
#   color = factor(1:6)
# )
#
#
# ggplot(df, aes(x, y, height = h, width = w, label = info, fill = color)) +
#   geom_tile() +
#   geom_text(color = "white", fontface = "bold") +
#   coord_fixed() +
#   scale_fill_brewer(type = "qual",palette = "Dark2") +
#   theme_void() +
#   guides(fill = F)


# if(!file.exists(vcf_gz_file) | file.size(vcf_gz_file) == 0){
#   rlogging::stop(paste0("File ",vcf_gz_file," does not exist or has zero size"))
# }
# rlogging::message(paste0("Reading and parsing VEP/vcfanno-annotated VCF file - ",vcf_gz_file))
# vcf_data_vr <- VariantAnnotation::readVcfAsVRanges(vcf_gz_file,genome = genome_assembly)
# n_all_unfiltered_calls <- length(vcf_data_vr)
# #vcf_data_vr <- vcf_data_vr[!is.na(vcf_data_vr$GT) & !(vcf_data_vr$GT == '.'),]
# vcf_data_vr <- vcf_data_vr[VariantAnnotation::called(vcf_data_vr)]
# vcf_data_vr <- pcgrr::postprocess_vranges_info(vcf_data_vr)
# vcf_data_df <- as.data.frame(vcf_data_vr)
# if(!is.null(sample_name) & nrow(vcf_data_df) > 0){
#   vcf_data_df$VCF_SAMPLE_ID <- sample_name
# }
# rlogging::message(paste0("Number of PASS variants: ",nrow(vcf_data_df)))
# if(any(grepl(paste0("VARIANT_CLASS$"),names(vcf_data_df)))){
#   n_snvs <- nrow(vcf_data_df[!is.na(vcf_data_df$VARIANT_CLASS) & vcf_data_df$VARIANT_CLASS == 'SNV',])
#   n_deletions <- nrow(vcf_data_df[!is.na(vcf_data_df$VARIANT_CLASS) & vcf_data_df$VARIANT_CLASS == 'deletion',])
#   n_insertions <- nrow(vcf_data_df[!is.na(vcf_data_df$VARIANT_CLASS) & vcf_data_df$VARIANT_CLASS == 'insertion',])
#   n_substitutions <- nrow(vcf_data_df[!is.na(vcf_data_df$VARIANT_CLASS) & vcf_data_df$VARIANT_CLASS == 'substitution',])
#   rlogging::message(paste0("Number of SNVs: ",n_snvs))
#   rlogging::message(paste0("Number of deletions: ",n_deletions))
#   rlogging::message(paste0("Number of insertions: ",n_insertions))
#   rlogging::message(paste0("Number of block substitutions: ",n_substitutions))
# }
# if(nrow(vcf_data_df) == 0){
#   rlogging::warning("Number of PASS variants in input VCF is 0 - no variants will be written to HTML report")
#   for(col in c('GENOME_VERSION','GENOMIC_CHANGE','PROTEIN_DOMAIN','PROTEIN_FEATURE','PROTEIN_CHANGE',
#                'DOCM_LITERATURE','DOCM_DISEASE','GENENAME','GENE_NAME','CLINVAR_TRAITS_ALL','CLINVAR',
#                'TCGA_FREQUENCY','COSMIC','DBSNP','KEGG_PATHWAY','TARGETED_DRUGS','CANCER_ASSOCIATIONS',
#                'CHEMBL_COMPOUND_TERMS','DISGENET_TERMS','CALL_CONFIDENCE')){
#     vcf_data_df[col] <- character(nrow(vcf_data_df))
#   }
#   for(col in c('DP_TUMOR','DP_NORMAL')){
#     vcf_data_df[col] <- integer(nrow(vcf_data_df))
#   }
#   for (col in c('AF_TUMOR','AF_NORMAL')){
#     vcf_data_df[col] <- numeric(nrow(vcf_data_df))
#   }
#   vcf_data_df <- dplyr::rename(vcf_data_df, CHROM = seqnames, POS = start, REF = ref, ALT = alt, CONSEQUENCE = Consequence)
#   return(vcf_data_df)
# }
#
#
#
# vcf_data_df <- dplyr::mutate(vcf_data_df, GENOME_VERSION = genome_assembly, PROTEIN_CHANGE = HGVSp)
# vcf_data_df <- dplyr::rename(vcf_data_df, CHROM = seqnames, POS = start, REF = ref, ALT = alt, CONSEQUENCE = Consequence)
# if(nrow(vcf_data_df[!is.na(vcf_data_df$PROTEIN_CHANGE) & stringr::str_detect(vcf_data_df$PROTEIN_CHANGE,":"),]) > 0){
#   vcf_data_df[!is.na(vcf_data_df$PROTEIN_CHANGE) & stringr::str_detect(vcf_data_df$PROTEIN_CHANGE,":"),]$PROTEIN_CHANGE <- stringr::str_split_fixed(vcf_data_df[!is.na(vcf_data_df$PROTEIN_CHANGE) & stringr::str_detect(vcf_data_df$PROTEIN_CHANGE,":"),]$PROTEIN_CHANGE,pattern = ":",2)[,2]
# }
# if(nrow(vcf_data_df[!is.na(vcf_data_df$PROTEIN_CHANGE) & stringr::str_detect(vcf_data_df$PROTEIN_CHANGE,"^ENSP"),]) > 0){
#   vcf_data_df[!is.na(vcf_data_df$PROTEIN_CHANGE) & stringr::str_detect(vcf_data_df$PROTEIN_CHANGE,"^ENSP"),]$PROTEIN_CHANGE <- NA
# }
# vcf_data_df$CODING_STATUS <- 'noncoding'
# coding_consequence_pattern <- "^(stop_|start_lost|frameshift_|missense_variant|splice_donor|splice_acceptor|inframe_)"
# if(nrow(vcf_data_df[!is.na(vcf_data_df$CONSEQUENCE) & stringr::str_detect(vcf_data_df$CONSEQUENCE,coding_consequence_pattern),]) > 0){
#   vcf_data_df[!is.na(vcf_data_df$CONSEQUENCE) & stringr::str_detect(vcf_data_df$CONSEQUENCE,coding_consequence_pattern),]$CODING_STATUS <- 'coding'
# }
#
# for(col in c('width','strand','totalDepth','refDepth','altDepth','sampleNames')){
#   if(col %in% colnames(vcf_data_df)){
#     vcf_data_df[,col] <- NULL
#   }
# }
# vcf_data_df$GENOMIC_CHANGE <- paste0(vcf_data_df$CHROM,":g.",vcf_data_df$POS,vcf_data_df$REF,">",vcf_data_df$ALT)
#
# if('EXON' %in% colnames(vcf_data_df)){
#   vcf_data_df$EXON <- as.integer(stringr::str_split_fixed(vcf_data_df$EXON,"/",2)[,1])
# }
# vcf_data_df <- pcgrr::add_pfam_domain_links(vcf_data_df, pcgr_data = pcgr_data)
# vcf_data_df <- pcgrr::add_swissprot_feature_descriptions(vcf_data_df, pcgr_data = pcgr_data)
# vcf_data_df <- pcgrr::add_read_support(vcf_data_df, pcgr_config)
#
