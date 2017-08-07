library(BSgenome.Hsapiens.UCSC.hg19)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(magrittr)
library(deconstructSigs)
library(data.table)
library(ggplot2)
library(plotly)
library(rlogging)


#' Function that plots a histogram of the the variant allelic support (tumor) - grouped by tiers
#'
#' @param tier_df data frame with somatic mutations
#' @param bin_size size of bins for allelic frequency
#' @return p geom_histogram plot from ggplot2
#'
#'
tier_af_distribution <- function(tier_df, bin_size = 0.1){
  af_bin_df <- data.frame()
  i <- 1
  num_bins <- as.integer(1 / bin_size)
  bin_start <- 0
  while(i <= num_bins){
    bin_end <- bin_start + bin_size
    bin_name <- as.character(paste0(bin_start,' - ',bin_end))
    j <- 1
    while(j <= 5){
      TIER <- paste0('TIER ',j)
      df <- data.frame(bin_name = bin_name, bin_start = bin_start, bin_end = bin_end, bin = as.integer(i), TIER = TIER, stringsAsFactors = F)
      af_bin_df <- rbind(af_bin_df, df)
      j <- j + 1
    }
    bin_start <- bin_end
    i <- i + 1
  }

  tier_df_trans <- transform(tier_df, bin = cut(AF_TUMOR, breaks = seq(from = 0, to = 1, by = bin_size), right = F, include.lowest = T, labels = F))
  tier_df_trans_bin <- as.data.frame(dplyr::group_by(tier_df_trans, TIER, bin) %>% dplyr::summarise(Count = n()))

  af_bin_df <- dplyr::left_join(af_bin_df, tier_df_trans_bin, by = c("bin","TIER"))

  p <- ggplot2::ggplot(data = af_bin_df) + ggplot2::geom_bar(mapping = ggplot2::aes(x = bin_name, y = Count, fill = TIER), stat = "identity") +
    ggplot2::scale_color_brewer(palette='Dark2') +
    ggplot2::theme_classic() +
    ggplot2::ylab("Count") +
    ggplot2::xlab("Variant allelic fraction - tumor") +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   #axis.text.x = element_blank(),
                   axis.text.x=ggplot2::element_text(angle=45,family="Helvetica",size=12,vjust = -0.1),
                   axis.title.x=ggplot2::element_text(family="Helvetica",size=12,vjust = -2),
                   axis.text.y=ggplot2::element_text(family="Helvetica",size=12),
                   axis.title.y=ggplot2::element_text(family="Helvetica",size=12,vjust=1.5),
                   plot.margin = (grid::unit(c(0.5, 2, 2, 0.5), "cm")),
                   legend.text=ggplot2::element_text(family="Helvetica",size=12))

  return(p)

}

#' Function that check for valid mutations/chromosomes in input
#'
#' @param vcf_data_df data frame
#' @param chromosome_column name of columns for which string replace is to be performed
#' @param pattern pattern to replace
#' @return vcf_data_df_valid data frame with valid mutations
#'
get_valid_chromosomes <- function(vcf_data_df, chromosome_column = 'CHROM', bsg = BSgenome.Hsapiens.UCSC.hg19){
  vcf_data_df_valid <- vcf_data_df
  vcf_data_df_valid[, chromosome_column] <- factor(vcf_data_df_valid[, chromosome_column])
  levels(vcf_data_df_valid[, chromosome_column]) <- sub("^([0-9XY])", "chr\\1", levels(vcf_data_df_valid[, chromosome_column]))
  levels(vcf_data_df_valid[, chromosome_column]) <- sub("^MT", "chrM", levels(vcf_data_df_valid[, chromosome_column]))
  levels(vcf_data_df_valid[, chromosome_column]) <- sub("^(GL[0-9]+).[0-9]", "chrUn_\\L\\1", levels(vcf_data_df_valid[, chromosome_column]), perl = T)
  unknown.regions <- levels(vcf_data_df_valid[, chromosome_column])[which(!(levels(vcf_data_df_valid[, chromosome_column]) %in% GenomeInfoDb::seqnames(bsg)))]
  if (length(unknown.regions) > 0) {
    unknown.regions <- paste(unknown.regions, collapse = ',\ ')
    rlogging::warning(paste('Check chr names -- not all match BSgenome.Hsapiens.UCSC.hg19::Hsapiens object:\n', unknown.regions, sep = ' '))
    vcf_data_df_valid <- vcf_data_df_valid[vcf_data_df_valid[, chromosome_column] %in% GenomeInfoDb::seqnames(bsg), ]
  }
  return(vcf_data_df_valid)

}

#' Function that performs stringr::str_replace on strings of multiple string columns of a dataframe
#'
#' @param df data frame
#' @param strings name of columns for which string replace is to be performed
#' @param pattern pattern to replace
#' @param replacement string to replace
#' @param replace_all logical - replace all occurrences
#' @return df
#'
#'
df_string_replace <- function(df, strings, pattern, replacement, replace_all = F){

  for(column_name in strings){
    if(column_name %in% colnames(df)){
      if(replace_all == F){
        df[,column_name] <- stringr::str_replace(df[,column_name],pattern = pattern, replacement = replacement)
      }else{
        df[,column_name] <- stringr::str_replace_all(df[,column_name],pattern = pattern, replacement = replacement)
      }
    }
  }
  return(df)
}

#' Function that transforms a tier-structured variant data frame into a MAF-like data frame (for input to 2020plus, MutSigCV)
#'
#' @param tier_df data frame with somatic mutations
#' @return maf_df
#'
#'

tier_to_maf <- function(tier_df){
  maf_df <- dplyr::select(tier_df, SYMBOL, GENOMIC_CHANGE, VCF_SAMPLE_ID, CONSEQUENCE, VARIANT_CLASS)
  maf_df$Hugo_Symbol <- maf_df$SYMBOL
  locus_info <- dplyr::select(tidyr::separate(maf_df,GENOMIC_CHANGE, c('chrom','pos','alleles'),sep=":",convert=T), chrom, pos, alleles)
  locus_info <- tidyr::separate(locus_info,alleles, c('ref','alt'),sep=">",convert=T)
  locus_info$chrom <- stringr::str_replace(locus_info$chrom, pattern = "g\\.chr", replacement = '')
  maf_df$Chromosome <- locus_info$chrom
  maf_df$Reference_Allele <- locus_info$ref
  maf_df$Tumor_Seq_Allele2 <- locus_info$alt
  maf_df$Start_Position <- locus_info$pos
  maf_df$End_Position <- maf_df$Start_Position + nchar(maf_df$Reference_Allele) - 1
  maf_df$Tumor_Sample_Barcode <- maf_df$VCF_SAMPLE_ID
  maf_df$NCBI_Build <- 'GRCh37'

  maf_df$Variant_Classification <- character(nrow(maf_df))
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant)"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant)"),]$Variant_Classification <- 'Splice_Site'
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"stop_gained"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"stop_gained"),]$Variant_Classification <- 'Nonsense_Mutation'
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"frameshift_variant") & maf_df$VARIANT_CLASS == 'deletion',]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"frameshift_variant") & maf_df$VARIANT_CLASS == 'deletion',]$Variant_Classification <- 'Frame_Shift_Del'
  }
  if(nrow( maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"frameshift_variant") & maf_df$VARIANT_CLASS == 'insertion',]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"frameshift_variant") & maf_df$VARIANT_CLASS == 'insertion',]$Variant_Classification <- 'Frame_Shift_Ins'
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"stop_lost"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"stop_lost"),]$Variant_Classification <- 'Nonstop_Mutation'
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"initiator_codon_variant|start_lost"),])){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"initiator_codon_variant|start_lost"),]$Variant_Classification <- 'Translation_Start_Site'
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"inframe_insertion") & maf_df$VARIANT_CLASS == 'insertion',]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"inframe_insertion") & maf_df$VARIANT_CLASS == 'insertion',]$Variant_Classification <- 'In_Frame_Ins'
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"inframe_deletion") & maf_df$VARIANT_CLASS == 'deletion',]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"inframe_deletion") & maf_df$VARIANT_CLASS == 'deletion',]$Variant_Classification <- 'In_Frame_Del'
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(missense_variant|coding_sequence_variant|conservative_missense_variant|rare_amino_acid_variant)"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(missense_variant|coding_sequence_variant|conservative_missense_variant|rare_amino_acid_variant)"),]$Variant_Classification <- 'Missense_Mutation'
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"3_prime_UTR_variant"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"3_prime_UTR_variant"),]$Variant_Classification <- "3'UTR"
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"5_prime_UTR_variant"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"5_prime_UTR_variant"),]$Variant_Classification <- "5'UTR"
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region)"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region)"),]$Variant_Classification <- "IGR"
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant)"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant)"),]$Variant_Classification <- "Silent"
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant)"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant)"),]$Variant_Classification <- "RNA"
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(transcript_amplification|intron_variant|INTRAGENIC|intragenic_variant)"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"^(transcript_amplification|intron_variant|INTRAGENIC|intragenic_variant)"),]$Variant_Classification <- "Intron"
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"downstream_gene_variant"),]) > 0){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"downstream_gene_variant"),]$Variant_Classification <- "3'Flank"
  }
  if(nrow(maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"upstream_gene_variant"),])){
    maf_df[stringr::str_detect(maf_df$CONSEQUENCE,"upstream_gene_variant"),]$Variant_Classification <- "5'Flank"
  }

  maf_df <- dplyr::select(maf_df, Hugo_Symbol, Chromosome, NCBI_Build, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode, Variant_Classification)
  chrom_order <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  maf_df$Chromosome <- factor(maf_df$Chromosome, levels = chrom_order)
  maf_df <- dplyr::arrange(maf_df, Chromosome, Start_Position, End_Position)
  return(maf_df)
}


#' Function that generates cancer genome report
#'
#' @param project_directory name of project directory
#' @param query_vcf name of VCF file with annotated query SNVs/InDels
#' @param logR_gain logR treshold for annotating copy number amplifications
#' @param logR_homdel logR threshold for annotating homozygous deletions
#' @param pcgr_data List of data frames with PCGR data annotations
#' @param signature_normalization_method deconstructSigs normalization method
#' @param cna_segments_tsv name of CNA segments file (tab-separated values)
#' @param sample_name sample identifier
#' @param signatures_limit Limit mutational signature analysis to a given number of signatures
#' @param print_biomarkers Logical indicating if biomarker data is to be written to file
#' @param print_tier_variants Logical indicating if tiered variant data is to be written to file
#' @param print_mutational_signatures Logical indicating if mutational signature data is to be written to file
#' @param print_cna_segments Logical indicating if annotated cnv segment data is to be written to file
#' @param print_html_report Logical indicating if HTML report is to be printed
#' @param predict_msi Logical indicating if MSI is to be predicted
#' @param identify_msigs Logical indicating if inference of known mutational signatures will be performed
#' @param tumor_dp_tag Column name for sequencing depth at variant site in tumor
#' @param tumor_af_tag Column name for variant allelic fraction (VAF) in tumor
#' @param normal_dp_tag Column name for sequencing depth in normal/control sample
#' @param normal_af_tag Column name for variant allelic fraction (VAF) in normal/control sample
#' @param call_conf_tag Column name for somatic call confidence status
#'
#' @return p
#'

generate_pcg_report <- function(project_directory, query_vcf, logR_gain, logR_homdel, pcgr_data, sample_name = 'SampleX', signature_normalization_method = 'default', cna_segments_tsv = NULL, signatures_limit = 6, print_biomarkers = TRUE, print_tier_variants = TRUE, print_mutational_signatures = TRUE, print_cna_segments = TRUE, print_maf = TRUE, print_html_report = TRUE,show_noncoding = FALSE,predict_MSI = FALSE, identify_msigs = FALSE, tumor_dp_tag = '_na', tumor_af_tag = '_na', normal_dp_tag = '_na', normal_af_tag = '_na', call_conf_tag = '_na', pcgr_version = '0.4.1'){

  report_data <- list(tier1_report = FALSE, tier2_report = FALSE, tier3_report = FALSE, tier4_report = FALSE, tier5_report = FALSE, msi_report = FALSE, missing_msi_data = FALSE, signature_report = FALSE, cna_report_oncogene_gain = FALSE, cna_report_tsgene_loss = FALSE, cna_report_biomarkers = FALSE, cna_report_segments = FALSE, missing_signature_data = FALSE)

  tier_tsv_fname <- paste0(project_directory, '/',sample_name,'.pcgr.snvs_indels.tiers.tsv')
  msig_tsv_fname <- paste0(project_directory, '/',sample_name,'.pcgr.mutational_signatures.tsv')
  biomarker_tsv_fname <- paste0(project_directory, '/',sample_name,'.pcgr.snvs_indels.biomarkers.tsv')
  cna_tsv_fname <- paste0(project_directory, '/',sample_name,'.pcgr.cna_segments.tsv')
  maf_fname <- paste0(project_directory, '/',sample_name,'.pcgr.maf')

  if(query_vcf != 'None'){
    sample_calls <- pcgrr2::get_calls(query_vcf, pcgr_data, sample_id = sample_name, tumor_dp_tag = tumor_dp_tag, tumor_af_tag = tumor_af_tag, normal_dp_tag = normal_dp_tag, normal_af_tag = normal_af_tag, call_conf_tag = call_conf_tag)
    report_data <- pcgrr2::generate_report_data(sample_calls, pcgr_data = pcgr_data, sample_name = sample_name, minimum_n_signature_analysis = 50, signature_normalization_method = signature_normalization_method, signatures_limit = signatures_limit, predict_MSI = predict_MSI, identify_msigs = identify_msigs)
    report_data$sample_name <- sample_name

    if(!is.null(report_data$tsv_variants) & print_tier_variants == TRUE){
      write.table(report_data$tsv_variants,file=tier_tsv_fname, sep="\t",col.names = T,row.names = F,quote = F)
    }
    if(!is.null(report_data$tsv_biomarkers) & print_biomarkers == TRUE){
      write.table(report_data$tsv_biomarkers,file=biomarker_tsv_fname, sep="\t",col.names = T,row.names = F,quote = F)
    }

    if(!is.null(report_data$signature_data) & print_mutational_signatures == TRUE){
      if(length(report_data$signature_data) > 1){
        sample_mutational_signatures <- report_data$signature_data$signatures_cancertypes_aetiologies
        sample_mutational_signatures$SampleID <- sample_name
        sample_mutational_signatures <- dplyr::select(sample_mutational_signatures,-Comments)
        write.table(sample_mutational_signatures,file=msig_tsv_fname, sep="\t",col.names = T,row.names = F,quote = F)
      }
    }
    if(!is.null(report_data$maf_df) & print_maf == TRUE){
      write.table(report_data$maf_df,file=maf_fname, sep="\t",col.names = T,row.names = F,quote = F)
    }

  }

  cna_data <- list(ranked_segments = data.frame(), oncogene_amplified = data.frame(), tsgene_homozygous_deletion = data.frame(),cna_df_for_print = data.frame(), cna_biomarkers = data.frame(), cna_biomarker_segments = data.frame())

  if(!is.null(cna_segments_tsv)){
    if(file.exists(cna_segments_tsv)){
      report_data$cna_report_tsgene_loss <- TRUE
      report_data$cna_report_oncogene_gain <- TRUE
      report_data$cna_report_biomarkers <- TRUE
      report_data$cna_report_segments <- TRUE

      cna_data <- pcgrr2::cna_segment_annotation(cna_segments_tsv, logR_gain, logR_homdel, pcgr_data)
      if(nrow(cna_data$cna_df_for_print) > 0 & print_cna_segments == TRUE){
        write.table(cna_data$cna_df_for_print,file=cna_tsv_fname,col.names = T,row.names = F,quote=F,sep="\t")
        gzip_command <- paste0('gzip -f ',cna_tsv_fname)
        system(gzip_command, intern=F)
      }
    }
  }

  if(print_html_report == TRUE){

    eval_tier1 <- report_data$tier1_report
    eval_tier2 <- report_data$tier2_report
    eval_tier3 <- report_data$tier3_report
    eval_tier4 <- report_data$tier4_report
    eval_tier5 <- report_data$tier5_report
    eval_signature_report <- report_data$signature_report
    eval_missing_signature_data <- report_data$missing_signature_data
    eval_cna_segments <- report_data$cna_report_segments
    eval_cna_loss <- report_data$cna_report_tsgene_loss
    eval_cna_gain <- report_data$cna_report_oncogene_gain
    eval_cna_biomarker <- report_data$cna_report_biomarkers
    eval_msi_report <- report_data$msi_report
    eval_missing_msi_data <- report_data$missing_msi_data

    if(show_noncoding == FALSE){
      eval_tier5 <- FALSE
    }


    rmarkdown::render(system.file("templates","report.Rmd", package="pcgrr2"), output_file = paste0(sample_name,'.pcgr.html'), output_dir = project_directory, intermediates_dir = project_directory, params = list(logR_gain = logR_gain, logR_homdel = logR_homdel, tier1_report = report_data$tier1_report, tier2_report = report_data$tier2_report, tier3_report = report_data$tier3_report, tier4_report = report_data$tier4_report, tier5_report = report_data$tier5_report, msi_report = report_data$msi_report, missing_msi_data = report_data$missing_msi_data, cna_report_tsgene_loss = report_data$cna_report_tsgene_loss, cna_report_oncogene_gain = report_data$cna_report_oncogene_gain, cna_report_biomarkers = report_data$cna_report_biomarkers, cna_report_segments = report_data$cna_report_segments, signature_report = report_data$signature_report, missing_signature_data = report_data$missing_signature_data), quiet=T)
  }

}

#' Function that generates a data frame with basic biomarker annotations from tier1 variants
#'
#' @param tier1_variants df with tier 1 variants
#' @param sample_id Sample identifier
#'
#' @return tsv_variants data frame with all tier 1 biomarkers for tab-separated output
#'
generate_biomarker_tsv <- function(tier1_variants, sample_name = 'test'){

  bm_tags <- c('BM_CLINICAL_SIGNIFICANCE','BM_EVIDENCE_LEVEL','BM_EVIDENCE_TYPE','BM_EVIDENCE_DIRECTION','BM_CANCER_TYPE', 'BM_THERAPEUTIC_CONTEXT','BM_RATING','BM_CITATION')
  all_biomarker_tags <- c(c('GENOMIC_CHANGE','GENOME_VERSION','VCF_SAMPLE_ID','SYMBOL','CONSEQUENCE'),bm_tags)
  tier1_tsv <- tier1_variants
  tsv_biomarkers <- NULL
  if(nrow(tier1_tsv) > 0){
    tier1_tsv$VCF_SAMPLE_ID <- sample_name
    tier1_tsv <- tier1_tsv %>% dplyr::select(dplyr::one_of(all_biomarker_tags))
    tier1_tsv$TIER <- 'TIER 1'
    tier1_tsv$TIER_DESCRIPTION <- 'Clinical biomarker - Predictive/prognostic/diagnostic/predisposing'

    tmp2 <- as.data.frame(tier1_tsv %>% dplyr::rowwise() %>% dplyr::mutate(BM_CITATION2 = paste(unlist(stringr::str_replace_all(stringr::str_match_all(BM_CITATION,">.+<"),"^>|<$","")),collapse =";")))

    tier1_tsv$BM_CITATION <- tmp2$BM_CITATION2
    tsv_biomarkers <- tier1_tsv %>% dplyr::distinct()
  }
  return(tsv_biomarkers)
}

#' Function that generates dense and tiered annotated variant datasets
#'
#' @param tier1_variants df with tier 1 variants
#' @param tier2_variants df with tier 2 variants
#' @param tier3_variants df with tier 3 variants
#' @param tier4_variants df with tier 4 variants
#' @param tier5_variants df with tier 5 variants
#' @param pcgr_data List of data frames with PCGR data annotations
#' @param sample_name Sample identifier
#'
#' @return tsv_variants data frame with tier-annotated list of variants for tab-separated output
#'
generate_tier_tsv <- function(tier1_variants, tier2_variants, tier3_variants, tier4_variants, tier5_variants, pcgr_data, sample_name = 'test'){

  rlogging::message("Generating tiered set of result variants for output in tab-separated values (TSV) file")
  bm_tags <- c('BM_CLINICAL_SIGNIFICANCE','BM_EVIDENCE_LEVEL','BM_EVIDENCE_TYPE','BM_EVIDENCE_DIRECTION','BM_CANCER_TYPE','BM_THERAPEUTIC_CONTEXT','BM_RATING','BM_CITATION')
  tier1_tsv <- tier1_variants
  tsv_variants <- NULL
  if(nrow(tier1_tsv) > 0){
    tier1_tsv <- as.data.frame(tier1_tsv %>% dplyr::select(-dplyr::one_of(bm_tags)))
    tier1_tsv$TIER <- 'TIER 1'
    tier1_tsv$TIER_DESCRIPTION <- 'Clinical biomarker - Predictive/prognostic/diagnostic/predisposing'
    tier1_tsv$VCF_SAMPLE_ID <- sample_name
    #tier1_tsv <- unique(tier1_tsv)
    tier1_tsv_unique <- tier1_tsv %>% dplyr::distinct()
    tsv_variants <- rbind(tsv_variants, dplyr::select(tier1_tsv_unique, dplyr::one_of(pcgr_data$pcgr_tsv_tiered_columns)))
  }
  tier2_tsv <- tier2_variants
  if(nrow(tier2_tsv) > 0){
    tier2_tsv$TIER <- 'TIER 2'
    tier2_tsv$TIER_DESCRIPTION <- 'Other cancer mutation hotspot/predicted driver mutation/curated cancer-associated mutation'
    tier2_tsv$VCF_SAMPLE_ID <- sample_name
    tsv_variants <- rbind(tsv_variants, dplyr::select(tier2_tsv, dplyr::one_of(pcgr_data$pcgr_tsv_tiered_columns)))
  }
  tier3_tsv <- tier3_variants
  if(nrow(tier3_tsv) > 0){
    tier3_tsv$TIER <- 'TIER 3'
    tier3_tsv$TIER_DESCRIPTION <- 'Other cancer census gene/proto-oncogene/tumor suppressor mutation'
    tier3_tsv$VCF_SAMPLE_ID <- sample_name
    tsv_variants <- rbind(tsv_variants, dplyr::select(tier3_tsv, dplyr::one_of(pcgr_data$pcgr_tsv_tiered_columns)))
  }
  tier4_tsv <- tier4_variants
  if(nrow(tier4_tsv) > 0){
    tier4_tsv$TIER <- 'TIER 4'
    tier4_tsv$VCF_SAMPLE_ID <- sample_name
    tier4_tsv$TIER_DESCRIPTION <- 'Other coding mutation'
    tsv_variants <- rbind(tsv_variants, dplyr::select(tier4_tsv, dplyr::one_of(pcgr_data$pcgr_tsv_tiered_columns)))
  }
  tier5_tsv <- tier5_variants
  if(nrow(tier5_tsv) > 0){
    tier5_tsv$TIER <- 'TIER 5'
    tier5_tsv$VCF_SAMPLE_ID <- sample_name
    tier5_tsv$TIER_DESCRIPTION <- 'Non-coding mutation'
    tsv_variants <- rbind(tsv_variants, dplyr::select(tier5_tsv, dplyr::one_of(pcgr_data$pcgr_tsv_tiered_columns)))
  }
  tsv_variants$COSMIC <- unlist(lapply(stringr::str_match_all(tsv_variants$COSMIC,"COSM[0-9]{1,}"),paste,collapse=","))
  tsv_variants$DBSNP <- unlist(lapply(stringr::str_match_all(tsv_variants$DBSNP,">rs[0-9]{1,}<"),paste,collapse=","))
  tsv_variants$DBSNP <- stringr::str_replace_all(tsv_variants$DBSNP,">|<","")
  tsv_variants$GENE_NAME <- unlist(lapply(stringr::str_match_all(tsv_variants$GENE_NAME,">.+<"),paste,collapse=","))
  tsv_variants$GENE_NAME <- stringr::str_replace_all(tsv_variants$GENE_NAME,">|<","")
  tsv_variants$CLINVAR <- unlist(lapply(stringr::str_match_all(tsv_variants$CLINVAR,">.+<"),paste,collapse=","))
  tsv_variants$CLINVAR <- stringr::str_replace_all(tsv_variants$CLINVAR,">|<","")
  tsv_variants$PROTEIN_DOMAIN <- unlist(lapply(stringr::str_match_all(tsv_variants$PROTEIN_DOMAIN,">.+<"),paste,collapse=","))
  tsv_variants$PROTEIN_DOMAIN <- stringr::str_replace_all(tsv_variants$PROTEIN_DOMAIN,">|<","")

  tsv_variants <- tsv_variants %>% dplyr::distinct()

  return(tsv_variants)
}


#' Function that generates report data for tiered precision oncology report
#'
#' @param sample_calls data frame with list of variant calls
#' @param pcgr_data List of data frames with PCGR annotations
#' @param sample_id sample identifier
#' @param minimum_n_signature_analysis minimum number of mutations for signature analysis
#' @param signatures_limit limit the number of possible mutational signatures
#'
#' @return report_data data frame with all report elements
#'
generate_report_data <- function(sample_calls, pcgr_data, sample_name = NULL,  minimum_n_signature_analysis = 50, signature_normalization_method = 'default', signatures_limit = 6, predict_MSI = FALSE, identify_msigs = FALSE){

  rlogging::message("Generating data for tiered cancer genome report")
  tier1_report <- TRUE
  tier2_report <- TRUE
  tier3_report <- TRUE
  tier4_report <- TRUE
  tier5_report <- TRUE
  msi_report <- FALSE
  missing_msi_data <- FALSE

  clinical_evidence_items_prognostic <- data.frame()
  clinical_evidence_items_diagnostic <- data.frame()
  clinical_evidence_items_predisposing <- data.frame()
  clinical_evidence_items_predictive <- data.frame()

  variants_tier1_display <- data.frame()
  variants_tier2_display <- data.frame()
  variants_tier2_hotspots <- data.frame()
  variants_tier2_curated_mutations <- data.frame()
  variants_tier2_predicted_drivers <- data.frame()
  variants_tier3_display <- data.frame()
  variants_tier4_display <- data.frame()
  variants_tier5_display <- data.frame()
  tsv_variants <- data.frame()
  tsv_biomarkers <- data.frame()
  maf_df <- data.frame()
  msi_prediction_data <- list()

  signature_report <- FALSE
  signature_data <- NULL
  missing_signature_data <- FALSE
  signature_call_set <- data.frame()
  min_variants_for_signature <- minimum_n_signature_analysis

  biomarker_descriptions <- data.frame()
  sample_calls_SNVs <- data.frame()
  sample_calls_INDELs <- data.frame()
  sample_calls_coding <- data.frame()
  sample_calls_noncoding <- data.frame()

  if(nrow(sample_calls) > 0){

    if(predict_MSI == TRUE){
      if(nrow(sample_calls) > 30){
        msi_prediction_data <- pcgrr2::predict_msi_status(sample_calls, simpleRepeats_gr = pcgr_data$simpleRepeats_gr, windowMasker_gr = pcgr_data$windowMasker_gr, msi_prediction_model = pcgr_data$msi_prediction_model, indelFracPlot_template = pcgr_data$indelFracPlot_msi_report, sample_name = sample_name)
        msi_report <- TRUE
      }
      else{
        missing_msi_data <- TRUE
        msi_report <- FALSE
      }
    }

    sample_calls_coding <- sample_calls %>% dplyr::filter(stringr::str_detect(CONSEQUENCE,"stop_gained|stop_lost|start_lost|frameshift_variant|missense_variant|splice_donor|splice_acceptor|inframe_deletion|inframe_insertion"))
    rlogging::message(paste0("Number of coding variants: ",nrow(sample_calls_coding)))
    sample_calls_noncoding <- sample_calls %>% dplyr::filter(!stringr::str_detect(CONSEQUENCE,"stop_gained|stop_lost|start_lost|frameshift_variant|missense_variant|splice_donor|splice_acceptor|inframe_deletion|inframe_insertion"))
    rlogging::message(paste0("Number of noncoding variants: ",nrow(sample_calls_noncoding)))
    sample_calls_SNVs <- sample_calls %>% dplyr::filter(VARIANT_CLASS == 'SNV')
    sample_calls_INDELs <- sample_calls %>% dplyr::filter(VARIANT_CLASS != 'SNV')

    if(identify_msigs == TRUE){
      if(any(grepl(paste0("VARIANT_CLASS$"),names(sample_calls)))){
        if(nrow(sample_calls[sample_calls$VARIANT_CLASS == 'SNV',]) >= minimum_n_signature_analysis){
          signature_call_set <- sample_calls[sample_calls$VARIANT_CLASS == 'SNV',]
          signature_call_set <- dplyr::filter(signature_call_set, CHROM != 'MT')
          signature_call_set$VCF_SAMPLE_ID <- sample_name
          signature_report <- TRUE

          mut_signature_contributions <- pcgrr2::signature_contributions_single_sample(signature_call_set, sample_name = sample_name, normalization_method = signature_normalization_method, signatures_limit = signatures_limit)

          signature_columns <- as.numeric(stringr::str_replace(as.character(mut_signature_contributions$which_signatures_df[mut_signature_contributions$which_signatures_df$signature_id != 'unknown',]$signature_id),"S",""))

          weight_df <- data.frame('Signature_ID' = as.character(mut_signature_contributions$which_signatures_df$signature_id), 'Weight' = round(as.numeric(mut_signature_contributions$which_signatures_df$weight),digits=3), stringsAsFactors = F)

          cancertypes_aetiologies <- pcgr_data$signatures_aetiologies[signature_columns,]
          signatures_cancertypes_aetiologies <- dplyr::left_join(cancertypes_aetiologies,weight_df,by=c("Signature_ID")) %>% dplyr::arrange(desc(Weight))
          signatures_cancertypes_aetiologies <- signatures_cancertypes_aetiologies[,c("Signature_ID","Weight","Cancer_types","Proposed_aetiology","Comments")]
          signatures_cancertypes_aetiologies$Trimer_normalization_method <- signature_normalization_method

          unknown_df <- dplyr::filter(mut_signature_contributions$which_signatures_df, signature_id == 'unknown')
          unknown_df <- data.frame(Signature_ID = unknown_df$signature_id, Weight = unknown_df$weight, stringsAsFactors = F)
          unknown_df$Weight <- round(as.numeric(unknown_df$Weight),digits=3)
          unknown_df$Cancer_types <- NA
          unknown_df$Proposed_aetiology <- NA
          unknown_df$Comments <- NA
          unknown_df$Trimer_normalization_method <- signature_normalization_method

          signatures_cancertypes_aetiologies <- rbind(signatures_cancertypes_aetiologies, unknown_df) %>% dplyr::arrange(desc(Weight))
          signature_data <- list('signature_call_set' = signature_call_set, 'mut_signature_contributions' = mut_signature_contributions, 'signatures_cancertypes_aetiologies' = signatures_cancertypes_aetiologies)

        }
        else{
          if(nrow(sample_calls[sample_calls$VARIANT_CLASS == 'SNV',]) > 0){
            signature_call_set <- sample_calls[sample_calls$VARIANT_CLASS == 'SNV',]
            signature_call_set <- dplyr::filter(signature_call_set, CHROM != 'MT')
          }
          rlogging::message(paste0("Too few variants (n = ",nrow(signature_call_set),") for reconstruction of mutational signatures by deconstructSigs"))
          missing_signature_data <- TRUE
          signature_data <- list('signature_call_set' = signature_call_set)

        }
      }
    }

    ## Analyze Tier1: actionable mutations and variants of clinical significance (diagnosis/prognosis etc)
    variants_tier1 <- pcgrr2::get_clinical_associations_civic_cbmdb(sample_calls_coding, pcgr_data = pcgr_data)
    if(nrow(variants_tier1) > 0){
      variants_tier1 <- dplyr::arrange(variants_tier1, EVIDENCE_LEVEL)
      clinical_evidence_items_prognostic <- dplyr::select(variants_tier1, dplyr::one_of(pcgr_data$tier1_tags_display)) %>% dplyr::filter(EVIDENCE_TYPE == 'Prognostic')
      clinical_evidence_items_diagnostic <- dplyr::select(variants_tier1, dplyr::one_of(pcgr_data$tier1_tags_display)) %>% dplyr::filter(EVIDENCE_TYPE == 'Diagnostic')
      clinical_evidence_items_predictive <- dplyr::select(variants_tier1, dplyr::one_of(pcgr_data$tier1_tags_display)) %>% dplyr::filter(EVIDENCE_TYPE == 'Predictive')
      clinical_evidence_items_predisposing <- dplyr::select(variants_tier1, dplyr::one_of(pcgr_data$tier1_tags_display)) %>% dplyr::filter(EVIDENCE_TYPE == 'Predisposing')
      variants_tier1_display <- variants_tier1 %>% dplyr::select(GENOMIC_CHANGE) %>% dplyr::distinct()
      tier1_report <- TRUE
      variants_tier1 <- dplyr::rename(variants_tier1, BM_CLINICAL_SIGNIFICANCE = CLINICAL_SIGNIFICANCE, BM_EVIDENCE_LEVEL = EVIDENCE_LEVEL, BM_EVIDENCE_TYPE = EVIDENCE_TYPE, BM_EVIDENCE_DIRECTION = EVIDENCE_DIRECTION, BM_CANCER_TYPE = CANCER_TYPE, BM_THERAPEUTIC_CONTEXT = THERAPEUTIC_CONTEXT, BM_CITATION = CITATION, BM_RATING = RATING)
    }

    ## Analyze Tier 2: curated mutations, cancer mutation hotspots and predicted driver mutations
    variants_tier2 <- dplyr::select(sample_calls_coding, dplyr::one_of(pcgr_data$pcgr_all_annotation_columns))
    variants_tier2 <- variants_tier2 %>% dplyr::filter(!is.na(INTOGEN_DRIVER_MUT) | !is.na(CANCER_MUTATION_HOTSPOT) | !is.na(OTHER_DISEASE_DOCM))
    if(nrow(variants_tier1) > 0){
      variants_tier2 <- dplyr::anti_join(variants_tier2, variants_tier1_display, by=c("GENOMIC_CHANGE"))
    }
    tier12 <- variants_tier1_display
    if(nrow(variants_tier2) > 0){
      variants_tier2 <- variants_tier2 %>% dplyr::arrange(desc(ONCOSCORE))
      tier12 <- rbind(variants_tier1_display,dplyr::select(variants_tier2,GENOMIC_CHANGE)) %>% dplyr::distinct()
      variants_tier2_display <- dplyr::select(variants_tier2, dplyr::one_of(pcgr_data$tier2_tags_display))
      variants_tier2_hotspots <- variants_tier2_display %>% dplyr::filter(!is.na(CANCER_MUTATION_HOTSPOT))
      variants_tier2_curated_mutations <- variants_tier2_display %>% dplyr::filter(is.na(CANCER_MUTATION_HOTSPOT) & !is.na(OTHER_DISEASE_DOCM))
      variants_tier2_predicted_drivers <- variants_tier2_display %>% dplyr::filter(is.na(CANCER_MUTATION_HOTSPOT) & is.na(OTHER_DISEASE_DOCM) & !is.na(INTOGEN_DRIVER_MUT))
    }

    ## Analyze Tier 3: coding mutations in oncogenes/tumor suppressors/cancer census genes
    variants_tier3 <- dplyr::select(sample_calls_coding, dplyr::one_of(pcgr_data$pcgr_all_annotation_columns))
    variants_tier3 <- variants_tier3 %>% dplyr::filter(!is.na(CANCER_CENSUS_SOMATIC) | ONCOGENE == TRUE | TUMOR_SUPPRESSOR == TRUE)
    if(nrow(tier12) > 0){
      variants_tier3 <- dplyr::anti_join(variants_tier3,tier12, by=c("GENOMIC_CHANGE"))
    }
    tier123 <- tier12
    if(nrow(variants_tier3) > 0){
      variants_tier3 <- variants_tier3 %>% dplyr::arrange(desc(ONCOSCORE))
      tier123 <- rbind(tier12,dplyr::select(variants_tier3,GENOMIC_CHANGE)) %>% dplyr::distinct()
      variants_tier3_display <- dplyr::select(variants_tier3, dplyr::one_of(pcgr_data$tier3_tags_display))
    }

    ## Analyze Tier 4: Other coding mutations
    variants_tier4 <- dplyr::select(sample_calls_coding, dplyr::one_of(pcgr_data$pcgr_all_annotation_columns))
    if(nrow(tier123) > 0){
      variants_tier4 <- dplyr::anti_join(variants_tier4,tier123, by=c("GENOMIC_CHANGE"))
    }
    if(nrow(variants_tier4) > 0){
      variants_tier4 <- variants_tier4 %>% dplyr::arrange(desc(ONCOSCORE))
      variants_tier4_display <- dplyr::select(variants_tier4, dplyr::one_of(pcgr_data$tier4_tags_display))
    }

    ## Analyze Tier 5: Non-coding mutations
    variants_tier5 <- dplyr::select(sample_calls_noncoding, dplyr::one_of(pcgr_data$pcgr_all_annotation_columns))
    if(nrow(variants_tier5) > 0){
      variants_tier5 <- variants_tier5 %>% dplyr::arrange(desc(ONCOSCORE))
      variants_tier5_display <- dplyr::select(variants_tier5, dplyr::one_of(pcgr_data$tier5_tags_display))
    }

    tsv_variants <- NULL
    tsv_biomarkers <- NULL
    tsv_variants <- pcgrr2::generate_tier_tsv(variants_tier1,variants_tier2,variants_tier3, variants_tier4, variants_tier5,pcgr_data = pcgr_data, sample_name = sample_name)

    maf_df <- pcgrr2::tier_to_maf(tsv_variants)
    tsv_biomarkers <- pcgrr2::generate_biomarker_tsv(variants_tier1, sample_name = sample_name)
  }
  else{
    missing_signature_data <- TRUE
    if(predict_MSI == T){
      msi_report <- FALSE
      missing_msi_data <- TRUE
    }

  }
  report_data <- list('tier1_report' = tier1_report, 'tier2_report' = tier2_report, 'tier3_report' = tier3_report, 'tier4_report' = tier4_report, 'tier5_report' = tier5_report, 'clinical_evidence_items_prognostic' = clinical_evidence_items_prognostic, 'clinical_evidence_items_diagnostic' = clinical_evidence_items_diagnostic, 'clinical_evidence_items_predisposing' = clinical_evidence_items_predisposing, 'clinical_evidence_items_predictive' = clinical_evidence_items_predictive, 'tsv_variants' = tsv_variants, 'tsv_biomarkers' = tsv_biomarkers, 'variants_tier1_display' = variants_tier1_display, 'variants_tier2_display' = variants_tier2_display, 'variants_tier2_hotspots' = variants_tier2_hotspots, 'variants_tier2_curated_mutations' = variants_tier2_curated_mutations, 'variants_tier2_predicted_drivers' = variants_tier2_predicted_drivers, 'variants_tier3_display' = variants_tier3_display, 'variants_tier4_display' = variants_tier4_display,'variants_tier5_display' = variants_tier5_display, 'msi_prediction_data' = msi_prediction_data, 'msi_report' = msi_report, 'missing_msi_data' = missing_msi_data, 'sample_calls'= sample_calls, 'sample_calls_coding' = sample_calls_coding, 'sample_calls_noncoding' = sample_calls_noncoding, 'sample_calls_SNVs' = sample_calls_SNVs, 'sample_calls_INDELs' = sample_calls_INDELs, 'signature_report' = signature_report, 'missing_signature_data' = missing_signature_data, 'signature_data' = signature_data, 'signatures_limit' = signatures_limit, 'maf_df' = maf_df, 'sample_name' = sample_name, 'cna_report_oncogene_gain' = FALSE, 'cna_report_tsgene_loss' = FALSE, 'cna_report_biomarkers' = FALSE, 'cna_report_segments' = FALSE)


  return(report_data)

}

#' Function that adds HTML links to different genetic variant identifiers
#'
#' @param var_df data frame with variants
#' @param vardb type of variant database
#' @param linktype type of link
#'
#' @return var_df
#'
annotate_variant_link <- function(var_df, vardb = "DBSNP", linktype = "dbsource"){
  if(vardb == 'DBSNP'){

    if(any(grepl(paste0("^DBSNPRSID$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df)))){
      var_df_unique_slim <- dplyr::select(var_df, DBSNPRSID, VAR_ID) %>% dplyr::filter(!is.na(DBSNPRSID)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        var_df_unique_slim_melted <- var_df_unique_slim %>% tidyr::separate_rows(DBSNPRSID,sep=",")
        if(linktype == "dbsource"){
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>% dplyr::mutate(tmp_dbsnp = paste0("<a href='http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",DBSNPRSID,"' target=\"_blank\">rs",DBSNPRSID,"</a>"))
        }

        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>% dplyr::summarise(DBSNPLINK = unlist(paste(tmp_dbsnp, collapse = ", ")))
        var_df_links <- dplyr::select(var_df_links, VAR_ID, DBSNPLINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$DBSNPLINK <- NA
      }
    }
    else{
      cat("WARNING: Could not generate DBSNP links - no DBSNPRSID provided in annotated VCF",sep="\n")
      var_df$DBSNPLINK <- NA
    }
  }

  if(vardb == 'DBNSFP'){
    if(any(grepl(paste0("EFFECT_PREDICTIONS"),names(var_df)))){
      var_df$PREDICTED_EFFECT <- var_df$EFFECT_PREDICTIONS
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"metalr:","<a href='https://www.ncbi.nlm.nih.gov/pubmed/25552646' target=\"_blank\">Ensembl-LogisticRegression</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"metasvm:","<a href='https://www.ncbi.nlm.nih.gov/pubmed/25552646' target=\"_blank\">Ensembl-SVM</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"mutationassessor:","<a href='http://mutationassessor.org' target=\"_blank\">MutationAssessor</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"mutationtaster:","<a href='http://www.mutationtaster.org' target=\"_blank\">MutationTaster</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"fathmm:","<a href='http://fathmm.biocompute.org.uk' target=\"_blank\">FATHMM</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"fathmm_mkl_coding:","<a href='http://fathmm.biocompute.org.uk/fathmmMKL.htm' target=\"_blank\">FATHMM-mkl</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"sift:","<a href='http://provean.jcvi.org/index.php' target=\"_blank\">SIFT</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"polyphen2_hdiv:","<a href='http://genetics.bwh.harvard.edu/pph2/' target=\"_blank\">PolyPhen2-HDIV</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"polyphen2_hvar:","<a href='http://genetics.bwh.harvard.edu/pph2/' target=\"_blank\">PolyPhen2-HVAR</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"sift:","<a href='http://genetics.bwh.harvard.edu/pph2/' target=\"_blank\">SIFT</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"lrt:","<a href='http://www.genetics.wustl.edu/jflab/lrt_query.html' target=\"_blank\">LRT</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"provean:","<a href='http://provean.jcvi.org/index.php' target=\"_blank\">PROVEAN</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"mutpred:","<a href='http://mutpred.mutdb.org' target=\"_blank\">MutPred</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"m-cap:","<a href='http://bejerano.stanford.edu/MCAP/' target=\"_blank\">M-CAP</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"splice_site_rf:","<a href='http://nar.oxfordjournals.org/content/42/22/13534' target=\"_blank\">Splice site effect (Random forest)</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"splice_site_ada:","<a href='http://nar.oxfordjournals.org/content/42/22/13534' target=\"_blank\">Splice site effect (Adaptive boosting)</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"cadd_phred:","<a href='http://cadd.gs.washington.edu' target=\"_blank\">CADD Phred</a>:")
      var_df$PREDICTED_EFFECT <- stringr::str_replace(var_df$PREDICTED_EFFECT,"gerp_rs:","<a href='http://mendel.stanford.edu/SidowLab/downloads/gerp/' target=\"_blank\">GERP++ RS score</a>:")
    }
    else{
      var_df$PREDICTED_EFFECT <- NA
    }
  }

  if(vardb == 'CLINVAR'){

    if(any(grepl(paste0("^CLINVAR_MSID$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df))) & any(grepl(paste0("^CLINVAR_TRAITS_ALL$"),names(var_df)))){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, CLINVAR_MSID, CLINVAR_TRAITS_ALL) %>% dplyr::filter(!is.na(CLINVAR_MSID)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        if(linktype == "dbsource"){
          var_df_unique_slim <- var_df_unique_slim %>% dplyr::mutate(CLINVARLINK = paste0("<a href='http://www.ncbi.nlm.nih.gov/clinvar/variation/",CLINVAR_MSID,"' target=\"_blank\">",CLINVAR_TRAITS_ALL,"</a>"))
        }
        var_df_links <- var_df_unique_slim %>% dplyr::select(VAR_ID, CLINVARLINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$CLINVARLINK <- NA
      }
    }
    else{
      cat("WARNING: Could not generate CLINVAR links - no CLINVAR_MSID provided in annotated VCF",sep="\n")
      var_df$CLINVARLINK <- NA
    }
  }

  if(vardb == 'NCBI_GENE'){

    if(any(grepl(paste0("^GENENAME$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df))) & any(grepl(paste0("^ENTREZ_ID$"),names(var_df)))){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, GENENAME, ENTREZ_ID) %>% dplyr::filter(!is.na(ENTREZ_ID)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        if(linktype == "dbsource"){
          var_df_unique_slim <- var_df_unique_slim %>% dplyr::mutate(NCBI_GENE_LINK = paste0("<a href='http://www.ncbi.nlm.nih.gov/gene/",ENTREZ_ID,"' target=\"_blank\">",GENENAME,"</a>"))
        }
        var_df_links <- var_df_unique_slim %>% dplyr::select(VAR_ID, NCBI_GENE_LINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$NCBI_GENE_LINK <- NA
      }
    }
    else{
      cat("WARNING: Could not generate NCBI_GENE links - no ENTREZ_ID provided in annotated VCF",sep="\n")
      var_df$NCBI_GENE_LINK <- NA
    }
  }

  if(vardb == 'COSMIC'){

    if(any(grepl(paste0("^COSMIC_MUTATION_ID$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df)))){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, COSMIC_MUTATION_ID) %>% dplyr::filter(!is.na(COSMIC_MUTATION_ID)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        var_df_unique_slim_melted <- var_df_unique_slim %>% tidyr::separate_rows(COSMIC_MUTATION_ID,sep="&")
        if(linktype == "dbsource"){
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>% dplyr::mutate(tmp_cosmic = paste0("<a href='http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=",stringr::str_replace(COSMIC_MUTATION_ID,"COSM",""),"' target=\"_blank\">",COSMIC_MUTATION_ID,"</a>"))
        }
        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>% dplyr::summarise(COSMICLINK = unlist(paste(tmp_cosmic, collapse = ", ")))
        var_df_links <- dplyr::select(var_df_links, VAR_ID, COSMICLINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$COSMICLINK <- NA
      }
    }
    else{
      cat("WARNING: Could not add COSMIC links - no COSMIC_MUTATION_ID provided in annotated VCF",sep="\n")
      var_df$COSMICLINK <- NA
    }
  }

  if(vardb == 'DGIDB'){
    if(any(grepl(paste0("^ANTINEOPLASTIC_DRUG_INTERACTION$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df)))){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, ANTINEOPLASTIC_DRUG_INTERACTION) %>% dplyr::filter(!is.na(ANTINEOPLASTIC_DRUG_INTERACTION)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        var_df_unique_slim_melted <- var_df_unique_slim %>% tidyr::separate_rows(ANTINEOPLASTIC_DRUG_INTERACTION,sep="&")
        if(linktype == "dbsource"){
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>% dplyr::mutate(tmp_dgidb = paste0("<a href='http://dgidb.genome.wustl.edu/drugs/",stringr::str_replace(ANTINEOPLASTIC_DRUG_INTERACTION,"\\S{1,}:",""),"' target=\"_blank\">",ANTINEOPLASTIC_DRUG_INTERACTION,"</a>"))
        }

        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>% dplyr::summarise(DGIDBLINK = unlist(paste(tmp_dgidb, collapse = ", ")))
        var_df_links <- dplyr::select(var_df_links, VAR_ID, DGIDBLINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$DGIDBLINK <- NA
      }
    }
    else{
      cat("WARNING: Could not generate DGIdb links - no DGIDB info provided in annotated VCF",sep="\n")
      var_df$DGIDBLINK <- NA
    }
  }

  return(var_df)

}

#' A function that converts the INFO tags in a VRanges object into a basic types
#'
#' @param vr VRanges object
#' @return vr Vranges object with INFO tags more simply formatted
#'

postprocess_vranges_info <- function(vr){
  ## Convert IntegerLists and CharacterLists to basic character lists
  vcf_annotations_df <- NULL
  for(tag in colnames(GenomicRanges::mcols(vr))){
    mcol_class <- class(GenomicRanges::mcols(vr)[,c(tag)])[1]

    if(mcol_class != "character" & mcol_class != "integer" & mcol_class != "logical" & mcol_class != "numeric"){
      annotation_track <- NULL
      #cat("TAG: ",tag, ', type:',mcol_class,'\n')
      if(mcol_class == "CompressedCharacterList"){
        annotation_track <- data.frame(val = as.character(Biostrings::unstrsplit(GenomicRanges::mcols(vr)[[tag]], sep=',')))
      }else{
        annotation_track <- data.frame(val = as.character(sapply(GenomicRanges::mcols(vr)[[tag]], paste, collapse=",")))
      }
      if(is.null(vcf_annotations_df)){
        vcf_annotations_df <- data.frame(annotation_track$val)
        names(vcf_annotations_df) <- c(tag)
        vcf_annotations_df[,c(tag)] <- as.character(vcf_annotations_df[,c(tag)])
      }
      else{
        vcf_annotations_df[,c(tag)] <- as.character(annotation_track$val)
      }
      ## add NA to empty values
      if(nrow(as.data.frame(vcf_annotations_df[nchar(vcf_annotations_df[,c(tag)]) == 0,])) != 0){
        if(dim(vcf_annotations_df)[2] == 1){
          vcf_annotations_df[nchar(vcf_annotations_df[,c(tag)]) == 0,] <- NA
        }
        else{
          vcf_annotations_df[nchar(vcf_annotations_df[,c(tag)]) == 0,][,c(tag)] <- NA
        }
      }
    }
    else{
      #cat("TAG: ",tag, ', type:',mcol_class,'\n')
      if(is.null(vcf_annotations_df)){
        vcf_annotations_df <- data.frame(GenomicRanges::mcols(vr)[,c(tag)])
        names(vcf_annotations_df) <- c(tag)
      }
      else{
        vcf_annotations_df[,c(tag)] <- GenomicRanges::mcols(vr)[,c(tag)]
      }
    }
  }

  ## add variant_id and sample_id
  position <- GenomicRanges::start(GenomicRanges::ranges(vr))
  vcf_annotations_df['VAR_ID'] <- paste(as.character(GenomeInfoDb::seqnames(vr)),position,VariantAnnotation::ref(vr),VariantAnnotation::alt(vr),sep="_")
  #vcf_annotations_df['VCF_SAMPLE_ID'] <- as.character(VariantAnnotation::sampleNames(vr))
  GenomicRanges::mcols(vr) <- S4Vectors::DataFrame(vcf_annotations_df)

  return(vr)

}

#' Function that adds PFAM name descriptions to PFAM identifiers
#'
#' @param vcf_data_df
#'
#' @return vcf_data_df_pfam
#'
add_pfam_domain_links <- function(vcf_data_df, pcgr_data){

  rlogging::message("Extending annotation descriptions related to PFAM protein domains")
  if("DOMAINS" %in% colnames(vcf_data_df) & "VAR_ID" %in% colnames(vcf_data_df)){
    pfam_df <- dplyr::select(vcf_data_df,DOMAINS,VAR_ID) %>% dplyr::filter(!is.na(DOMAINS))
    if(nrow(pfam_df) == 0){
      vcf_data_df$PROTEIN_DOMAIN <- NA
      return(vcf_data_df)
    }
    pfam_df <- pfam_df %>% dplyr::distinct() %>% tidyr::separate_rows(DOMAINS,sep="&") %>% dplyr::filter(stringr::str_detect(DOMAINS,"Pfam_domain"))
    pfam_df$DOMAINS <- stringr::str_replace(pfam_df$DOMAINS,"Pfam_domain:","")
    pfam_df <- dplyr::left_join(pfam_df,pcgr_data$pfam_domains,by=c("DOMAINS" = "pfam_id")) %>% dplyr::select(VAR_ID,url)
    pfam_df <- dplyr::rename(pfam_df, PD = url)
    pfam_ret <- as.data.frame(dplyr::group_by(pfam_df, VAR_ID) %>% dplyr::summarise(PROTEIN_DOMAIN = paste(PD, collapse=", ")))

    if(nrow(pfam_ret) > 0){
      vcf_data_df <- dplyr::left_join(vcf_data_df,pfam_ret,by=c("VAR_ID" = "VAR_ID"))
    }
    else{
      vcf_data_df$PROTEIN_DOMAIN <- NA
    }
  }
  else{
    vcf_data_df$PROTEIN_DOMAIN <- NA
  }

  return(vcf_data_df)
}


#' Function that appends clinical annotations for somatic cancer variants
#'
#' @param vcf_data_df data frame with variants
#' @param mapping - one of 'exact' (allele-specific) or 'approximate' (codon or exon-level biomarkers)
#' @param ncgc - logical indicating whether NCGC-specific tags are to be appended
#'
#' @return vcf_data_df
#'
get_clinical_associations_civic_cbmdb <- function(sample_calls, pcgr_data){

  rlogging::message("Looking up biomarkers for precision oncology")
  civic_biomarkers <- pcgr_data$civic_biomarkers
  cbmdb_biomarkers <- pcgr_data$cbmdb_biomarkers
  if("pubmed_html_link" %in% colnames(civic_biomarkers)){
    civic_biomarkers <- dplyr::rename(civic_biomarkers, citation = pubmed_html_link)
  }
  if("evidence_description" %in% colnames(civic_biomarkers)){
    civic_biomarkers <- dplyr::rename(civic_biomarkers, description = evidence_description)
  }
  if("pubmed_html_link" %in% colnames(cbmdb_biomarkers)){
    cbmdb_biomarkers <- dplyr::rename(cbmdb_biomarkers, citation = pubmed_html_link)
  }
  if("evidence_description" %in% colnames(cbmdb_biomarkers)){
    cbmdb_biomarkers <- dplyr::rename(cbmdb_biomarkers, description = evidence_description)
  }
  biomarker_descriptions <- data.frame()

  #cbmdb_biomarkers <- dplyr::filter(cbmdb_biomarkers)
  #civic_biomarkers <- dplyr::filter(civic_biomarkers)
  civic_biomarkers <- dplyr::select(civic_biomarkers, -c(status,variant_id))
  all_eitems <- data.frame()

  for(mapping in c('exact','codon','exon')){
    clinical_evidence_items <- data.frame()
    if(mapping == 'exact'){
      sample_calls_civic <- sample_calls %>% dplyr::filter(!is.na(CIVIC_ID))
      if(nrow(sample_calls_civic) > 0){
        tmp <- dplyr::select(sample_calls_civic,CIVIC_ID,VAR_ID) %>% tidyr::separate_rows(CIVIC_ID,sep=",")
        sample_calls_civic <- merge(tmp,dplyr::select(sample_calls_civic,-c(CIVIC_ID)),by.x = "VAR_ID",by.y = "VAR_ID")
        civic_calls <- dplyr::select(sample_calls_civic,dplyr::one_of(pcgr_data$pcgr_all_annotation_columns))
        eitems <- dplyr::left_join(civic_calls,dplyr::filter(dplyr::select(civic_biomarkers,-c(civic_exon,civic_consequence,civic_codon,transvar_id,civic_id)),alteration_type == 'MUT'),by=c("CIVIC_ID" = "evidence_id"))
        names(eitems) <- toupper(names(eitems))
        clinical_evidence_items <- rbind(clinical_evidence_items, eitems)
      }
      sample_calls_cbmdb <- sample_calls %>% dplyr::filter(is.na(CIVIC_ID) & !is.na(CBMDB_ID))
      if(nrow(sample_calls_cbmdb) > 0){
        tmp <- dplyr::select(sample_calls_cbmdb,CBMDB_ID,VAR_ID) %>% tidyr::separate_rows(CBMDB_ID,sep=",")
        tmp$CBMDB_ID <- as.integer(tmp$CBMDB_ID)
        sample_calls_cbmdb <- merge(tmp,dplyr::select(sample_calls_cbmdb,-c(CBMDB_ID)),by.x = "VAR_ID",by.y = "VAR_ID")
        cbmdb_calls <- dplyr::select(sample_calls_cbmdb,dplyr::one_of(pcgr_data$pcgr_all_annotation_columns))
        eitems <- dplyr::left_join(cbmdb_calls,dplyr::filter(dplyr::select(cbmdb_biomarkers,-c(drug_family,transvar_id)),alteration_type == 'MUT'),by=c("CBMDB_ID" = "CBMDB_ID"))
        names(eitems) <- toupper(names(eitems))
        eitems <- eitems %>% dplyr::distinct()
        clinical_evidence_items <- rbind(clinical_evidence_items, eitems)
      }
    }
    else{
      sample_calls_civic <- sample_calls %>% dplyr::filter(!is.na(CIVIC_ID_2))
      if(nrow(sample_calls_civic) > 0){
        tmp <- dplyr::select(sample_calls_civic,CIVIC_ID_2,VAR_ID) %>% tidyr::separate_rows(CIVIC_ID_2,sep=",")
        sample_calls_civic <- merge(tmp,dplyr::select(sample_calls_civic,-c(CIVIC_ID_2)),by.x = "VAR_ID",by.y = "VAR_ID")
        civic_calls <- dplyr::select(sample_calls_civic,dplyr::one_of(pcgr_data$pcgr_all_annotation_columns))
        if(any(stringr::str_detect(civic_calls$CIVIC_ID_2,"EID"))){
          clinical_evidence_items <- dplyr::left_join(civic_calls,dplyr::filter(dplyr::select(civic_biomarkers,-civic_id),alteration_type == 'MUT'),by=c("CIVIC_ID_2" = "evidence_id"))
        }
        else{
          clinical_evidence_items <- dplyr::left_join(civic_calls,dplyr::filter(civic_biomarkers,alteration_type == 'MUT'),by=c("CIVIC_ID_2" = "civic_id"))
        }
        clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(mapping_category != 'gene')
        names(clinical_evidence_items) <- toupper(names(clinical_evidence_items))
        if(nrow(clinical_evidence_items) > 0){
          if(mapping == 'codon' & 'CIVIC_CODON' %in% colnames(clinical_evidence_items)){
            if(nrow(clinical_evidence_items[!is.na(clinical_evidence_items$CIVIC_CODON),]) > 0){
              clinical_evidence_items$CIVIC_CODON <- as.character(clinical_evidence_items$CIVIC_CODON)
              clinical_evidence_items <- dplyr::filter(clinical_evidence_items, !is.na(CIVIC_CODON))
              clinical_evidence_items$CODON <- as.character(stringr::str_replace_all(clinical_evidence_items$PROTEIN_CHANGE,"p\\.[A-Z]{1,}|[A-Z]|[A-Z]{1,}$|fs$|del$|dup$|delins[A-Z]{1,}$",""))
              clinical_evidence_items$CIVIC_CODON <- as.character(clinical_evidence_items$CIVIC_CODON)
              clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(startsWith(CODON,CIVIC_CODON))
              if(nrow(clinical_evidence_items) > 0){
                clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(!is.na(CIVIC_CONSEQUENCE) & startsWith(CONSEQUENCE,CIVIC_CONSEQUENCE) | is.na(CIVIC_CONSEQUENCE))
                if(nrow(clinical_evidence_items) > 0){
                  clinical_evidence_items <- clinical_evidence_items %>% dplyr::select(-c(EXON,CODON,CIVIC_CONSEQUENCE,MAPPING_CATEGORY,CIVIC_CODON,CIVIC_EXON))
                  names(clinical_evidence_items) <- toupper(names(clinical_evidence_items))
                }
                else{
                  clinical_evidence_items <- data.frame()
                }
              }
              else{
                clinical_evidence_items <- data.frame()
              }
            }
            else{
              clinical_evidence_items <- data.frame()
            }
          }
          if(mapping == 'exon' & 'CIVIC_EXON' %in% colnames(clinical_evidence_items)){
            if(nrow(clinical_evidence_items[!is.na(clinical_evidence_items$CIVIC_EXON),]) > 0){
              clinical_evidence_items <- dplyr::filter(clinical_evidence_items, !is.na(CIVIC_EXON))
              clinical_evidence_items$EXON <- as.integer(stringr::str_split_fixed(clinical_evidence_items$EXON,"/",2)[,1])
              clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(EXON == CIVIC_EXON)
              if(nrow(clinical_evidence_items) > 0){
                clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(!is.na(CIVIC_CONSEQUENCE) & startsWith(CONSEQUENCE,CIVIC_CONSEQUENCE) | is.na(CIVIC_CONSEQUENCE))
                if(nrow(clinical_evidence_items) > 0){
                  clinical_evidence_items <- clinical_evidence_items %>% dplyr::select(-c(EXON,CIVIC_CONSEQUENCE,MAPPING_CATEGORY,CIVIC_CODON,CIVIC_EXON))
                  names(clinical_evidence_items) <- toupper(names(clinical_evidence_items))
                }
                else{
                  clinical_evidence_items <- data.frame()
                }
              }
              else{
                clinical_evidence_items <- data.frame()
              }
            }
            else{
              clinical_evidence_items <- data.frame()
            }
          }
        }
      }

    }

    if(nrow(clinical_evidence_items) > 0){
      pcgr_all_annotation_columns_reduced <- pcgr_data$pcgr_all_annotation_columns[-which(pcgr_data$pcgr_all_annotation_columns == 'EXON' | pcgr_data$pcgr_all_annotation_columns == 'CIVIC_ID' | pcgr_data$pcgr_all_annotation_columns == 'CIVIC_ID_2' | pcgr_data$pcgr_all_annotation_columns == 'CBMDB_ID')]

      clinical_evidence_items$BIOMARKER_MAPPING <- mapping
      all_tier1_tags <- c(pcgr_all_annotation_columns_reduced,c("CLINICAL_SIGNIFICANCE","EVIDENCE_LEVEL","EVIDENCE_TYPE","EVIDENCE_DIRECTION","VARIANT_ORIGIN","CANCER_TYPE","DESCRIPTION","CITATION","THERAPEUTIC_CONTEXT","RATING","BIOMARKER_MAPPING"))
      clinical_evidence_items <- dplyr::select(clinical_evidence_items, dplyr::one_of(all_tier1_tags))
      unique_variants <- clinical_evidence_items %>% dplyr::select(SYMBOL,CONSEQUENCE,PROTEIN_CHANGE,CDS_CHANGE) %>% dplyr::distinct()
      clinical_evidence_items <- clinical_evidence_items %>% dplyr::arrange(EVIDENCE_LEVEL,RATING)
      rlogging::message(paste0(nrow(clinical_evidence_items),' clinical evidence item(s) found .. (',nrow(unique_variants),' unique variant(s)), mapping = ',mapping))
      rlogging::message('Underlying variant(s):')
      for(i in 1:nrow(unique_variants)){
        rlogging::message(paste(unique_variants[i,],collapse=" "))
      }
      all_eitems <- rbind(all_eitems, clinical_evidence_items)
    }
    else{
      rlogging::message(paste0(nrow(clinical_evidence_items),' clinical evidence item(s) found .. mapping = ',mapping))
    }

  }

  return(all_eitems)

}

#' Function that adds read support (depth, allelic fraction) for tumor and normal
#'
#' @param vcf_data_df
#'
#' @return vcf_data_df
#'
add_read_support <- function(vcf_data_df,tumor_dp_tag = '_na', tumor_af_tag = '_na', normal_dp_tag = '_na', normal_af_tag = '_na', call_conf_tag = '_na'){
  for(v in c('DP_TUMOR','AF_TUMOR','DP_NORMAL','AF_NORMAL','CALL_CONFIDENCE')){
    vcf_data_df[v] <- NA
  }
  if(tumor_dp_tag != '_na' & tumor_dp_tag %in% colnames(vcf_data_df)){
    vcf_data_df[,'DP_TUMOR'] <- as.integer(vcf_data_df[,tumor_dp_tag])
  }
  if(tumor_af_tag != '_na' & tumor_af_tag %in% colnames(vcf_data_df)){
    vcf_data_df[,'AF_TUMOR'] <- round(as.numeric(vcf_data_df[,tumor_af_tag]), digits=3)
  }
  if(normal_dp_tag != '_na' & normal_dp_tag %in% colnames(vcf_data_df)){
    vcf_data_df[,'DP_NORMAL'] <- as.integer(vcf_data_df[,normal_dp_tag])
  }
  if(normal_af_tag != '_na' & normal_af_tag %in% colnames(vcf_data_df)){
    vcf_data_df[,'AF_NORMAL'] <- round(as.numeric(vcf_data_df[,normal_af_tag]), digits=3)
  }
  if(call_conf_tag != '_na' & call_conf_tag %in% colnames(vcf_data_df)){
    vcf_data_df[,'CALL_CONFIDENCE'] <- as.character(vcf_data_df[,call_conf_tag])
  }

  return(vcf_data_df)
}


#' Function that adds SwissProt feature descriptions based on keys coming from pcgr pipeline
#'
#' @param vcf_data_df
#'
#' @return vcf_data_df
#'
add_swissprot_feature_descriptions <- function(vcf_data_df, pcgr_data){
  rlogging::message("Extending annotation descriptions related to UniprotKB/SwissProt protein features")
  swissprot_features <- pcgr_data$swissprot_features
  swissprot_features$UNIPROT_FEATURE <- paste(paste(swissprot_features$uniprot_id,swissprot_features$feature_type,sep=":"),paste(swissprot_features$aa_start,swissprot_features$aa_stop,sep="-"),sep=":")
  swissprot_features$PF <- paste(paste(swissprot_features$type_description,paste(swissprot_features$aa_start,swissprot_features$aa_stop,sep="-"),sep=":"),swissprot_features$description,sep=":")

  if("UNIPROT_FEATURE" %in% colnames(vcf_data_df) & "VAR_ID" %in% colnames(vcf_data_df)){
    feature_df <- dplyr::select(vcf_data_df,UNIPROT_FEATURE,VAR_ID) %>% dplyr::distinct()
    if(nrow(feature_df) == 0){
      vcf_data_df$PROTEIN_FEATURE <- NA
      return(vcf_data_df)
    }
    feature_df <- feature_df %>% tidyr::separate_rows(UNIPROT_FEATURE,sep="&")
    feature_df <- as.data.frame(dplyr::left_join(feature_df,dplyr::select(swissprot_features,UNIPROT_FEATURE,PF),by=c("UNIPROT_FEATURE")))
    feature_df <- as.data.frame(dplyr::group_by(feature_df, VAR_ID) %>% dplyr::summarise(PROTEIN_FEATURE = paste(PF, collapse=", ")))
    if(nrow(feature_df) > 0){
      if(nrow(feature_df[feature_df$PROTEIN_FEATURE == "NA",]) > 0){
        feature_df[feature_df$PROTEIN_FEATURE == "NA",]$PROTEIN_FEATURE <- NA
      }
      vcf_data_df <- dplyr::left_join(dplyr::select(vcf_data_df,-UNIPROT_FEATURE),feature_df, by=c("VAR_ID" = "VAR_ID"))
    }
    else{
      vcf_data_df$PROTEIN_FEATURE <- NA
    }
  }
  else{
    vcf_data_df$PROTEIN_FEATURE <- NA
  }

  return(vcf_data_df)

}

#' Function that reads a VCF from pcgr pipeline
#'
#' @param vcf_gz_file
#' @param pcgr_data Lisf ot data frames with PCGR data annotations
#' @param sample_id Sample identifier
#'
#' @return vcf_data_df
#'
get_calls <- function(vcf_gz_file, pcgr_data, sample_id = NULL, tumor_dp_tag = '_na', tumor_af_tag = '_na', normal_dp_tag = '_na', normal_af_tag = '_na', call_conf_tag = '_na'){
  if(!file.exists(vcf_gz_file) | file.size(vcf_gz_file) == 0){
    rlogging::stop(paste0("File ",vcf_gz_file," does not exist or has zero size"))
  }
  rlogging::message(paste0("Reading and parsing VEP/vcfanno-annotated VCF file - ",vcf_gz_file))
  vcf_data_vr <- VariantAnnotation::readVcfAsVRanges(vcf_gz_file,genome = "hg19")
  n_all_unfiltered_calls <- length(vcf_data_vr)
  #vcf_data_vr <- vcf_data_vr[!is.na(vcf_data_vr$GT) & !(vcf_data_vr$GT == '.'),]
  vcf_data_vr <- vcf_data_vr[VariantAnnotation::called(vcf_data_vr)]
  vcf_data_vr <- pcgrr2::postprocess_vranges_info(vcf_data_vr)
  vcf_data_df <- as.data.frame(vcf_data_vr)
  if(!is.null(sample_id) & nrow(vcf_data_df) > 0){
    vcf_data_df$VCF_SAMPLE_ID <- sample_id
  }
  rlogging::message(paste0("Number of PASS variants: ",nrow(vcf_data_df)))
  rlogging::message(paste0("Number of non-PASS variants: ",n_all_unfiltered_calls - nrow(vcf_data_df)))
  rlogging::message("Now considering PASS variants only")
  if(any(grepl(paste0("VARIANT_CLASS$"),names(vcf_data_df)))){
    n_snvs <- nrow(vcf_data_df[vcf_data_df$VARIANT_CLASS == 'SNV',])
    n_deletions <- nrow(vcf_data_df[vcf_data_df$VARIANT_CLASS == 'deletion',])
    n_insertions <- nrow(vcf_data_df[vcf_data_df$VARIANT_CLASS == 'insertion',])
    n_substitutions <- nrow(vcf_data_df[vcf_data_df$VARIANT_CLASS == 'substitution',])
    rlogging::message(paste0("Number of SNVs: ",n_snvs))
    rlogging::message(paste0("Number of deletions: ",n_deletions))
    rlogging::message(paste0("Number of insertions ",n_insertions))
    rlogging::message(paste0("Number of block substitutions ",n_substitutions))
  }
  if(nrow(vcf_data_df) == 0){
    rlogging::warning("Number of PASS variants in input VCF is 0 - no variants will be written to HTML report")
    for(col in c('GENOME_VERSION','GENOMIC_CHANGE','PROTEIN_DOMAIN','PROTEIN_FEATURE','OTHER_LITERATURE_DOCM','OTHER_DISEASE_DOCM','GENENAME','GENE_NAME','CLINVAR_TRAITS_ALL','CLINVAR','COSMIC','DBSNP','KEGG_PATHWAY','ANTINEOPLASTIC_DRUG_INTERACTIONS','DELETION_MECHANISM','CALL_CONFIDENCE')){
      vcf_data_df[col] <- character(nrow(vcf_data_df))
    }
    for(col in c('DP_TUMOR','DP_NORMAL')){
      vcf_data_df[col] <- integer(nrow(vcf_data_df))
    }
    for (col in c('AF_TUMOR','AF_NORMAL')){
      vcf_data_df[col] <- numeric(nrow(vcf_data_df))
    }
    vcf_data_df <- dplyr::rename(vcf_data_df, CHROM = seqnames, POS = start, REF = ref, ALT = alt, CONSEQUENCE = Consequence, PROTEIN_CHANGE = HGVSp_short)
    return(vcf_data_df)
  }

  vcf_data_df$GENOME_VERSION <- 'GRCh37'
  vcf_data_df <- dplyr::rename(vcf_data_df, CHROM = seqnames, POS = start, REF = ref, ALT = alt, CONSEQUENCE = Consequence, PROTEIN_CHANGE = HGVSp_short)
  vcf_data_df$GENOMIC_CHANGE <- paste(paste(paste(paste0("g.chr",vcf_data_df$CHROM),vcf_data_df$POS,sep=":"),vcf_data_df$REF,sep=":"),vcf_data_df$ALT,sep=">")
  vcf_data_df <- pcgrr2::get_deletion_mechanism(vcf_data_df)
  vcf_data_df <- pcgrr2::add_pfam_domain_links(vcf_data_df, pcgr_data = pcgr_data)
  vcf_data_df <- pcgrr2::add_swissprot_feature_descriptions(vcf_data_df, pcgr_data = pcgr_data)
  vcf_data_df <- pcgrr2::add_read_support(vcf_data_df, tumor_dp_tag = tumor_dp_tag, tumor_af_tag = tumor_af_tag, normal_dp_tag = normal_dp_tag, normal_af_tag = normal_af_tag, call_conf_tag = call_conf_tag)
  rlogging::message("Extending annotation descriptions related to Database of Curated Mutations (DoCM)")
  vcf_data_df <- dplyr::left_join(vcf_data_df, pcgr_data$docm_literature, by=c("VAR_ID"))

  gencode_xref <- dplyr::rename(pcgr_data$gene_xref, Gene = ensembl_gene_id, GENENAME = name, ENTREZ_ID = entrezgene)
  gencode_xref <- gencode_xref %>% dplyr::filter(!is.na(Gene)) %>% dplyr::select(Gene, GENENAME, ENTREZ_ID) %>% dplyr::distinct()
  gencode_xref$GENENAME <- stringr::str_replace(gencode_xref$GENENAME," \\[.{1,}$","")
  gencode_xref$ENTREZ_ID <- as.character(gencode_xref$ENTREZ_ID)
  gencode_xref <- dplyr::filter(gencode_xref, !is.na(GENENAME) & !is.na(ENTREZ_ID))

  vcf_data_df <- dplyr::left_join(vcf_data_df,gencode_xref,by=c("ENTREZ_ID","Gene"))
  vcf_data_df <- dplyr::left_join(vcf_data_df,pcgr_data$kegg_gene_pathway_links, by=c("ENTREZ_ID" = "gene_id"))
  vcf_data_df <- dplyr::rename(vcf_data_df, KEGG_PATHWAY = kegg_pathway_urls)

  tmp <- dplyr::select(pcgr_data$clinvar, CLINVAR_TRAITS_ALL, CLINVAR_MSID, var_id)
  tmp$var_id <- stringr::str_replace(tmp$var_id,"chr","")
  tmp <- dplyr::rename(tmp, VAR_ID = var_id)
  if ("CLINVAR_MSID" %in% colnames(vcf_data_df)){
    rlogging::message("Extending annotation descriptions related to ClinVar")
    vcf_data_df <- dplyr::left_join(vcf_data_df,tmp,by=c("CLINVAR_MSID","VAR_ID"))
  }

  if("EFFECT_PREDICTIONS" %in% colnames(vcf_data_df)){
    vcf_data_df$EFFECT_PREDICTIONS <- stringr::str_replace_all(vcf_data_df$EFFECT_PREDICTIONS,"\\.&|\\.$","NA&")
    vcf_data_df$EFFECT_PREDICTIONS <- stringr::str_replace_all(vcf_data_df$EFFECT_PREDICTIONS,"&$","")
    vcf_data_df$EFFECT_PREDICTIONS <- stringr::str_replace_all(vcf_data_df$EFFECT_PREDICTIONS,"&",", ")
    vcf_data_df <- pcgrr2::annotate_variant_link(vcf_data_df, vardb = 'DBNSFP')

  }
  if("ONCOSCORE" %in% colnames(vcf_data_df)){
    if(nrow(vcf_data_df[is.na(vcf_data_df$ONCOSCORE),]) > 0){
      vcf_data_df[is.na(vcf_data_df$ONCOSCORE),]$ONCOSCORE <- 0
    }
  }

  vcf_data_df <- pcgrr2::df_string_replace(vcf_data_df, strings = c('COSMIC_SITE_HISTOLOGY','INTOGEN_DRIVER_MUT','CONSEQUENCE','CANCER_CENSUS_SOMATIC','CANCER_CENSUS_GERMLINE','COSMIC_DRUG_RESISTANCE'), pattern = "&", replacement = ", ", replace_all = T)
  vcf_data_df <- pcgrr2::df_string_replace(vcf_data_df, strings = c('VEP_ALL_CONSEQUENCE','DOCM_DISEASE'), pattern = ",", replacement = ", ", replace_all = T)

  ## Add HTML links for COSMIC, DBSNP, DGIDB and CLINVAR entries
  if(!("COSMIC" %in% colnames(vcf_data_df))){
    vcf_data_df <- pcgrr2::annotate_variant_link(vcf_data_df, vardb = 'COSMIC')
    vcf_data_df <- dplyr::rename(vcf_data_df, COSMIC = COSMICLINK)
  }
  if(!("DBSNP" %in% colnames(vcf_data_df))){
    vcf_data_df <- pcgrr2::annotate_variant_link(vcf_data_df, vardb = 'DBSNP')
    vcf_data_df <- dplyr::rename(vcf_data_df, DBSNP = DBSNPLINK)
  }
  if(!("CLINVAR" %in% colnames(vcf_data_df))){
    vcf_data_df <- pcgrr2::annotate_variant_link(vcf_data_df, vardb = 'CLINVAR')
    vcf_data_df <- dplyr::rename(vcf_data_df, CLINVAR = CLINVARLINK)
  }
  if(!("ANTINEOPLASTIC_DRUG_INTERACTIONS" %in% colnames(vcf_data_df))){
    vcf_data_df <- pcgrr2::annotate_variant_link(vcf_data_df, vardb = 'DGIDB')
    vcf_data_df <- dplyr::rename(vcf_data_df, ANTINEOPLASTIC_DRUG_INTERACTIONS = DGIDBLINK)
  }
  if(!("GENE_NAME" %in% colnames(vcf_data_df))){
    vcf_data_df <- pcgrr2::annotate_variant_link(vcf_data_df, vardb = 'NCBI_GENE')
    vcf_data_df <- dplyr::rename(vcf_data_df, GENE_NAME = NCBI_GENE_LINK)
  }
  return(vcf_data_df)

}

#' Function that annotates underlying mechanism (repeat-mediated, micro-homology etc) of indels
#'
#' @param vcf_data_df Data frame with somatic variants
#'
#' @return vcf_data_df
#'

get_deletion_mechanism <- function(vcf_data_df, genome = 'hg19'){

  rlogging::message("Determination of deletion mechanisms (repeat, microhomology etc) by investigation of junction sequence")
  vcf_data_df_valid <- pcgrr2::get_valid_chromosomes(vcf_data_df, chromosome_column = 'CHROM', bsg = BSgenome.Hsapiens.UCSC.hg19)
  deletions <- vcf_data_df_valid[vcf_data_df$VARIANT_CLASS == 'deletion',]
  if(nrow(deletions) > 0){
    deletions <- dplyr::select(deletions, CHROM, POS, end, REF, ALT, VARIANT_CLASS, GENOMIC_CHANGE)
    deletions$DELETED_BASES <- substring(deletions$REF,2)
    deletions <- dplyr::filter(deletions, nchar(DELETED_BASES) > 1)
    if(nrow(deletions) > 0){
      genome_seq <- BSgenome.Hsapiens.UCSC.hg19
      seqinfo <- GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), genome = 'hg19')
      if(genome == 'hg38'){
        seqinfo <- GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38)), genome = 'hg38')
        genome_seq <- BSgenome.Hsapiens.UCSC.hg38
      }
      vcf_df_gr <- GenomicRanges::makeGRangesFromDataFrame(deletions, keep.extra.columns = T, seqinfo = seqinfo, seqnames.field = 'CHROM',start.field = 'POS', end.field = 'end', ignore.strand = T, starts.in.df.are.0based = F)

      downstream_gr <- GenomicRanges::flank(vcf_df_gr, width = nchar(mcols(vcf_df_gr)$DELETED_BASES), start=F)
      downstream_seq <- Biostrings::getSeq(genome_seq, downstream_gr)
      upstream_gr <- GenomicRanges::flank(vcf_df_gr, width = nchar(mcols(vcf_df_gr)$DELETED_BASES), start=T)
      upstream_seq <- Biostrings::getSeq(genome_seq, upstream_gr)
      df_up <- data.frame('FLANK_UPSTREAM'=toupper(unlist(strsplit(toString(upstream_seq),", "))),stringsAsFactors=F)
      df_down <- data.frame('FLANK_DOWNSTREAM'=toupper(unlist(strsplit(toString(downstream_seq),", "))),stringsAsFactors=F)
      deletions_flank_df <- cbind(deletions, df_up, df_down)
      #deletions_flank_df$stringdist_downstream_flank <- stringdist::stringdist(deletions_flank_df$DELETED_BASES, deletions_flank_df$FLANK_DOWNSTREAM)
      #deletions_flank_df$stringdist_upstream_flank <- stringdist::stringdist(deletions_flank_df$DELETED_BASES, deletions_flank_df$FLANK_UPSTREAM)
      deletions_flank_df$downstream_startswith_del <- FALSE
      if(nrow(deletions_flank_df[substring(deletions_flank_df$DELETED_BASES,1,2) == substring(deletions_flank_df$FLANK_DOWNSTREAM,1,2),]) > 0){
        deletions_flank_df[substring(deletions_flank_df$DELETED_BASES,1,2) == substring(deletions_flank_df$FLANK_DOWNSTREAM,1,2),]$downstream_startswith_del <- TRUE
      }
      deletions_flank_df$DELETION_MECHANISM <- NA
      if(nrow(deletions_flank_df[deletions_flank_df$DELETED_BASES == deletions_flank_df$FLANK_DOWNSTREAM,]) > 0){
        deletions_flank_df[deletions_flank_df$DELETED_BASES == deletions_flank_df$FLANK_DOWNSTREAM,]$DELETION_MECHANISM <- 'Repeat-mediated'
      }
      if(nrow(deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) & deletions_flank_df$downstream_startswith_del == TRUE,]) > 0){
        deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) & deletions_flank_df$downstream_startswith_del == TRUE,]$DELETION_MECHANISM <- 'Microhomology'
      }
      if(nrow(deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) & stringr::str_detect(deletions_flank_df$DELETED_BASES,'^(TTT|AAA)') & stringr::str_detect(deletions_flank_df$FLANK_UPSTREAM,'(TTT|AAA)$'),]) > 0){
        deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) & stringr::str_detect(deletions_flank_df$DELETED_BASES,'^(TTT|AAA)') & stringr::str_detect(deletions_flank_df$FLANK_UPSTREAM,'(TTT|AAA)$'),]$DELETION_MECHANISM <- 'Mononucleotide-repeat-mediated'
      }
      if(nrow(deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM),]) > 0){
        deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM),]$DELETION_MECHANISM <- 'Other'
      }

      deletions_flank_df <- dplyr::select(deletions_flank_df, GENOMIC_CHANGE, DELETION_MECHANISM)
      vcf_data_df <- dplyr::left_join(vcf_data_df, deletions_flank_df, by=c("GENOMIC_CHANGE"))

    }
  }

  return(vcf_data_df)

}


#' Function that filters a data frame with variants according to population-specific germline frequencies
#'
#' @param var_df data frame with variants
#' @param pop population ('european' or 'global')
#' @param dbquery '1KG' or 'ExAC'
#' @param min_af minimum allele frequency required for variant to be filtered
#'
#' @return var_df
#'

filter_db_germline_variants <- function(var_df, pop='european',dbquery = '1KG', min_af = 0.05){

  if(pop == 'nor'){
    if(any(grepl(paste0("^AF_NOR$"),names(var_df)))){
      var_df <- var_df %>% dplyr::filter(is.na(AF_NOR) | AF_NOR < min_af)
    }
  }
  else{
    if(pop == 'european' | pop == 'global'){
      if(dbquery == '1KG'){
        if(any(grepl(paste0("^EUR_AF_1KG$"),names(var_df)))){
          var_df <- var_df %>% dplyr::filter(is.na(EUR_AF_1KG) | EUR_AF_1KG < min_af)
        }
      }
      else{
        if(any(grepl(paste0("^NFE_AF_EXAC$"),names(var_df)))){
          var_df <- var_df %>% dplyr::filter(is.na(NFE_AF_EXAC) | NFE_AF_EXAC < min_af)
        }
      }
    }
    if(pop == 'global'){
      if(dbquery == '1KG'){
        pop_tags <- c('EAS_AF_1KG','AMR_AF_1KG','AFR_AF_1KG','SAS_AF_1KG')
        for(poptag in pop_tags){
          if(any(grepl(paste0("^",poptag,"$"),names(var_df)))){
            var_df <- var_df[is.na(var_df[poptag]) | var_df[poptag] < min_af,]
          }
        }
      }
      else{
        pop_tags <- c('SAS_AF_EXAC','EAS_AF_EXAC','AMR_AF_EXAC','AFR_AF_EXAC','FIN_AF_EXAC','OTH_AF_EXAC')
        for(poptag in pop_tags){
          if(any(grepl(paste0("^",poptag,"$"),names(var_df)))){
            var_df <- var_df[is.na(var_df[poptag]) | var_df[poptag] < min_af,]
          }
        }
      }
    }
  }
  return(var_df)
}


#' Function that assigns a category ('Rare','Common' etc) to population-specific germline frequencies
#'
#' @param var_df data frame with variants
#' @param pop_af_column population_column
#'
#' @return var_df
#'
assign_poplevel_frequency_class <- function(var_df, pop_af_column){

  if(any(grepl(paste0("^",pop_af_column,"$"),names(var_df)))){
    if(nrow(var_df[!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.05,]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.05),]$pop_common <- 'Common'
    }
    if(nrow(var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.01 & var_df[pop_af_column] < 0.05),]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.01 & var_df[pop_af_column] < 0.05),]$pop_lowfreq <- 'LowFreq'
    }
    if(nrow(var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.001 & var_df[pop_af_column] < 0.01),]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.001 & var_df[pop_af_column] < 0.01),]$pop_rare <- 'Rare'
    }
    if(nrow(var_df[!is.na(var_df[pop_af_column]) & var_df[pop_af_column] < 0.001 & var_df[pop_af_column] > 0,]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] < 0.001 & var_df[pop_af_column] > 0),]$pop_veryrare <- 'VeryRare'
    }
    if(nrow(var_df[!is.na(var_df[pop_af_column]) & var_df[pop_af_column] == 0.00,]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] == 0.00),]$pop_monomorphic <- 'Monomorphic'
    }
  }
  return(var_df)

}

#' Function that assigns a category ('Rare','Common' etc) to population-specific germline frequencies
#'
#' @param var_df data frame with variants
#' @param dbquery 1KG or ExAC
#' @param pop population
#' @param result_tag name of result column
#'
#' @return var_df
#'
assign_poplevel_frequency <- function(var_df, dbquery='1KG', pop='european', result_tag = 'FREQ_EXAC_EUROPEAN'){

  var_df$pop_monomorphic <- rep(NA,nrow(var_df))
  var_df$pop_common <- rep(NA,nrow(var_df))
  var_df$pop_rare <- rep(NA,nrow(var_df))
  var_df$pop_veryrare <- rep(NA,nrow(var_df))
  var_df$pop_lowfreq <- rep(NA,nrow(var_df))

  pop_db <- data.frame('population' = 'american','db' = '1KG', 'tag' = 'AMR_AF_1KG', stringsAsFactors = F)
  pop_db <- rbind(pop_db, data.frame('population' = 'african', 'db' = '1KG', 'tag' = 'AFR_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'global', 'db' = '1KG', 'tag' = 'GLOBAL_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'european', 'db' = '1KG', 'tag' = 'EUR_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'east_asian', 'db' = '1KG', 'tag' = 'EAS_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'south_asian', 'db' = '1KG','tag' = 'SAS_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'african', 'db' = 'ExAC', 'tag' = 'AFR_AF_EXAC'))
  pop_db <- rbind(pop_db, data.frame('population' = 'american', 'db' = 'ExAC', 'tag' = 'AMR_AF_EXAC'))
  pop_db <- rbind(pop_db, data.frame('population' = 'european', 'db' = 'ExAC', 'tag' = 'NFE_AF_EXAC'))
  pop_db <- rbind(pop_db, data.frame('population' = 'east_asian', 'db' = 'ExAC', 'tag' = 'EAS_AF_EXAC'))
  pop_db <- rbind(pop_db, data.frame('population' = 'south_asian', 'db' = 'ExAC','tag' = 'SAS_AF_EXAC'))
  pop_db <- rbind(pop_db, data.frame('population' = 'global', 'db' = 'ExAC','tag' = 'GLOBAL_AF_EXAC'))

  tags <- character()
  tags <- c(dplyr::filter(pop_db, population == pop & db == dbquery)$tag)
  for(i in 1:length(tags)){
    var_df <- assign_poplevel_frequency_class(var_df, tags[i])
  }
  var_df[result_tag] <- stringr::str_replace_all(stringr::str_replace_all(paste(var_df$pop_monomorphic,var_df$pop_rare,var_df$pop_veryrare,var_df$pop_lowfreq,var_df$pop_common,sep=","), "(,{0,}NA(,){0,}){1,}",","),"(^,)|(,$)","")
  var_df[nchar(var_df[,result_tag]) == 0,result_tag] <- NA

  var_df <- dplyr::select(var_df, -pop_monomorphic, -pop_common, -pop_rare, -pop_veryrare, -pop_lowfreq)
  return(var_df)
}

