#' Function that predicts MSI status based on fraction of indels among calls
#'
#' @param vcf_data_df data frame with somatic mutations/indels
#' @param simpleRepeats_gr Genomic Ranges object with sequence repeats
#' @param windowMasker_gr Genomic Ranges object with
#' @param msi_prediction_model statistical model for MSI prediction/classification
#' @param msi_prediction_dataset underlying dataset from TCGA used for development of statistical classifier
#' @param target_size_mb size of targeted genomic region (coding)
#' @param bsg BSgenome data object in Biostrings objects (e.g. BSgenome.Hsapiens.UCSC.hg19)
#' @param genome_assembly hg19/hg38
#' @param sample_name name of sample
#' @return msi_data
#'
#'
predict_msi_status <- function(vcf_data_df, simpleRepeats_gr, windowMasker_gr, msi_prediction_model, msi_prediction_dataset, target_size_mb, bsg = BSgenome.Hsapiens.UCSC.hg19, genome_assembly = 'hg19',sample_name = 'Test'){

  mutations_valid <- pcgrr::get_valid_chromosomes(vcf_data_df, chromosome_column = 'CHROM', bsg = bsg)
  mutations_valid <- dplyr::select(mutations_valid, CHROM,POS,REF,ALT,CONSEQUENCE,SYMBOL,GENOMIC_CHANGE,VARIANT_CLASS,PROTEIN_DOMAIN,GENE_NAME,PROTEIN_CHANGE,PROTEIN_FEATURE,MUTATION_HOTSPOT,DOCM_DISEASE,DOCM_LITERATURE,CLINVAR,TCGA_FREQUENCY,CANCER_ASSOCIATIONS,AF_TUMOR, DP_TUMOR,AF_NORMAL,DP_NORMAL,CALL_CONFIDENCE)
  seqinfo <- GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(bsg)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(bsg)), genome = genome_assembly)
  vcf_df_gr <- GenomicRanges::makeGRangesFromDataFrame(mutations_valid, keep.extra.columns = T, seqinfo = seqinfo, seqnames.field = 'CHROM',start.field = 'POS', end.field = 'POS', ignore.strand = T, starts.in.df.are.0based = F)

  variant_repeat_hits <- GenomicRanges::findOverlaps(vcf_df_gr, simpleRepeats_gr, type="any", select="all", ignore.strand = T)
  ranges <- simpleRepeats_gr[subjectHits(variant_repeat_hits)]
  mcols(ranges) <- c(mcols(ranges),mcols(vcf_df_gr[queryHits(variant_repeat_hits)]))
  variants_in_repeats <- as.data.frame(mcols(ranges))

  vcf_df_repeatAnnotated <- mutations_valid
  vcf_df_repeatAnnotated$repeatStatus <- NA

  ## a single variant may intersect overlapping repeats; make unique repeat variants
  variants_in_repeats_unique <- dplyr::select(variants_in_repeats, GENOMIC_CHANGE) %>% dplyr::distinct()
  if(nrow(variants_in_repeats_unique) > 0){
    variants_in_repeats_unique$repeatStatus <- 'simpleRepeat'
    vcf_df_repeatAnnotated <- dplyr::left_join(mutations_valid, variants_in_repeats_unique,by="GENOMIC_CHANGE")
  }

  variant_winmask_hits <- GenomicRanges::findOverlaps(vcf_df_gr, windowMasker_gr, type="any", select="all", ignore.strand = T)
  ranges <- windowMasker_gr[subjectHits(variant_winmask_hits)]
  mcols(ranges) <- c(mcols(ranges),mcols(vcf_df_gr[queryHits(variant_winmask_hits)]))
  variants_in_winmask <- as.data.frame(mcols(ranges))

  ## a single variant may intersect overlapping repeats; make unique repeat variants per patient tumor
  variants_in_winmask_unique <- dplyr::select(variants_in_winmask, GENOMIC_CHANGE) %>% dplyr::distinct()
  if(nrow(variants_in_winmask_unique) > 0){
    variants_in_winmask_unique$winMaskStatus <- 'winMaskDust'
    vcf_df_repeatAnnotated <- dplyr::left_join(vcf_df_repeatAnnotated, variants_in_winmask_unique,by="GENOMIC_CHANGE")
  }
  else{
    vcf_df_repeatAnnotated$winMaskStatus <- NA
  }

  msi_stats <- data.frame('sample_name' = sample_name, stringsAsFactors = F)

  msi_stats1 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(repeatStatus) & (VARIANT_CLASS == 'insertion' | VARIANT_CLASS == 'deletion')) %>%
    dplyr::summarise(repeat_indels = n())

  msi_stats2 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(repeatStatus) & VARIANT_CLASS == 'SNV') %>%
    dplyr::summarise(repeat_SNVs = n())

  msi_stats3 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(repeatStatus)) %>%
    dplyr::summarise(repeat_indelSNVs = n())

  winmask_indels <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(winMaskStatus) & (VARIANT_CLASS == 'insertion' | VARIANT_CLASS == 'deletion')) %>%
    dplyr::summarise(winmask_indels = n())

  winmask_snvs <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(winMaskStatus) & VARIANT_CLASS == 'SNV') %>%
    dplyr::summarise(winmask_SNVs = n())

  winmask_tot <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(winMaskStatus)) %>%
    dplyr::summarise(winmask_indelSNVs = n())

  msi_stats4 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(is.na(repeatStatus) & (VARIANT_CLASS == 'insertion' | VARIANT_CLASS == 'deletion')) %>%
    dplyr::summarise(nonRepeat_indels = n())

  msi_stats5 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(is.na(repeatStatus) & VARIANT_CLASS == 'SNV') %>%
    dplyr::summarise(nonRepeat_SNVs = n())

  msi_stats6 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(is.na(repeatStatus)) %>%
    dplyr::summarise(nonRepeat_indelSNVs = n())

  msi_stats7 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(VARIANT_CLASS == 'insertion' | VARIANT_CLASS == 'deletion') %>%
    dplyr::summarise(indels = n())

  msi_stats8 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(VARIANT_CLASS == 'SNV') %>%
    dplyr::summarise(SNVs = n())

  msi_stats9<- vcf_df_repeatAnnotated %>%
    dplyr::summarise(indelSNVs = n())

  msi_stats10 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'MLH1' & stringr::str_detect(CONSEQUENCE,"frameshift_|missense_|splice_|stop_|inframe_")) %>%
    dplyr::summarise(MLH1 = n())

  msi_stats11 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'MLH3' & stringr::str_detect(CONSEQUENCE,"frameshift_|missense_|splice_|stop_|inframe_")) %>%
    dplyr::summarise(MLH3 = n())

  msi_stats12 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'MSH2' & stringr::str_detect(CONSEQUENCE,"frameshift_|missense_|splice_|stop_|inframe_")) %>%
    dplyr::summarise(MSH2 = n())

  msi_stats13 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'MSH3' & stringr::str_detect(CONSEQUENCE,"frameshift_|missense_|splice_|stop_|inframe_")) %>%
    dplyr::summarise(MSH3 = n())

  msi_stats14 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'MSH6' & stringr::str_detect(CONSEQUENCE,"frameshift_|missense_|splice_|stop_|inframe_")) %>%
    dplyr::summarise(MSH6 = n())

  msi_stats15 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'PMS1' & stringr::str_detect(CONSEQUENCE,"frameshift_|missense_|splice_|stop_|inframe_")) %>%
    dplyr::summarise(PMS1 = n())

  msi_stats16 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'PMS2' & stringr::str_detect(CONSEQUENCE,"frameshift_|missense_|splice_|stop_|inframe_")) %>%
    dplyr::summarise(PMS2 = n())

  msi_stats17 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'POLE' & stringr::str_detect(CONSEQUENCE,"frameshift_|missense_|splice_|stop_|inframe_")) %>%
    dplyr::summarise(POLE = n())

  msi_stats18 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'POLD1' & stringr::str_detect(CONSEQUENCE,"frameshift_|missense_|splice_|stop_|inframe_")) %>%
    dplyr::summarise(POLD1 = n())


  msi_stats1$sample_name <- sample_name
  msi_stats2$sample_name <- sample_name
  msi_stats3$sample_name <- sample_name
  msi_stats4$sample_name <- sample_name
  msi_stats5$sample_name <- sample_name
  msi_stats6$sample_name <- sample_name
  msi_stats7$sample_name <- sample_name
  msi_stats8$sample_name <- sample_name
  msi_stats9$sample_name <- sample_name
  msi_stats10$sample_name <- sample_name
  msi_stats11$sample_name <- sample_name
  msi_stats12$sample_name <- sample_name
  msi_stats13$sample_name <- sample_name
  msi_stats14$sample_name <- sample_name
  msi_stats15$sample_name <- sample_name
  msi_stats16$sample_name <- sample_name
  msi_stats17$sample_name <- sample_name
  msi_stats18$sample_name <- sample_name
  winmask_indels$sample_name <- sample_name
  winmask_snvs$sample_name <- sample_name
  winmask_tot$sample_name <- sample_name

  msi_stats <- dplyr::left_join(msi_stats,msi_stats1,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats2,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats3,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats4,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats5,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats6,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats7,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats8,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats9,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats10,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats11,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats12,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats13,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats14,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats15,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats16,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats17,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,msi_stats18,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,winmask_tot,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,winmask_snvs,by="sample_name")
  msi_stats <- dplyr::left_join(msi_stats,winmask_indels,by="sample_name")

  msi_stats$fracWinMaskIndels <- msi_stats$winmask_indels / msi_stats$indels
  msi_stats$fracWinMaskSNVs <- msi_stats$winmask_SNVs / msi_stats$SNVs
  msi_stats$fracRepeatIndels <- msi_stats$repeat_indels / msi_stats$repeat_indelSNVs
  msi_stats$fracNonRepeatIndels <- msi_stats$nonRepeat_indels / msi_stats$nonRepeat_indelSNVs
  msi_stats$fracIndels <- msi_stats$indels / msi_stats$indelSNVs
  msi_stats$tmb <- as.numeric(msi_stats$indelSNVs) / target_size_mb
  msi_stats$tmb_indel <- as.numeric(msi_stats$indels) / target_size_mb
  msi_stats$tmb_snv <- as.numeric(msi_stats$SNVs) / target_size_mb
  for(stat in c('fracWinMaskIndels','fracWinMaskSNVs','fracRepeatIndels','fracRepeatIndels','fracNonRepeatIndels','fracIndels','tmb','tmb_snv','tmb_indel')){
    if(nrow(msi_stats[is.na(msi_stats[stat]),]) > 0){
      msi_stats[is.na(msi_stats[stat]),][stat] <- 0
    }
  }

  mmr_pol_df <- mutations_valid %>% dplyr::filter(stringr::str_detect(SYMBOL, "^(MLH1|MLH3|MSH2|MSH3|MSH6|PMS1|PMS2|POLD1|POLE)$") & stringr::str_detect(CONSEQUENCE,"frameshift_|missense_|splice_|stop_|inframe_"))
  #mmr_pol_df <- dplyr::select(mmr_pol_df, -c(CHROM,POS,REF,ALT,start_field,end_field))
  mmr_pol_df <- dplyr::select(mmr_pol_df, -c(CHROM,POS,REF,ALT))
  mmr_pol_df <- dplyr::rename(mmr_pol_df, GENE = SYMBOL)
  mmr_pol_df <- mmr_pol_df %>% dplyr::select(GENE, CONSEQUENCE, PROTEIN_CHANGE, GENE_NAME, VARIANT_CLASS, PROTEIN_DOMAIN, PROTEIN_FEATURE, dplyr::everything())

  msi_predictors <- c('fracWinMaskIndels','fracWinMaskSNVs','fracRepeatIndels','fracNonRepeatIndels','fracIndels','MLH1','MLH3','MSH2','MSH3','MSH6','PMS1','PMS2','POLD1','POLE','tmb','tmb_indel','tmb_snv')
  msi_class <- predict(msi_prediction_model, dplyr::select(msi_stats,msi_predictors))
  if(msi_class == 'MSS'){
    msi_stats$predicted_class <- 'MSS (Microsatellite stable)'
    msi_stats$vb <- 'MSS'
  }
  else{
    msi_stats$predicted_class <- 'MSI.H (Microsatellite instability - high)'
    msi_stats$vb <- 'MSI - High'
  }
  rlogging::message(paste0("Predicted MSI status: ", msi_stats$predicted_class))
  rlogging::message(paste0("MSI - Indel fraction: ", round(msi_stats$fracNonRepeatIndels, digits = 3)))
  #rlogging::message('------')
  msi_data <- list('mmr_pol_variants' = mmr_pol_df, 'msi_stats' = msi_stats, 'tcga_dataset' = msi_prediction_dataset)

  return(msi_data)

}

#' Function that generates MSI prediction data for PCGR report
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param pcgr_version PCGR software version
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param genome_seq BSgenome object
#' @param genome_assembly human genome assembly version
#'
generate_report_data_msi <- function(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, genome_assembly){

  pcg_report_msi <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = 'msi')
  rlogging::message('------')
  rlogging::message("Predicting microsatellite instability status")

  #msi_sample_calls <- sample_calls %>% dplyr::filter(stringr::str_detect(CONSEQUENCE,"^(frameshift_|missense_|splice_|synonymous_|stop_|start_lost|inframe_)"))
  msi_sample_calls <- sample_calls %>% dplyr::filter(CODING_STATUS == "coding")
  rlogging::message(paste0("n = ",nrow(msi_sample_calls)," coding variants used for MSI prediction"))
  if(nrow(msi_sample_calls) > 30){
    pcg_report_msi[['prediction']] <- pcgrr::predict_msi_status(msi_sample_calls, simpleRepeats_gr = pcgr_data$simpleRepeats_gr, windowMasker_gr = pcgr_data$windowMasker_gr, msi_prediction_model = pcgr_data$msi_model$model, msi_prediction_dataset = pcgr_data$msi_model$tcga_dataset, target_size_mb = pcgr_config$mutational_burden$target_size_mb, bsg = genome_seq, genome_assembly = genome_assembly, sample_name = sample_name)
    pcg_report_msi[['eval']] <- TRUE
  }
  else{
    rlogging::message(paste0("Too few variants (n < 30) for robust MSI prediction:"))
    pcg_report_msi[['missing_data']] <- TRUE
  }

  return(pcg_report_msi)
}

#' Function that plots the indel fraction for a given sample and contrasts this with the distribution for MSI-H/MSS samples from TCGA
#'
#' @param tcga_msi_dataset underlying dataset from TCGA used for development of statistical classifier
#' @param indel_fraction fraction of indels of all mutations (SNVs + indels)
#' @return p
#'
#'

msi_indel_fraction_plot <- function(tcga_msi_dataset, indel_fraction){

  p <- ggplot2::ggplot(data = tcga_msi_dataset) +
    ggplot2::geom_histogram(mapping = ggplot2::aes(x = fracIndels, color = MSI_status, fill = MSI_status), position = "dodge", binwidth = 0.01) +
    ggplot2::ylab("Number of TCGA samples") +
    ggplot2::theme(
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(colour = "grey50"),
      panel.grid.minor = element_line(colour = "grey50"),
      panel.ontop = F
    ) +
    ggplot2::scale_color_brewer(palette='Dark2') +
    ggplot2::scale_fill_brewer(palette='Dark2') +
    ggplot2::xlab("InDel fraction among somatic calls") +
    ggplot2::theme(plot.title = ggplot2::element_text(family="Helvetica",size=16,hjust=0.5, face="bold"),
                   axis.text.x= ggplot2::element_text(family="Helvetica",size=14,face="bold"),
                   axis.title.x = ggplot2::element_text(family="Helvetica",size=14),
                   legend.title = ggplot2::element_blank(),
                   legend.text=ggplot2::element_text(family="Helvetica",size=14),
                   axis.text.y=ggplot2::element_text(family="Helvetica",size=14),
                   axis.title.y=ggplot2::element_text(family="Helvetica",size=14,vjust=1.5),
                   plot.margin = (grid::unit(c(0.5, 2, 2, 0.5), "cm"))) +
    ggplot2::geom_vline(xintercept=as.numeric(indel_fraction), size=1.4,linetype='dashed')

  return(p)

}


