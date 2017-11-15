#' Function that predicts MSI status based on fraction of indels among calls
#'
#' @param vcf_data_df data frame with somatic mutations/indels
#' @param simpleRepeats_gr Genomic Ranges object with sequence repeats
#' @param windowMasker_gr Genomic Ranges object with
#' @param msi_prediction_model statistical model for MSI prediction/classification
#' @param indel_frac_plot_tcga distribution plot of indel fraction for TCGA samples with known MSI status
#' @param sample_name name of sample
#' @return msi_data
#'
#'
predict_msi_status <- function(vcf_data_df, simpleRepeats_gr, windowMasker_gr, msi_prediction_model, indelFracPlot_template, sample_name = 'Test'){
  rlogging::message("Predicting microsatellite instability status")
  mutations_valid <- pcgrr::get_valid_chromosomes(vcf_data_df, chromosome_column = 'CHROM', bsg = BSgenome.Hsapiens.UCSC.hg19)
  mutations_valid <- dplyr::select(mutations_valid, CHROM,POS,end,REF,ALT,CONSEQUENCE,SYMBOL,GENOMIC_CHANGE,VARIANT_CLASS,PROTEIN_DOMAIN,GENE_NAME,PROTEIN_CHANGE,PROTEIN_FEATURE,CANCER_MUTATION_HOTSPOT,OTHER_DISEASE_DOCM,OTHER_LITERATURE_DOCM,CLINVAR,TCGA_FREQUENCY,CANCER_ASSOCIATIONS,AF_TUMOR, DP_TUMOR,AF_NORMAL,DP_NORMAL,CALL_CONFIDENCE)
  mutations_valid$start_field <- mutations_valid$POS
  mutations_valid$end_field <- NA
  mutations_valid <- dplyr::filter(mutations_valid, !is.na(VARIANT_CLASS))
  if(nrow(mutations_valid[mutations_valid$VARIANT_CLASS == 'deletion',]) > 0){
    mutations_valid[mutations_valid$VARIANT_CLASS == 'deletion',]$end_field <- mutations_valid[mutations_valid$VARIANT_CLASS == 'deletion',]$end
  }
  if(nrow(mutations_valid[mutations_valid$VARIANT_CLASS == 'substitution',]) > 0){
    mutations_valid[mutations_valid$VARIANT_CLASS == 'substitution',]$end_field <- mutations_valid[mutations_valid$VARIANT_CLASS == 'substitution',]$end
  }
  if(nrow(mutations_valid[mutations_valid$VARIANT_CLASS == 'insertion',]) > 0){
    mutations_valid[mutations_valid$VARIANT_CLASS == 'insertion',]$end_field <- mutations_valid[mutations_valid$VARIANT_CLASS == 'insertion',]$POS
  }
  if(nrow(mutations_valid[mutations_valid$VARIANT_CLASS == 'SNV',]) > 0){
    mutations_valid[mutations_valid$VARIANT_CLASS == 'SNV',]$end_field <- mutations_valid[mutations_valid$VARIANT_CLASS == 'SNV',]$POS
  }
  if(nrow(mutations_valid[mutations_valid$VARIANT_CLASS == 'sequence_alteration',]) > 0){
    mutations_valid[mutations_valid$VARIANT_CLASS == 'sequence_alteration',]$end_field <- mutations_valid[mutations_valid$VARIANT_CLASS == 'sequence_alteration',]$POS
  }

  mutations_valid <- dplyr::select(mutations_valid, -end)
  mutations_valid <- dplyr::filter(mutations_valid, !is.na(end_field) & !is.na(start_field))

  ## encode as GenomicRanges to enable intersection with repeats
  seqinfo_hg19 <- GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), genome = 'hg19')
  vcf_df_gr <- GenomicRanges::makeGRangesFromDataFrame(mutations_valid, keep.extra.columns = T, seqinfo = seqinfo_hg19, seqnames.field = 'CHROM',start.field = 'start_field', end.field = 'end_field', ignore.strand = T, starts.in.df.are.0based = F)

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

  #vcf_df_repeatAnnotated <- dplyr::left_join(mutations_valid, variants_in_repeats_unique,by="GENOMIC_CHANGE")
  #vcf_df_repeatAnnotated <- dplyr::left_join(vcf_df_repeatAnnotated, variants_in_winmask_unique,by="GENOMIC_CHANGE")

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
    dplyr::filter(SYMBOL == 'MLH1' & stringr::str_detect(CONSEQUENCE,"frameshift_variant|missense_variant|splice_donor|splice_acceptor|stop_gained|stop_lost|start_lost|inframe_deletion|inframe_insertion")) %>%
    dplyr::summarise(MLH1 = n())

  msi_stats11 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'MLH3' & stringr::str_detect(CONSEQUENCE,"frameshift_variant|missense_variant|splice_donor|splice_acceptor|stop_gained|stop_lost|start_lost|inframe_deletion|inframe_insertion")) %>%
    dplyr::summarise(MLH3 = n())

  msi_stats12 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'MSH2' & stringr::str_detect(CONSEQUENCE,"frameshift_variant|missense_variant|splice_donor|splice_acceptor|stop_gained|stop_lost|start_lost|inframe_deletion|inframe_insertion")) %>%
    dplyr::summarise(MSH2 = n())

  msi_stats13 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'MSH3' & stringr::str_detect(CONSEQUENCE,"frameshift_variant|missense_variant|splice_donor|splice_acceptor|stop_gained|stop_lost|start_lost|inframe_deletion|inframe_insertion")) %>%
    dplyr::summarise(MSH3 = n())

  msi_stats14 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'MSH6' & stringr::str_detect(CONSEQUENCE,"frameshift_variant|missense_variant|splice_donor|splice_acceptor|stop_gained|stop_lost|start_lost|inframe_deletion|inframe_insertion")) %>%
    dplyr::summarise(MSH6 = n())

  msi_stats15 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'PMS1' & stringr::str_detect(CONSEQUENCE,"frameshift_variant|missense_variant|splice_donor|splice_acceptor|stop_gained|stop_lost|start_lost|inframe_deletion|inframe_insertion")) %>%
    dplyr::summarise(PMS1 = n())

  msi_stats16 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'PMS2' & stringr::str_detect(CONSEQUENCE,"frameshift_variant|missense_variant|splice_donor|splice_acceptor|stop_gained|stop_lost|start_lost|inframe_deletion|inframe_insertion")) %>%
    dplyr::summarise(PMS2 = n())

  msi_stats17 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'POLE' & stringr::str_detect(CONSEQUENCE,"frameshift_variant|missense_variant|splice_donor|splice_acceptor|stop_gained|stop_lost|start_lost|inframe_deletion|inframe_insertion")) %>%
    dplyr::summarise(POLE = n())

  msi_stats18 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(SYMBOL == 'POLD1' & stringr::str_detect(CONSEQUENCE,"frameshift_variant|missense_variant|splice_donor|splice_acceptor|stop_gained|stop_lost|start_lost|inframe_deletion|inframe_insertion")) %>%
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
  for(stat in c('fracWinMaskIndels','fracWinMaskSNVs','fracRepeatIndels','fracRepeatIndels','fracNonRepeatIndels','fracIndels')){
    if(nrow(msi_stats[is.na(msi_stats[stat]),]) > 0){
      msi_stats[is.na(msi_stats[stat]),][stat] <- 0
    }
  }

  mmr_pol_df <- mutations_valid %>% dplyr::filter(stringr::str_detect(SYMBOL, "^(MLH1|MLH3|MSH2|MSH3|MSH6|PMS1|PMS2|POLD1|POLE)$") & stringr::str_detect(CONSEQUENCE,"frameshift_variant|missense_variant|splice_donor|splice_acceptor|stop_gained|stop_lost|start_lost|inframe_deletion|inframe_insertion"))
  mmr_pol_df <- dplyr::select(mmr_pol_df, -c(CHROM,POS,REF,ALT,start_field,end_field))
  mmr_pol_df <- dplyr::rename(mmr_pol_df, GENE = SYMBOL, DOCM_DISEASE = OTHER_DISEASE_DOCM, DOCM_LITERATURE = OTHER_LITERATURE_DOCM)
  mmr_pol_df <- mmr_pol_df %>% dplyr::select(GENE, CONSEQUENCE, PROTEIN_CHANGE, GENE_NAME, VARIANT_CLASS, PROTEIN_DOMAIN, PROTEIN_FEATURE, dplyr::everything())

  msi_predictors <- c('fracWinMaskIndels','fracWinMaskSNVs','fracRepeatIndels','fracRepeatIndels','fracNonRepeatIndels','fracIndels','MLH1','MLH3','MSH2','MSH3','MSH6','PMS1','PMS2','POLD1','POLE')
  msi_class <- predict(msi_prediction_model, dplyr::select(msi_stats,msi_predictors))
  if(msi_class == 'MSS'){
    msi_stats$predicted_class <- 'MSS (Microsatellite stable)'
  }
  else{
    msi_stats$predicted_class <- 'MSI.H (Microsatellite instability - high)'
  }
  rlogging::message(paste0("Predicted MSI status: ", msi_stats$predicted_class))
  rlogging::message(paste0("MSI - Indel fraction: ", round(msi_stats$fracNonRepeatIndels, digits = 3)))
  msi_data <- list('mmr_pol_df' = mmr_pol_df, 'msi_stats' = msi_stats, 'msi_plot' = indelFracPlot_template)
  return(msi_data)

}

