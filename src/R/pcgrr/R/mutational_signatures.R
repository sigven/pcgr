#' Function that retrieves relative estimates of known somatic signatures from a single tumor
#'
#' @param mut_data data frame with somatic mutations (VCF_SAMPLE_ID, CHROM, POS, REF, ALT)
#' @param sample_name sample name
#' @param normalization_method metod for normalization of context counts (deconstructSigs)
#' @param cosmic_cancertypes_aetiologies list of known signaturea and associated etiologies/cancertypes
#' @param signature_limit max number of contributing signatures
#' @param skip_signatures vector of signatures that are skipped in search space, e.g. c('Signature.29','Signature.23')
#' @param signature_cutoff discard any signature contributions with weight less than this amount
#' @param bsg genome sequence object (BSgenome.Hsapiens.UCSC)
#'
#'
signature_contributions_single_sample <- function(sample_calls, sample_name, normalization_method = 'default',  cosmic_signatures_aetiologies = NULL, signature_limit = 6, skip_signatures = NULL, signature_cutoff = 0.06, bsg = BSgenome.Hsapiens.UCSC.hg19){

  n_muts <- nrow(sample_calls)
  sample_calls <- pcgrr::get_valid_chromosomes(sample_calls)
  sigs.input <- deconstructSigs::mut.to.sigs.input(mut.ref = sample_calls, sample.id = "VCF_SAMPLE_ID",chr = "CHROM",pos = "POS", ref = "REF", alt = "ALT")
  all_signatures <- paste0("Signature.",rep(1:30))
  if(!is.null(skip_signatures))
    all_signatures <- setdiff(all_signatures,skip_signatures)

  sample_1 <- deconstructSigs::whichSignatures(tumor.ref = sigs.input, sample.id = sample_name, associated = all_signatures, signature.cutoff = signature_cutoff, signatures.limit = signature_limit, signatures.ref = signatures.cosmic,contexts.needed = T,tri.counts.method = normalization_method)

  nonzero_signatures <- sample_1$weights[which(colSums(sample_1$weights != 0) > 0)]
  n <- 1
  signature_contributions <- NULL
  while(n <= ncol(nonzero_signatures)){
    df <- data.frame("sample_name" = sample_name, "signature_id" = stringr::str_replace(colnames(nonzero_signatures)[n],"ignature\\.",""), "weight" = as.numeric(nonzero_signatures[,n]))
    signature_contributions <- rbind(signature_contributions,df)
    rlogging::message(paste0("Inferred weighted contribution of ",df$signature_id,": ",round(df$weight,digits = 3)))
    n <- n + 1
  }

  signature_contributions <- rbind(signature_contributions, data.frame("sample_name" = sample_name, "signature_id" = "unknown", "weight" = sample_1$unknown))
  signature_columns <- as.numeric(stringr::str_replace(as.character(signature_contributions[signature_contributions$signature_id != 'unknown',]$signature_id),"S",""))
  weight_df <- data.frame('Signature_ID' = as.character(signature_contributions$signature_id), 'Weight' = round(as.numeric(signature_contributions$weight),digits=3), stringsAsFactors = F)
  cancertypes_aetiologies <- cosmic_signatures_aetiologies[signature_columns,]
  signatures_cancertypes_aetiologies <- dplyr::left_join(cancertypes_aetiologies,weight_df,by=c("Signature_ID")) %>% dplyr::arrange(desc(Weight))
  signatures_cancertypes_aetiologies <- signatures_cancertypes_aetiologies[,c("Signature_ID","Weight","Cancer_types","Proposed_aetiology","Comments","Keyword")]
  signatures_cancertypes_aetiologies$Trimer_normalization_method <- normalization_method
  signatures_cancertypes_aetiologies$Sample_Name <- sample_name

  if(nrow(signature_contributions[signature_contributions$signature_id == 'unknown',]) > 0){
    unknown_df <- as.data.frame(dplyr::filter(signature_contributions, signature_id == 'unknown') %>% dplyr::rename(Weight = weight))
    unknown_df$Weight <- round(as.numeric(unknown_df$Weight), digits=3)
    unknown_df <- dplyr::mutate(unknown_df, Sample_Name = as.character(sample_name), Trimer_normalization_method = as.character(normalization_method))
    unknown_df <- unknown_df %>% dplyr::mutate(Signature_ID = as.character(signature_id)) %>% dplyr::select(-c(signature_id,sample_name))
    signatures_cancertypes_aetiologies <- dplyr::bind_rows(signatures_cancertypes_aetiologies, unknown_df)
  }
  signatures_cancertypes_aetiologies <- signatures_cancertypes_aetiologies %>% dplyr::arrange(desc(Weight))
  rlogging::message('------')

  return(list(deconstructsigs_which_signatures = sample_1, cancertypes_aetiologies = signatures_cancertypes_aetiologies))
}

#' Function that generates mutational signatures data for PCGR report
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param pcgr_version PCGR software version
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param genome_seq BSgenome object
#' @param genome_assembly human genome assembly version
#'
generate_report_data_signatures <- function(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, genome_assembly){

  pcg_report_signatures <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = 'm_signature')

  rlogging::message('------')
  rlogging::message(paste0("Identifying weighted contributions of known mutational signatures using deconstructSigs"))
  rlogging::message(paste0("deconstructSigs normalization method ('tri.counts.method'): ",pcgr_config$mutational_signatures$mutsignatures_normalization))
  if(any(grepl(paste0("VARIANT_CLASS$"),names(sample_calls)))){
    if(nrow(sample_calls[!is.na(sample_calls$VARIANT_CLASS) & sample_calls$VARIANT_CLASS == 'SNV',]) >= pcgr_config[['mutational_signatures']][['mutsignatures_mutation_limit']]){
      pcg_report_signatures[['variant_set']][['all']] <- sample_calls[!is.na(sample_calls$VARIANT_CLASS) & sample_calls$VARIANT_CLASS == 'SNV',]
      pcg_report_signatures[['variant_set']][['all']] <- dplyr::filter(pcg_report_signatures[['variant_set']][['all']], CHROM != 'MT')
      pcg_report_signatures[['variant_set']][['all']]$VCF_SAMPLE_ID <- sample_name
      pcg_report_signatures[['variant_set']][['all']] <- dplyr::select(pcg_report_signatures[['variant_set']][['all']], CHROM, POS, REF, ALT, VCF_SAMPLE_ID)
      pcg_report_signatures[['eval']] <- TRUE
      mut_signature_contributions <- pcgrr::signature_contributions_single_sample(pcg_report_signatures[['variant_set']][['all']],
                                                                                  sample_name = sample_name,
                                                                                  normalization_method = pcgr_config$mutational_signatures$mutsignatures_normalization,
                                                                                  cosmic_signatures_aetiologies = pcgr_data$signatures_aetiologies,
                                                                                  signature_limit = pcgr_config$mutational_signatures$mutsignatures_signature_limit,
                                                                                  signature_cutoff = pcgr_config$mutational_signatures$mutsignatures_cutoff)
      pcg_report_signatures[['result']][['deconstructsigs_which_signatures']] <- mut_signature_contributions$deconstructsigs_which_signatures
      pcg_report_signatures[['result']][['cancertypes_aetiologies']] <- mut_signature_contributions$cancertypes_aetiologies
    }else{
      if(nrow(sample_calls[!is.na(sample_calls$VARIANT_CLASS) & sample_calls$VARIANT_CLASS == 'SNV',]) > 0){
        pcg_report_signatures[['variant_set']][['all']] <- sample_calls[!is.na(sample_calls$VARIANT_CLASS) & sample_calls$VARIANT_CLASS == 'SNV',]
        pcg_report_signatures[['variant_set']][['all']] <- dplyr::filter(pcg_report_signatures[['variant_set']][['all']], CHROM != 'MT')
        pcg_report_signatures[['variant_set']][['all']]$VCF_SAMPLE_ID <- sample_name
        pcg_report_signatures[['variant_set']][['all']] <- dplyr::select(pcg_report_signatures[['variant_set']][['all']], CHROM, POS, REF, ALT, VCF_SAMPLE_ID)
      }else{
        pcg_report_signatures[['variant_set']][['all']] <- data.frame()
      }
      rlogging::message(paste0("Too few SNVs (n = ",nrow(pcg_report_signatures[['variant_set']][['all']]),") for reconstruction of mutational signatures by deconstructSigs, limit set to ", pcgr_config[['mutational_signatures']][['mutsignatures_mutation_limit']]))
      pcg_report_signatures[['missing_data']] <- TRUE

      pcg_report_signatures[['result']][['cancertypes_aetiologies']] <- data.frame()
      pcg_report_signatures[['result']][['deconstructsigs_which_signatures']] <- NULL
    }
  }

  return(pcg_report_signatures)
}
