#' Function that retrieves relative estimates of known somatic signatures from a single tumor
#'
#' @param mut_data data frame with somatic mutations (VCF_SAMPLE_ID, CHROM, POS, REF, ALT)
#' @param sample_name sample name
#' @param normalization_method metod for normalization of context counts (deconstructSigs)
#' @param cosmic_cancertypes_aetiologies list of known signaturea and associated etiologies/cancertypes
#' @param signatures_limit max number of contributing signatures
#'
#'
signature_contributions_single_sample <- function(mut_data, sample_name, normalization_method = 'default',  cosmic_signatures_aetiologies = NULL, signatures_limit = 6, bsg = BSgenome.Hsapiens.UCSC.hg19){
  n_muts = nrow(mut_data)
  rlogging::message(paste0("Identifying weighted contributions of known mutational signatures using deconstructSigs (n = ",n_muts," SNVs)"))
  rlogging::message(paste0("deconstructSigs normalization method ('tri.counts.method'): ",normalization_method))
  sigs.input <- deconstructSigs::mut.to.sigs.input(mut.ref = mut_data, sample.id = "VCF_SAMPLE_ID",chr = "CHROM",pos = "POS", ref = "REF", alt = "ALT")
  sample_1 <- deconstructSigs::whichSignatures(tumor.ref = sigs.input, sample.id = sample_name, signatures.limit = signatures_limit, signatures.ref = signatures.cosmic,contexts.needed = T,tri.counts.method = normalization_method)

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
  signatures_cancertypes_aetiologies <- signatures_cancertypes_aetiologies[,c("Signature_ID","Weight","Cancer_types","Proposed_aetiology","Comments")]
  signatures_cancertypes_aetiologies$Trimer_normalization_method <- normalization_method

  unknown_df <- dplyr::filter(signature_contributions, signature_id == 'unknown')
  unknown_df <- data.frame(Signature_ID = unknown_df$signature_id, Weight = unknown_df$weight, stringsAsFactors = F)
  unknown_df$Weight <- round(as.numeric(unknown_df$Weight),digits=3)
  unknown_df$Cancer_types <- NA
  unknown_df$Proposed_aetiology <- NA
  unknown_df$Comments <- NA
  unknown_df$Trimer_normalization_method <- normalization_method

  signatures_cancertypes_aetiologies <- rbind(signatures_cancertypes_aetiologies, unknown_df) %>% dplyr::arrange(desc(Weight))

  return(list(whichSignatures_object = sample_1, cancertypes_aetiologies = signatures_cancertypes_aetiologies))
}
