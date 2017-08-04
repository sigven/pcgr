#' Function that retrieves relative estimates of known somatic signatures from a single tumor
#'
#' @param mut_data data frame with somatic mutations (VCF_SAMPLE_ID, CHROM, POS, REF, ALT)
#' @param sample_name sample name
#' @param normalization_method metod for normalization of context counts (deconstructSigs)
#' @param signatures_limit max number of contributing signatures
#'
#'
signature_contributions_single_sample <- function(mut_data, sample_name, normalization_method = 'default',  signatures_limit = 6, bsg = BSgenome.Hsapiens.UCSC.hg19){
  n_muts = nrow(mut_data)
  rlogging::message(paste0("Identifying weighted contributions of known mutational signatures using deconstructSigs (n = ",n_muts," SNVs)"))
  rlogging::message(paste0("deconstuctSigs normalization method ('tri.counts.method'): ",normalization_method))
  sigs.input <- deconstructSigs::mut.to.sigs.input(mut.ref = mut_data, sample.id = "VCF_SAMPLE_ID",chr = "CHROM",pos = "POS", ref = "REF", alt = "ALT")
  sample_1 <- deconstructSigs::whichSignatures(tumor.ref = sigs.input, sample.id = sample_name, signatures.limit = signatures_limit, signatures.ref = signatures.cosmic,contexts.needed = T,tri.counts.method = normalization_method)

  nonzero_signatures <- sample_1$weights[which(colSums(sample_1$weights != 0) > 0)]
  n <- 1
  signature_contributions_df <- NULL
  while(n <= ncol(nonzero_signatures)){
    df <- data.frame("sample_name" = sample_name, "signature_id" = stringr::str_replace(colnames(nonzero_signatures)[n],"ignature\\.",""), "weight" = as.numeric(nonzero_signatures[,n]))
    signature_contributions_df <- rbind(signature_contributions_df,df)
    rlogging::message(paste0("Inferred weighted contribution of ",df$signature_id,": ",round(df$weight,digits = 3)))
    n <- n + 1
  }
  signature_contributions_df <- rbind(signature_contributions_df, data.frame("sample_name" = sample_name, "signature_id" = "unknown", "weight" = sample_1$unknown))
  return(list(which_signatures_obj = sample_1, which_signatures_df = signature_contributions_df))
}
