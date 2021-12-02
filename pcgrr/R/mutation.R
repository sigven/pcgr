#' Function that assigns one of six mutation types to a list of mutations
#'
#' @param var_df data frame with variants
#'
#' @return var_df
#'
#' @export
assign_mutation_type <- function(var_df) {

  invisible(
    assertthat::assert_that(
      is.data.frame(var_df),
      msg = "Argument 'var_df' must be a valid data.frame"))
  assertable::assert_colnames(var_df, c("VARIANT_CLASS", "REF", "ALT"),
                              only_colnames = F, quiet = T)
  var_df <- var_df %>%
    dplyr::mutate(
      MUTATION_TYPE =
        dplyr::case_when(VARIANT_CLASS == "SNV" & REF == "G" &
                           ALT == "A" ~ "C>T:G>A",
                         VARIANT_CLASS == "SNV" & REF == "C" &
                           ALT == "T" ~ "C>T:G>A",
                         VARIANT_CLASS == "SNV" & REF == "G" &
                           ALT == "C" ~ "C>G:G>C",
                         VARIANT_CLASS == "SNV" & REF == "C" &
                           ALT == "G" ~ "C>G:G>C",
                         VARIANT_CLASS == "SNV" & REF == "C" &
                           ALT == "A" ~ "C>A:G>T",
                         VARIANT_CLASS == "SNV" & REF == "G" &
                           ALT == "T" ~ "C>A:G>T",
                         VARIANT_CLASS == "SNV" & REF == "A" &
                           ALT == "G" ~ "A>G:T>C",
                         VARIANT_CLASS == "SNV" & REF == "T" &
                           ALT == "C" ~ "A>G:T>C",
                         VARIANT_CLASS == "SNV" & REF == "A" &
                           ALT == "T" ~ "A>T:T>A",
                         VARIANT_CLASS == "SNV" & REF == "T" &
                           ALT == "A" ~ "A>T:T>A",
                         VARIANT_CLASS == "SNV" & REF == "A" &
                           ALT == "C" ~ "A>C:T>G",
                         VARIANT_CLASS == "SNV" & REF == "T" &
                           ALT == "G" ~ "A>C:T>G",
                         FALSE ~ as.character(NA)))
  return(var_df)
}


#' Function that transforms a tier-structured variant data frame
#' into a MAF-like data frame (for input to 2020plus, MutSigCV)
#'
#' @param maf_df data frame with somatic mutations
#' @param genome_seq BSgenome object
#' @param seqinfo seqinfo object

#' @return maf_all
#'
#' @export
get_proper_maf_alleles <- function(maf_df, genome_seq, seqinfo) {

  maf_df_valid <-
    pcgrr::get_valid_chromosomes(maf_df,
                                 chromosome_column = "Chromosome",
                                 bsg = genome_seq)
  if ("end" %in% colnames(maf_df_valid)) {
    maf_df_valid <- dplyr::select(maf_df_valid, -.data$end)
  }

  maf_snv <- maf_df_valid %>%
    dplyr::filter(.data$Variant_Type == "SNP") %>%
    dplyr::mutate(REF = .data$Reference_Allele,
                  ALT = .data$Tumor_Seq_Allele2, POS = .data$Start_Position)

  maf_all <- maf_snv
  maf_ins <- dplyr::filter(maf_df_valid, .data$Variant_Type == "INS")
  maf_del <- dplyr::filter(maf_df_valid, .data$Variant_Type == "DEL")

  if (nrow(maf_del) > 0) {
    ## get appropriate alleles (VCF-like) of reference and alternate (DELETIONS)
    maf_del_gr <-
      GenomicRanges::makeGRangesFromDataFrame(maf_del, keep.extra.columns = T,
                                              seqinfo = seqinfo,
                                              seqnames.field = "Chromosome",
                                              start.field = "Start_Position",
                                              end.field = "End_Position",
                                              ignore.strand = T,
                                              starts.in.df.are.0based = F)

    maf_del_flank_gr <- GenomicRanges::flank(maf_del_gr, width = 1, start = T)
    maf_del_flank_seq <- Biostrings::getSeq(genome_seq, maf_del_flank_gr)
    maf_del_seq <- Biostrings::getSeq(genome_seq, maf_del_gr)
    vcf_alleles_alt <-
      data.frame(ALT =
                   toupper(unlist(strsplit(toString(maf_del_flank_seq), ", "))),
                 stringsAsFactors = F)
    vcf_alleles_ref <-
      data.frame(REF =
                   toupper(unlist(strsplit(toString(maf_del_seq), ", "))),
                 stringsAsFactors = F)
    vcf_alleles <- cbind(vcf_alleles_ref, vcf_alleles_alt)
    vcf_alleles$REF <- paste0(vcf_alleles$ALT, vcf_alleles$REF)
    maf_del <- cbind(maf_del, vcf_alleles)
    maf_del$POS <- maf_del$Start_Position - 1

    maf_all <- rbind(maf_all, maf_del)
  }

  if (nrow(maf_ins) > 0) {
    ## get appropriate alleles (VCF-like) of reference and alternate (INSERTIONS)
    maf_ins_gr <-
      GenomicRanges::makeGRangesFromDataFrame(maf_ins,
                                              keep.extra.columns = T,
                                              seqinfo = seqinfo,
                                              seqnames.field = "Chromosome",
                                              start.field = "Start_Position",
                                              end.field = "Start_Position",
                                              ignore.strand = T,
                                              starts.in.df.are.0based = F)
    maf_ins_seq <- Biostrings::getSeq(genome_seq, maf_ins_gr)
    vcf_alleles_alt <-
      data.frame(REF =
                   toupper(unlist(strsplit(toString(maf_ins_seq), ", "))),
                 stringsAsFactors = F)
    maf_ins <- cbind(maf_ins, vcf_alleles_alt)
    maf_ins$ALT <- paste0(maf_ins$REF, maf_ins$Tumor_Seq_Allele2)
    maf_ins$POS <- maf_ins$Start_Position

    maf_all <- rbind(maf_all, maf_ins)
  }


  maf_all$CHROM <- stringr::str_replace(maf_all$Chromosome, "chr", "")
  maf_all$GENOMIC_CHANGE <-
    paste(paste(paste(paste0("g.chr", maf_all$CHROM),
                      maf_all$POS, sep = ":"), maf_all$REF, sep = ":"),
          maf_all$ALT, sep = ">")
  return(maf_all)

}
