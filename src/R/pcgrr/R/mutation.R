
#' Function that annotates underlying mechanism (repeat-mediated, micro-homology etc) of indels
#'
#' @param vcf_data_df Data frame with somatic variants
#'
#' @return vcf_data_df
#'

get_deletion_mechanism <- function(vcf_data_df, genome_assembly = "hg19") {

  rlogging::message("Determination of deletion mechanisms (repeat, microhomology etc) by investigation of junction sequence")
  genome_seq <- BSgenome.Hsapiens.UCSC.hg38
  if (genome_assembly == "hg19") {
    genome_seq <- BSgenome.Hsapiens.UCSC.hg19
  }
  vcf_data_df_valid <- pcgrr::get_valid_chromosomes(vcf_data_df, chromosome_column = "CHROM", bsg = genome_seq)
  deletions <- vcf_data_df_valid[vcf_data_df_valid$VARIANT_CLASS == "deletion", ]
  if (nrow(deletions) > 0) {
    deletions <- dplyr::select(deletions, CHROM, POS, end, REF, ALT, VARIANT_CLASS, GENOMIC_CHANGE)
    deletions$DELETED_BASES <- substring(deletions$REF, 2)
    deletions <- dplyr::filter(deletions, nchar(DELETED_BASES) > 1)
    if (nrow(deletions) > 0) {
      seqinfo <- GenomeInfoDb::Seqinfo(seqnames =
                                         GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(genome_seq)),
                                       seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(genome_seq)),
                                       genome = genome_assembly)
      vcf_df_gr <- GenomicRanges::makeGRangesFromDataFrame(deletions,
                                                           keep.extra.columns = T,
                                                           seqinfo = seqinfo, seqnames.field = "CHROM",
                                                           start.field = "POS", end.field = "end",
                                                           ignore.strand = T, starts.in.df.are.0based = F)

      downstream_gr <- GenomicRanges::flank(vcf_df_gr, width = nchar(mcols(vcf_df_gr)$DELETED_BASES), start = F)
      downstream_seq <- Biostrings::getSeq(genome_seq, downstream_gr)
      upstream_gr <- GenomicRanges::flank(vcf_df_gr, width = nchar(mcols(vcf_df_gr)$DELETED_BASES), start = T)
      upstream_seq <- Biostrings::getSeq(genome_seq, upstream_gr)
      df_up <- data.frame("FLANK_UPSTREAM" = toupper(unlist(strsplit(toString(upstream_seq), ", "))), stringsAsFactors = F)
      df_down <- data.frame("FLANK_DOWNSTREAM" = toupper(unlist(strsplit(toString(downstream_seq), ", "))), stringsAsFactors = F)
      deletions_flank_df <- cbind(deletions, df_up, df_down)
      deletions_flank_df$downstream_startswith_del <- FALSE
      if (nrow(deletions_flank_df[substring(deletions_flank_df$DELETED_BASES, 1, 2) == substring(deletions_flank_df$FLANK_DOWNSTREAM, 1, 2), ]) > 0) {
        deletions_flank_df[substring(deletions_flank_df$DELETED_BASES, 1, 2) == substring(deletions_flank_df$FLANK_DOWNSTREAM, 1, 2), ]$downstream_startswith_del <- TRUE
      }
      deletions_flank_df$DELETION_MECHANISM <- NA
      if (nrow(deletions_flank_df[deletions_flank_df$DELETED_BASES == deletions_flank_df$FLANK_DOWNSTREAM, ]) > 0) {
        deletions_flank_df[deletions_flank_df$DELETED_BASES == deletions_flank_df$FLANK_DOWNSTREAM, ]$DELETION_MECHANISM <- "Repeat-mediated"
      }
      if (nrow(deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) & deletions_flank_df$downstream_startswith_del == TRUE, ]) > 0) {
        deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) & deletions_flank_df$downstream_startswith_del == TRUE, ]$DELETION_MECHANISM <- "Microhomology"
      }
      if (nrow(deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) &
                                  stringr::str_detect(deletions_flank_df$DELETED_BASES, "^(TTT|AAA)") &
                                  stringr::str_detect(deletions_flank_df$FLANK_UPSTREAM, "(TTT|AAA)$"), ]) > 0) {
        deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM) &
                             stringr::str_detect(deletions_flank_df$DELETED_BASES, "^(TTT|AAA)") &
                             stringr::str_detect(deletions_flank_df$FLANK_UPSTREAM, "(TTT|AAA)$"), ]$DELETION_MECHANISM <- "Mononucleotide-repeat-mediated"
      }
      if (nrow(deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM), ]) > 0) {
        deletions_flank_df[is.na(deletions_flank_df$DELETION_MECHANISM), ]$DELETION_MECHANISM <- "Other"
      }

      deletions_flank_df <- dplyr::select(deletions_flank_df, GENOMIC_CHANGE, DELETION_MECHANISM)
      vcf_data_df <- dplyr::left_join(vcf_data_df, deletions_flank_df, by = c("GENOMIC_CHANGE"))

    }
  }

  return(vcf_data_df)

}


#' Function that assigns one of six mutation types to a list of mutations
#'
#' @param var_df data frame with variants
#'
#' @return var_df
#'
assign_mutation_type <- function(var_df) {

  invisible(assertthat::assert_that(is.data.frame(var_df), msg = "Argument 'var_df' must be a valid data.frame"))
  assertable::assert_colnames(var_df, c("VARIANT_CLASS", "REF", "ALT"), only_colnames = F, quiet = T)
  var_df <- var_df %>%
    dplyr::mutate(MUTATION_TYPE = dplyr::case_when(VARIANT_CLASS == "SNV" & REF == "G" & ALT == "A" ~ "C>T:G>A",
                                                   VARIANT_CLASS == "SNV" & REF == "C" & ALT == "T" ~ "C>T:G>A",
                                                   VARIANT_CLASS == "SNV" & REF == "G" & ALT == "C" ~ "C>G:G>C",
                                                   VARIANT_CLASS == "SNV" & REF == "C" & ALT == "G" ~ "C>G:G>C",
                                                   VARIANT_CLASS == "SNV" & REF == "C" & ALT == "A" ~ "C>A:G>T",
                                                   VARIANT_CLASS == "SNV" & REF == "G" & ALT == "T" ~ "C>A:G>T",
                                                   VARIANT_CLASS == "SNV" & REF == "A" & ALT == "G" ~ "A>G:T>C",
                                                   VARIANT_CLASS == "SNV" & REF == "T" & ALT == "C" ~ "A>G:T>C",
                                                   VARIANT_CLASS == "SNV" & REF == "A" & ALT == "T" ~ "A>T:T>A",
                                                   VARIANT_CLASS == "SNV" & REF == "T" & ALT == "A" ~ "A>T:T>A",
                                                   VARIANT_CLASS == "SNV" & REF == "A" & ALT == "C" ~ "A>C:T>G",
                                                   VARIANT_CLASS == "SNV" & REF == "T" & ALT == "G" ~ "A>C:T>G",
                                                   FALSE ~ as.character(NA)))
  return(var_df)
}


#' Function that transforms a tier-structured variant data frame into a MAF-like data frame (for input to 2020plus, MutSigCV)
#'
#' @param tier_df data frame with somatic mutations
#' @return maf_df
#'
#'
get_proper_maf_alleles <- function(maf_df, genome_seq, seqinfo) {

  maf_df_valid <- pcgrr::get_valid_chromosomes(maf_df,
                                               chromosome_column = "Chromosome",
                                               bsg = genome_seq)
  if ("end" %in% colnames(maf_df_valid)) {
    maf_df_valid <- dplyr::select(maf_df_valid, -end)
  }

  maf_SNV <- maf_df_valid %>%
    dplyr::filter(Variant_Type == "SNP") %>%
    dplyr::mutate(REF = Reference_Allele, ALT = Tumor_Seq_Allele2, POS = Start_Position)

  maf_ALL <- maf_SNV
  maf_INS <- dplyr::filter(maf_df_valid, Variant_Type == "INS")
  maf_DEL <- dplyr::filter(maf_df_valid, Variant_Type == "DEL")

  if (nrow(maf_DEL) > 0) {
    ## get appropriate alleles (VCF-like) of reference and alternate (DELETIONS)
    maf_del_gr <- GenomicRanges::makeGRangesFromDataFrame(maf_DEL, keep.extra.columns = T,
                                                          seqinfo = seqinfo,
                                                          seqnames.field = "Chromosome",
                                                          start.field = "Start_Position",
                                                          end.field = "End_Position",
                                                          ignore.strand = T,
                                                          starts.in.df.are.0based = F)

    maf_del_flank_gr <- GenomicRanges::flank(maf_del_gr, width = 1, start = T)
    maf_del_flank_seq <- Biostrings::getSeq(genome_seq, maf_del_flank_gr)
    maf_del_seq <- Biostrings::getSeq(genome_seq, maf_del_gr)
    vcf_alleles_alt <- data.frame(ALT = toupper(unlist(strsplit(toString(maf_del_flank_seq), ", "))), stringsAsFactors = F)
    vcf_alleles_ref <- data.frame(REF = toupper(unlist(strsplit(toString(maf_del_seq), ", "))), stringsAsFactors = F)
    vcf_alleles <- cbind(vcf_alleles_ref, vcf_alleles_alt)
    vcf_alleles$REF <- paste0(vcf_alleles$ALT, vcf_alleles$REF)
    maf_DEL <- cbind(maf_DEL, vcf_alleles)
    maf_DEL$POS <- maf_DEL$Start_Position - 1

    maf_ALL <- rbind(maf_ALL, maf_DEL)
  }

  if (nrow(maf_INS) > 0) {
    ## get appropriate alleles (VCF-like) of reference and alternate (INSERTIONS)
    maf_ins_gr <- GenomicRanges::makeGRangesFromDataFrame(maf_INS,
                                                          keep.extra.columns = T,
                                                          seqinfo = seqinfo,
                                                          seqnames.field = "Chromosome",
                                                          start.field = "Start_Position",
                                                          end.field = "Start_Position",
                                                          ignore.strand = T,
                                                          starts.in.df.are.0based = F)
    maf_ins_seq <- Biostrings::getSeq(genome_seq, maf_ins_gr)
    vcf_alleles_alt <- data.frame(REF = toupper(unlist(strsplit(toString(maf_ins_seq), ", "))), stringsAsFactors = F)
    maf_INS <- cbind(maf_INS, vcf_alleles_alt)
    maf_INS$ALT <- paste0(maf_INS$REF, maf_INS$Tumor_Seq_Allele2)
    maf_INS$POS <- maf_INS$Start_Position

    maf_ALL <- rbind(maf_ALL, maf_INS)
  }


  maf_ALL$CHROM <- stringr::str_replace(maf_ALL$Chromosome, "chr", "")
  maf_ALL$GENOMIC_CHANGE <- paste(paste(paste(paste0("g.chr", maf_ALL$CHROM), maf_ALL$POS, sep = ":"), maf_ALL$REF, sep = ":"), maf_ALL$ALT, sep = ">")
  return(maf_ALL)

}


dna_kmer_distribution <- function() {

  trimer_data <- data.frame(context = character(), k = integer(), combined_context = character(), stringsAsFactors = F)
  bases <- c("A", "C", "G", "T")
  for (b in bases) {
    for (c in bases) {
      for (d in bases) {
        triplet <- Biostrings::DNAStringSet(paste0(b, c, d))
        rev_triplet <- Biostrings::reverseComplement(triplet)
        first <- toString(triplet)
        sec <- toString(rev_triplet)
        if (sec <= first) {
          tmp <- first
          first <- sec
          sec <- tmp
        }
        combined_context <- paste(first, "/", sec, sep = "")
        trimer_data <- rbind(trimer_data, data.frame(context = as.character(triplet), k = 3,
                                                     combined_context = as.character(combined_context),
                                                     stringsAsFactors = F))
      }
    }
  }

  pentamer_data <- data.frame(context = character(), k = integer(),
                              combined_context = character(),
                              stringsAsFactors = F)
  for (b in bases) {
    for (c in bases) {
      for (d in bases) {
        for (e in bases) {
          for (f in bases) {
            pentamer <- Biostrings::DNAStringSet(paste0(b, c, d, e, f))
            rev_pentamer <- Biostrings::reverseComplement(pentamer)

            first <- toString(pentamer)
            sec <- toString(rev_pentamer)
            if (sec <= first) {
              tmp <- first
              first <- sec
              sec <- tmp
            }

            combined_context <- paste(first, "/", sec, sep = "")
            pentamer_data <- rbind(pentamer_data, data.frame(context = as.character(pentamer),
                                                             k = 5, combined_context = as.character(combined_context),
                                                             stringsAsFactors = F))
          }
        }
      }
    }
  }

  all_kmer_data <- rbind(trimer_data, pentamer_data)
  return(all_kmer_data)

}

#' A function that retrieves the nucleotide context of SNVs in the human genome
#'
#' @param vr VRanges object with SNVs
#' @param ref Reference genome (BSGenome.Hsapiens.UCSC.hg19)
#' @param k Size of odd k-mer at SNV (3, 5 etc)
#' @return df_context Data frame with four variables: variant_id, original_alteration (C>A etc.), combined_context (ACA/TGT etc.), combined_alteration (C>A/G>T)
#'
#' @export
#'

mutation_context <- function(vr, ref, k = 3, kmer_set = NULL) {

  if (k %% 2 != 1)
    stop("'k' must be odd.")
  mid <- (k + 1) / 2

  ###NEW
  vr_snvs <- vr[width(VariantAnnotation::ref(vr)) == 1 & width(VariantAnnotation::alt(vr)) == 1, ]
  mcols(vr_snvs)$var_id <- paste(as.character(seqnames(vr_snvs)),
                                 start(vr_snvs),
                                 VariantAnnotation::ref(vr_snvs),
                                 VariantAnnotation::alt(vr_snvs), sep = "_")

  gr <- GenomicRanges::granges(vr_snvs)
  s <- strand(gr)
  if (any(s == "*"))
    stop("The strand must be explicit, in order to read the correct strand.")

  ranges <- GenomicRanges::resize(gr, k, fix = "center")
  context <- Biostrings::getSeq(ref, ranges)
  original_con <- context


  ### GET COMBINED CONTEXT STRING, ORDER BY
  ### CONTEXT
  f_context <- data.frame(context = unlist(strsplit(toString(context), ", ")), stringsAsFactors = F)
  kmer_contexts <- kmer_set %>% dplyr::filter(k == k)
  combined_context <- dplyr::inner_join(f_context, kmer_contexts, by = c("context"))$combined_context

  ref_base <- Biostrings::DNAStringSet(VariantAnnotation::ref(vr_snvs))
  alt_base <- Biostrings::DNAStringSet(VariantAnnotation::alt(vr_snvs))

  unique_alterations_combined <- rep(c("C>T:G>A", "A>G:T>C", "A>T:T>A", "C>G:G>C", "A>C:T>G", "C>A:G>T"), 2)
  unique_alterations_single <- c("C>T", "A>G", "A>T", "C>G", "A>C", "C>A", "G>A", "T>C", "T>A", "G>C", "T>G", "G>T")
  unique_alterations <- data.frame(alteration = unique_alterations_single,
                                   combined_alteration = unique_alterations_combined,
                                   stringsAsFactors = F)

  alterations <- data.frame(alteration = paste(unlist(strsplit(toString(ref_base), ", ")),
                                               rep(">", length(vr_snvs)),
                                               unlist(strsplit(toString(alt_base), ", ")),
                                               sep = ""), stringsAsFactors = F)
  combined_alteration <- dplyr::inner_join(alterations, unique_alterations,
                                           by = c("alteration"))$combined_alteration
  ## CHCK: assign individually ?
  vr_snvs$original_alteration <- alterations$alteration
  vr_snvs$combined_context <- combined_context
  vr_snvs$combined_alteration <- combined_alteration

  df_context <- unique(dplyr::select(as(vr_snvs, "data.frame"), var_id, original_alteration,
                                     combined_context, combined_alteration))
  df_context$var_id <- as.character(df_context$var_id)
  df_context$transi_transv <- rep("NA", nrow(df_context))
  if (nrow(df_context[df_context$combined_alteration == "C>T:G>A" | df_context$combined_alteration == "A>G:T>C", ]) > 0) {
    df_context[df_context$combined_alteration == "C>T:G>A" | df_context$combined_alteration == "A>G:T>C", ]$transi_transv <- "transition"
  }
  if (nrow(df_context[df_context$combined_alteration != "C>T:G>A" & df_context$combined_alteration != "A>G:T>C", ]) > 0) {
    df_context[df_context$combined_alteration != "C>T:G>A" & df_context$combined_alteration != "A>G:T>C", ]$transi_transv <- "transversion"
  }

  return(df_context)
}
