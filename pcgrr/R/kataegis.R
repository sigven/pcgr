#' Function that detects kataegis events from a data frame
#' with genomic cooordinates of mutations
#'
#' @param data data frame with somatic mutations, as produced by kataegis_input
#' @param sample_id sample identifier
#' @param build genomoe assembly build
#' @param min.mut minimum number of mutations in localized hypermutated region
#' @param max.dis maximum distance of kataegis event (basepairs)
#' @param chr column name in data that denotes chromosome
#' @param pos column name in data that denotes position
#' @param txdb transcript database (txdb)
#'
#' @return kataegis_df data frame with potential kataegis events
#'
#' @export
kataegis_detect <- function(data, sample_id = "sample_id",
                            build = "grch37", min.mut = 6,
                            max.dis = 1000,
                            chr = "chr", pos = "pos",
                            txdb = NULL) {
  pcgrr::log4r_info(paste0("Detecting possible kataegis events (clusters of C>T ",
                    "(APOBEC enzyme family) and C>T/G ",
                    "(TLS DNA polymerase) mutations"))
  assertable::assert_colnames(
    data, c(chr, pos), only_colnames = F, quiet = T)
  invisible(
    assertthat::assert_that(
      build == "grch37" | build == "grch38",
      msg = paste0("Value for argument build ('",
                   build,
                   "') not allowed")))

  chr.arm <- c(1.25e+08, 93300000, 9.1e+07, 50400000, 48400000,
              6.1e+07, 59900000, 45600000, 4.9e+07, 40200000, 53700000,
              35800000, 17900000, 17600000, 1.9e+07, 36600000, 2.4e+07,
              17200000, 26500000, 27500000, 13200000, 14700000, 60600000,
              12500000)
  if (build == "grch38") {
    chr.arm <- c(123400000, 93900000, 90900000, 5e+07, 48800000,
                59800000, 60100000, 45200000, 4.3e+07, 39800000, 53400000,
                35500000, 17700000, 17200000, 1.9e+07, 36800000, 25100000,
                18500000, 26200000, 28100000, 1.2e+07, 1.5e+07, 6.1e+07,
                10400000)
  }
  num <- dim(data)[1] - 5
  katPoint <- matrix(nrow = num, ncol = 8)
  i <- 1
  mutnum <- 1
  Cmutnum <- 0
  for (i in 1:num) {
    if (data$ref[i] %in% c("C", "G")) {
      Cmutnum <- Cmutnum + 1
    }
    if (data$dis[i + 1] <= max.dis) {
      mutnum <- mutnum + 1
    } else {
      if (mutnum >= min.mut) {
        len <- data$pos[i] - data$pos[i - mutnum + 1] + 1
        chr.n <- gsub(pattern = "chr", replacement = "", x = data$chr[i],
                     fixed = TRUE)
        chr.n <- gsub(pattern = "X", replacement = "23", x = chr.n,
                     fixed = TRUE)
        chr.n <- gsub(pattern = "Y", replacement = "24", x = chr.n,
                     fixed = TRUE)
        chr.n <- as.numeric(chr.n)
        if (data$pos[i] <= chr.arm[chr.n]) {
          arm <- paste(chr.n, "p", sep = "")
        } else if (data$pos[i - mutnum + 1] >= chr.arm[chr.n]) {
          arm <- paste(chr.n, "q", sep = "")
        } else {
          arm <- paste(chr.n, "p, ", chr.n, "q", sep = "")
        }
        katPoint[i, 1:8] <-
          c(sample_id, data$chr[i], data$pos[i - mutnum + 1],
            data$pos[i], arm, len, mutnum, round(Cmutnum / mutnum, 3))
      }
      mutnum <- 1
      Cmutnum <- 0
    }
  }
  kataegis_df <- data.frame(stats::na.omit(katPoint))
  names(kataegis_df) <- c("sample_id", "chrom", "start",
                          "end", "chrom.arm", "length", "number.mut",
                          "weight.C>X")
  if (NROW(kataegis_df) > 0) {
    for (i in 1:dim(kataegis_df)[1]) {
      if (as.numeric(as.character(kataegis_df$"weight.C>X"[i])) < 0.8) {
        kataegis_df$confidence[i] <- 0
      } else {
        chrom_i <- kataegis_df$chrom[i]
        kataegis_df$confidence[i] <-
          length(which(
            subset(
              kataegis_df,
              as.numeric(
                as.character(kataegis_df$"weight.C>X")) >= 0.8)$chrom == chrom_i
            ))
        if (kataegis_df$confidence[i] > 3) {
          kataegis_df$confidence[i] <- 3
        }
      }
    }
    kataegis_df <- kataegis_df |>
      dplyr::arrange(dplyr::desc(.data$confidence))

    # if (!is.null(txdb)) {
    #   gr <-
    #     GenomicRanges::GRanges(
    #       seqnames = S4Vectors::Rle(kataegis_df$chrom),
    #       ranges = IRanges::IRanges(start = as.numeric(as.character(kataegis_df$start)),
    #                        end = as.numeric(as.character(kataegis_df$end))))
    #   peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 3000),
    #                            TxDb = txdb, annoDb = "org.Hs.eg.db")
    #   kataegis_df$annotation <- peakAnno@anno$annotation
    #   kataegis_df$distanceToTSS <- peakAnno@anno$distanceToTSS
    #   kataegis_df$geneName <- peakAnno@anno$SYMBOL
    #   kataegis_df$geneID <- peakAnno@anno$geneId
    # }
  }
  pcgrr::log4r_info(paste(dim(kataegis_df)[1],
                          "potential kataegis events identified",
                sep = " "))
  return(kataegis_df)
}

#' Function that detects kataegis events from a data frame
#' with genomic cooordinates of mutations
#'
#' @param variant_set data frame with raw set of somatic mutations
#' @param chr column name in data that denotes chromosome
#' @param pos column name in data that denotes position
#' @param ref column name in data that denotes reference allele
#' @param alt column name in data that denotes alternate allele
#' @param build genome build (grch37 or hg38)
#' @param context_size size of neighbouring sequence context
#'
#' @export
kataegis_input <- function(variant_set, chr = "chr", pos = "pos", ref = "ref",
                           alt = "alt", build = NULL, context_size = 10) {

  invisible(assertthat::assert_that(
    is.data.frame(variant_set),
    msg = paste0("Argument 'variant_set' needs be of type data.frame")))
  invisible(assertthat::assert_that(
    build == "grch37" | build == "grch38",
    msg = paste0("Value for argument build ('", build,
                 "') not allowed, allowed reference builds are: 'grch37' or 'grch38'")))
  assertable::assert_colnames(variant_set,
                              c(chr, pos, ref, alt),
                              only_colnames = F, quiet = T)

  mut_data <- variant_set[, c(chr, pos, ref, alt)]
  names(mut_data) <- c("chr", "pos", "ref", "alt")
  mut_data <- mut_data |>
    dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1 &
                    stringr::str_detect(ref, "^(A|C|T|G)$") &
                    stringr::str_detect(alt, "^(A|C|G|T)$"))

  #context_size <- 10
  if (nrow(mut_data) >= 100) {
    bsg <- get_genome_obj(build)
    chr.lens <- as.integer(utils::head(GenomeInfoDb::seqlengths(bsg), 24))
    mut_data$build <- build
    ref_base <-  Biostrings::DNAStringSet(mut_data$ref)
    alt_base <-  Biostrings::DNAStringSet(mut_data$alt)
    conv.start <-  mut_data$pos - context_size
    conv.end <-  mut_data$pos + context_size
    context <-  Biostrings::getSeq(bsg, mut_data$chr,
                                   start = conv.start, end = conv.end)
    if (TRUE) {
      idx <-  mut_data$ref %in% c("A", "G")
      context[idx] <-  Biostrings::reverseComplement(context[idx])
      ref_base[idx] <-  Biostrings::reverseComplement(ref_base[idx])
      alt_base[idx] <-  Biostrings::reverseComplement(alt_base[idx])
    }
    mut_data$alteration <-  paste(ref_base, alt_base, sep = ">")
    mut_data$context <- context
    # Replace chr x and y with numeric value (23 and 24) for better
    # ordering
    seq <-  gsub(pattern = "chr", replacement = "", x = mut_data$chr,
               fixed = TRUE)
    seq <-  gsub(pattern = "X", replacement = "23", x = seq, fixed = TRUE)
    seq <-  gsub(pattern = "Y", replacement = "24", x = seq, fixed = TRUE)
    mut_data$seq <-  as.numeric(seq)
    mut_data <-  mut_data[order(mut_data$seq, mut_data$pos), ]
    chr.lens.sum <-  cumsum(as.numeric(chr.lens))
    chr.lens.sum <-  c(0, chr.lens.sum)
    mut_data$dis <-  c(mut_data$pos[1],
                       diff(mut_data$pos + chr.lens.sum[mut_data$seq]))
  } else {
    mut_data <- NULL
  }
  return(mut_data)
}

#' Function that generates data frame with potential kataegis events
#'
#' @param variant_set data frame with SNVs/InDels (must contain 'CHROM',
#' 'POS','REF','ALT')
#' @param sample_name name of tumor sample
#' @param build genome assembly (grch37/grch38)
#'
#' @export
generate_report_data_kataegis <- function(variant_set,
                                          sample_name = "SampleX",
                                          build = "grch37") {

  pcg_report_kataegis <-
    pcgrr::init_kataegis_content()
  if (NROW(variant_set) == 0) {
    return(pcg_report_kataegis)
  }

  pcgrr::log4r_info("------")
  pcgrr::log4r_info(
    paste0("Kataegis detection from genomic distribution of SNVs"))


  invisible(assertthat::assert_that(
    is.data.frame(variant_set),
    msg = paste0("Argument 'variant_set' needs be of type data.frame")))
  assertable::assert_colnames(
    variant_set, c("CHROM", "REF", "ALT", "POS"), only_colnames = F, quiet = T)
  invisible(assertthat::assert_that(
    build == "grch37" | build == "grch38",
    msg =
      paste0("Value for argument build ('", build,
             "') not allowed, allowed reference builds are: 'grch37' or 'grch38'")))

  chr_prefix <- FALSE
  chromosome_names <- unique(variant_set[, "CHROM"])
  for (m in chromosome_names) {
    if (startsWith(m, "chr")) {
      chr_prefix <- TRUE
    }
  }
  if (chr_prefix == F) {
    variant_set <- variant_set |>
      dplyr::mutate(CHROM = paste0("chr", .data$CHROM))
  }

  kataegis_data <- pcgrr::kataegis_input(variant_set, chr = "CHROM",
                                         pos = "POS", ref = "REF",
                                         alt = "ALT",
                                         build = build)
  if (!is.null(kataegis_data)) {

    if (nrow(kataegis_data) > 100) {
      pcg_report_kataegis[["eval"]] <- TRUE

      pcg_report_kataegis[["events"]] <-
        pcgrr::kataegis_detect(kataegis_data,
                               sample_id = sample_name,
                               build = build)
    }
  }else{
    pcgrr::log4r_info(
      paste0(
        "No or too few SNVs (< 100) found in input - skipping kataegis detection"))
    #pcg_report_kataegis[["eval"]] <- FALSE
  }
  return(pcg_report_kataegis)

}
