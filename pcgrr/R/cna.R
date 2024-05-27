#' Function that prepares somatic copy number and mutation data for CNAqc plot
#'
#' @param callset_snv Somatic SNV callset list object - PCGR
#' @param callset_cna Somatic CNA callset list object - PCGR
#' @param settings PCGR run/configuration settings
#'
#' @export
make_cnaqc_object <- function(
    callset_snv = NULL,
    callset_cna = NULL,
    settings = NULL) {

  if (!("variant" %in% names(callset_cna) &
       "variant" %in% names(callset_snv))) {
    pcgrr::log4r_fatal("'callset_snv' or 'callset_cna' lacks variant set")
  }

  assertthat::assert_that(
    is.data.frame(callset_cna[['variant']]))
  assertthat::assert_that(
    is.data.frame(callset_snv[['variant']]))

  assertable::assert_colnames(
    callset_snv[['variant']],
    c("CHROM","POS",
      "REF","ALT",
      "DP_TUMOR",
      "VAF_TUMOR",
      "SYMBOL",
      "CONSEQUENCE",
      "ONCOGENICITY",
      "HGVSP"),
    only_colnames = F,
    quiet = T)

  assertable::assert_colnames(
    callset_cna[['variant']],
    c("CHROM",
      "SEGMENT_START",
      "SEGMENT_END",
      "CN_MAJOR",
      "CN_MINOR"),
    only_colnames = F,
    quiet = T)

  mutations <- data.frame()

  if(NROW(callset_snv$variant) > 0){

    mutations <- callset_snv$variant |>
      dplyr::mutate(
        from = .data$POS,
        to = .data$POS,
        ref = .data$REF,
        alt = .data$ALT,
        chr = paste0('chr', .data$CHROM),
        FILTER = "PASS",
        DP = .data$DP_TUMOR,
        VAF = .data$VAF_TUMOR,
        NV = as.integer(
          round(
            .data$DP_TUMOR *
              .data$VAF_TUMOR,
            digits = 0)),
        GENE = .data$SYMBOL,
        is_driver = FALSE,
        driver_label = paste0(
          .data$SYMBOL," ",
          .data$HGVSP)
      ) |>
      dplyr::mutate(
        to = dplyr::if_else(
          .data$VARIANT_CLASS == "insertion",
          .data$POS + nchar(.data$ALT) - 1,
          as.integer(.data$to)
        )
      ) |>
      dplyr::mutate(is_driver = dplyr::if_else(
        stringr::str_detect(
          .data$ONCOGENICITY, "Oncogenic") &
          .data$ONCOGENE == TRUE &
          !stringr::str_detect(
            .data$CONSEQUENCE, "frameshift"
          ),
        as.logical(TRUE),
        as.logical(.data$is_driver)
      )) |>
      dplyr::select(c("chr", "from", "to", "ref",
                    "alt", "FILTER", "DP", "NV",
                    "VAF", "CONSEQUENCE", "GENE",
                    "is_driver", "driver_label"))
  }else{
    ## Make a single dummy mutation if mutations are absent
    data('example_dataset_CNAqc', package = 'CNAqc')
    mutations <- example_dataset_CNAqc$mutations[3,]
  }

  cna <- callset_cna$variant |>
    dplyr::mutate(
      chr = paste0('chr', .data$CHROM),
      from = .data$SEGMENT_START,
      to = .data$SEGMENT_END,
      length = .data$SEGMENT_END -
        .data$SEGMENT_START,
      covRatio = as.numeric(NA),
      Major = .data$CN_MAJOR,
      minor = .data$CN_MINOR) |>
    dplyr::select(
      c("chr",
        "from",
        "to",
        "length",
        "covRatio",
        "Major",
        "minor")) |>
    dplyr::distinct()

  purity <- NULL
  ref <- "hg38"
  if (settings$conf$sample_properties$purity != "NA") {
    purity <- settings$conf$sample_properties$purity
  }
  if (settings$genome_assembly == "grch37") {
    ref <- "hg19"
  }

  cnaqc_data <- list()
  cnaqc_data[['mutations']] <- mutations
  cnaqc_data[['cna']] <- cna
  cnaqc_data[['purity']] <- purity
  cnaqc_data[['ref']] <- ref
  cnaqc_data[['sample']] <- settings$sample_id

  cnaqc <- suppressMessages(CNAqc::init(
   mutations = cnaqc_data$mutations,
   cna = cnaqc_data$cna,
   purity = cnaqc_data$purity,
   sample = cnaqc_data$sample,
   ref = cnaqc_data$ref
  ))

  return(cnaqc)

}
