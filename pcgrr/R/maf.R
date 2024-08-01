#' Function that takes a MAF file generated with vcf2maf
#' and filters out variants that are presumably germline (tumor-only run)
#'
#' @param callset Callset with pre-processed somatic SNV/InDel variants
#' @param settings PCGR run/configuration settings
#'
#' @export
filter_maf_file <- function(callset, settings) {

  ## check if vcf2maf is TRUE
  if (as.logical(settings[['conf']][['other']][['vcf2maf']]) == FALSE) {
    return(0)
  }

  filtered_vars_maf_like <- data.frame()

  if("variant" %in% names(callset)) {
    if(NROW(callset[['variant']]) == 0) {
      return(0)
    }

    pcgrr::log4r_info(paste0(
      "Updating MAF file with filtered somatic SNV/InDels"))

    if(all(c("CHROM",
             "POS",
             "REF",
             "ALT") %in% colnames(callset[['variant']]))) {
      filtered_vars_maf_like <- callset[['variant']] |>
        dplyr::select(c("CHROM", "POS", "REF", "ALT")) |>
        dplyr::mutate(Tumor_Seq_Allele1 = dplyr::case_when(
          nchar(.data$REF) > 1 &
            nchar(.data$ALT) == 1 &
            substr(.data$REF, 1, 1) == .data$ALT ~ substr(.data$REF, 2, nchar(.data$REF)),
          nchar(.data$ALT) > 1 &
            nchar(.data$REF) == 1 &
            substr(.data$ALT, 1, 1) == .data$REF ~ "-", # insertion
          TRUE ~ .data$REF
        )) |>
        dplyr::mutate(Tumor_Seq_Allele2 = dplyr::case_when(
          nchar(.data$REF) > 1 &
            nchar(.data$ALT) == 1 &
            substr(.data$REF, 1, 1) == .data$ALT ~ "-",
          nchar(.data$ALT) > 1 &
            nchar(.data$REF) == 1 &
            substr(.data$ALT, 1, 1) == .data$REF ~ substr(.data$ALT, 2, nchar(.data$ALT)),
          TRUE ~ .data$ALT
        )) |>
        dplyr::mutate(Variant_Type = dplyr::case_when(
          nchar(.data$REF) == 1 &
            nchar(.data$ALT) == 1 ~ "SNP",
          nchar(.data$REF) > 1 &
            nchar(.data$ALT) == 1 ~ "DEL",
          nchar(.data$ALT) > 1 &
            nchar(.data$REF) == 1 ~ "INS",
          TRUE ~ "MNP"
        )) |>
        dplyr::mutate(
          Chromosome = .data$CHROM,
          Start_Position = dplyr::case_when(
            .data$Variant_Type == "DEL" &
              substr(.data$REF, 1, 1) == .data$ALT ~ .data$POS + 1,
            TRUE ~ .data$POS
          )
        ) |>
        dplyr::select(
          c("Chromosome",
            "Start_Position",
            "Tumor_Seq_Allele1",
            "Tumor_Seq_Allele2",
            "Variant_Type")
        )
    }

  }

  maf_data_unfiltered <- data.frame()
  maf_data_header <- NULL

  ## check if unfiltered MAF file exists and read it - if not, return 0
  if (file.exists(settings[['molecular_data']][['fname_maf_tsv']])) {
    if(!(file.size(settings[['molecular_data']][['fname_maf_tsv']]) == 0)) {
      maf_data_header <- readLines(
        settings[['molecular_data']][['fname_maf_tsv']], n = 1)

      maf_data_unfiltered <- readr::read_tsv(
        settings[['molecular_data']][['fname_maf_tsv']],
        show_col_types = F, col_names = T,
        comment = "#", na = ""
      )
    } else {
      pcgrr::log4r_warn("MAF file is empty - no filtering will be performed")
      return(0)
    }
  }

  if(NROW(maf_data_unfiltered) > 0 &
     NROW(filtered_vars_maf_like) > 0) {
    if(all(c("Chromosome",
             "Start_Position",
             "Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2",
             "Variant_Type") %in% colnames(maf_data_unfiltered)) &
       all(c("Chromosome",
             "Start_Position",
             "Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2",
             "Variant_Type") %in% colnames(filtered_vars_maf_like))) {
      maf_data_filtered <- maf_data_unfiltered |>
        dplyr::semi_join(
          filtered_vars_maf_like,
          by = c("Chromosome", "Start_Position",
                 "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
                 "Variant_Type")
        )

      if(!is.null(maf_data_header) &
         NROW(maf_data_filtered) > 0) {
        file.remove(settings[['molecular_data']][['fname_maf_tsv']])
        writeLines(maf_data_header,
                   settings[['molecular_data']][['fname_maf_tsv']])
        options(scipen = 999)
        readr::write_tsv(
          maf_data_filtered,
          settings[['molecular_data']][['fname_maf_tsv']],
          append = TRUE, col_names = T, quote = "none", na = "")
      }
    }
  }

  return(0)
}
