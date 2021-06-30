
max_af_gnomad <- function(sample_calls){
  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## set maximum AF from gnomAD (all populations)
  gnomad_cols <- c("GLOBAL_AF_GNOMAD", "NFE_AF_GNOMAD",
                   "AMR_AF_GNOMAD",
                   "AFR_AF_GNOMAD", "SAS_AF_GNOMAD",
                   "EAS_AF_GNOMAD", "ASJ_AF_GNOMAD",
                   "FIN_AF_GNOMAD", "OTH_AF_GNOMAD")
  sample_calls$MAX_AF_GNOMAD <- 0
  for (c in gnomad_cols) {
    if(c %in% colnames(sample_calls)){
      if (nrow(
        sample_calls[!is.na(sample_calls[, c]) &
                     sample_calls[, c] > sample_calls$MAX_AF_GNOMAD, ]) > 0) {
        sample_calls[!is.na(sample_calls[, c]) &
                       sample_calls[, c] > sample_calls$MAX_AF_GNOMAD,
                     "MAX_AF_GNOMAD"] <-
          sample_calls[!is.na(sample_calls[, c]) &
                         sample_calls[, c] > sample_calls$MAX_AF_GNOMAD, c]
      }
    }

  }
  return(sample_calls)
}

max_af_onekg <- function(sample_calls){
  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## set maximum AF from 1000 Genomes Project (all populations)
  onekg_cols <- c("GLOBAL_AF_1KG", "AMR_AF_1KG", "AFR_AF_1KG",
                  "EAS_AF_1KG", "SAS_AF_1KG", "EUR_AF_1KG")

  sample_calls$MAX_AF_1KG <- 0
  for (c in onekg_cols) {
    if (c %in% colnames(sample_calls)){
      if (nrow(
        sample_calls[!is.na(sample_calls[, c]) &
                     sample_calls[, c] > sample_calls$MAX_AF_1KG,]) > 0) {
        sample_calls[!is.na(sample_calls[, c]) &
                     sample_calls[, c] > sample_calls$MAX_AF_1KG,
                     "MAX_AF_1KG"] <-
        sample_calls[!is.na(sample_calls[, c]) &
                       sample_calls[, c] > sample_calls$MAX_AF_1KG, c]
      }
    }

  }

  return(sample_calls)
}

clinvar_germline_status <- function(sample_calls){

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## assign STATUS_CLINVAR_GERMLINE status to all calls recorded
  ## in ClinVar with a "germline" variant-of-origin
  if (("CLINVAR_MSID" %in% colnames(sample_calls)) &
      ("CLINVAR_VARIANT_ORIGIN" %in% colnames(sample_calls))) {
    sample_calls <- sample_calls %>%
      dplyr::mutate(
        STATUS_CLINVAR_GERMLINE =
          dplyr::if_else(
            !is.na(CLINVAR_MSID) &
              stringr::str_detect(CLINVAR_VARIANT_ORIGIN, "germline") &
              !stringr::str_detect(CLINVAR_VARIANT_ORIGIN, "somatic"),
            TRUE, FALSE))
  }
  return(sample_calls)
}


dbsnp_germline_status <- function(sample_calls){

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## assign STATUS_DBSNP_GERMLINE status to all calls recorded in
  ## dbSNP (except relevant in a somatic setting, as defined by ClinVar/DoCM)
  if ("DBSNPRSID" %in% colnames(sample_calls) &
      "DOCM_PMID" %in% colnames(sample_calls) &
      "CLINVAR_MSID" %in% colnames(sample_calls) &
      "CLINVAR_VARIANT_ORIGIN" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls %>%
      dplyr::mutate(
        STATUS_DBSNP_GERMLINE =
          dplyr::if_else(!is.na(DBSNPRSID), TRUE, FALSE)) %>%
      dplyr::mutate(
        STATUS_DBSNP_GERMLINE =
          dplyr::if_else(STATUS_DBSNP_GERMLINE == T &
                           !is.na(DOCM_PMID),
                         FALSE, STATUS_DBSNP_GERMLINE)) %>%
      dplyr::mutate(
        STATUS_DBSNP_GERMLINE =
          dplyr::if_else(
            STATUS_DBSNP_GERMLINE == T &
              !is.na(CLINVAR_MSID) &
              stringr::str_detect(CLINVAR_VARIANT_ORIGIN, "somatic"),
            FALSE, STATUS_DBSNP_GERMLINE))
  }
  return(sample_calls)
}

tcga_somatic_status <- function(sample_calls){

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## assign STATUS_TCGA_SOMATIC to calls with presence
  ## in any of the TCGA cohorts (TCGA_PANCANCER_COUNT > 0)
  if ("TCGA_PANCANCER_COUNT" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls %>%
      dplyr::mutate(
        STATUS_TCGA_SOMATIC =
          dplyr::if_else(!is.na(TCGA_PANCANCER_COUNT), TRUE, FALSE))
  }
  return(sample_calls)

}

cosmic_somatic_status <- function(sample_calls){

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )

  ## assign STATUS_COSMIC to all calls with an identifier in COSMIC
  if ("COSMIC_MUTATION_ID" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls %>%
      dplyr::mutate(
        STATUS_COSMIC =
          dplyr::if_else(!is.na(COSMIC_MUTATION_ID), TRUE, FALSE))
  }
  return(sample_calls)

}

hom_af_status <- function(sample_calls){


  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                          msg = paste0("Argument 'sample_calls' must be of ",
                                       "type data.frame"))
  )
  ## assign STATUS_LIKELY_GERMLINE_HOMOZYGOUS to all calls
  ## with 100% allelic fraction of alternative allele
  if ("AF_TUMOR" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls %>%
      dplyr::mutate(
        STATUS_LIKELY_GERMLINE_HOMOZYGOUS =
          dplyr::if_else(!is.na(AF_TUMOR) & AF_TUMOR == 1,
                         TRUE, FALSE))
  }
  return(sample_calls)
}

pon_status <- function(sample_calls){

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## assign STATUS_PON to all calls overlapping the
  ## user-defined panel-of-normals VCF ("PANEL_OF_NORMALS" == T)
  if ("PANEL_OF_NORMALS" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls %>%
      dplyr::mutate(
        STATUS_PON =
          dplyr::if_else(PANEL_OF_NORMALS == TRUE,
                         TRUE, FALSE))
  }
  return(sample_calls)
}

het_af_germline_status <- function(sample_calls){

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## assign STATUS_LIKELY_GERMLINE_HETEROZYGOUS to all calls
  ## that i) have the alternative allele
  ## in the [0.40,0.60] AF range, ii) are registered in dbSNP,
  ## iii) in gnomAD
  ## and iv) not present in COSMIC/TCGA
  if ("AF_TUMOR" %in% colnames(sample_calls) &
      "MAX_AF_GNOMAD" %in% colnames(sample_calls) &
      "STATUS_COSMIC" %in% colnames(sample_calls) &
      "STATUS_TCGA_SOMATIC" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls %>%
      dplyr::mutate(
        STATUS_LIKELY_GERMLINE_HETEROZYGOUS =
          dplyr::if_else(!is.na(MAX_AF_GNOMAD) &
                           STATUS_DBSNP_GERMLINE == TRUE &
                           !is.na(AF_TUMOR) &
                           AF_TUMOR >= 0.40 & AF_TUMOR <= 0.60 &
                           STATUS_TCGA_SOMATIC == FALSE &
                           STATUS_COSMIC == FALSE, TRUE, FALSE))
  }
  return(sample_calls)
}


#' Function that assigns a SOMATIC_CLASSIFICATION to variants
#' based on evidence found in variant set,
#' potentially limited by user-defined options
#'
#' @param sample_calls data frame with variants
#' @param config configuration object
#'
#' @return sample_calls
#'

assign_somatic_classification <- function(sample_calls, config) {

  sample_calls$SOMATIC_CLASSIFICATION <- "SOMATIC"

  ## Assign non-somatic classification based on various evidence criteria
  ## 1) Frequency of variant in any of the five 1000 Genomes
  ##    superpopulations is greater than the defined thresholds by the user
  ## 2) Frequency of variant in any of the gnomAD populations is
  ##    greater than the defined thresholds by the user
  ## 3) Variant is recorded in ClinVar as germline
  ## 4) Variant is found in the user-defined panel-of-normals VCF
  ## 5) Evidence for a likely homozygous germline variant -
  ##      allelic fraction in tumor sample (AF_TUMOR) is 100%
  ##      (vary rare for true somatic variants)
  ## 6) Evidence for a likely heterozygous germline variant must
  ##    satisfy three criteria:
  ##    i) Allelic fraction of alternative allele in tumor sample
  ##        (AF_TUMOR) is 40-60%,
  ##    ii) Variant is present in dbSNP AND gnomAD
  ##    iii) Variant is neither in COSMIC nor TCGA
  ## 7) Variant is recorded in dbSNP (non-somatic ClinVar/DoCM/COSMIC/TCGA)

  sample_calls <- sample_calls %>%
    dplyr::mutate(
      SOMATIC_CLASSIFICATION =
        dplyr::if_else(STATUS_POPFREQ_1KG_ABOVE_TOLERATED == TRUE,
                       "GERMLINE_1KG", SOMATIC_CLASSIFICATION)) %>%
    dplyr::mutate(
      SOMATIC_CLASSIFICATION =
        dplyr::if_else(STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED == TRUE &
                         SOMATIC_CLASSIFICATION == "SOMATIC",
                       "GERMLINE_GNOMAD", SOMATIC_CLASSIFICATION)) %>%
    dplyr::mutate(
      SOMATIC_CLASSIFICATION =
        dplyr::if_else(STATUS_CLINVAR_GERMLINE == TRUE &
                         SOMATIC_CLASSIFICATION == "SOMATIC",
                       "GERMLINE_CLINVAR", SOMATIC_CLASSIFICATION)) %>%
    dplyr::mutate(
      SOMATIC_CLASSIFICATION =
        dplyr::if_else(STATUS_PON == TRUE &
                         config[["tumor_only"]][["exclude_pon"]] == TRUE &
                         SOMATIC_CLASSIFICATION == "SOMATIC",
                       "GERMLINE_PON", SOMATIC_CLASSIFICATION)) %>%

    dplyr::mutate(
      SOMATIC_CLASSIFICATION =
        dplyr::if_else(
          STATUS_LIKELY_GERMLINE_HOMOZYGOUS == TRUE &
            config[["tumor_only"]][["exclude_likely_hom_germline"]] == TRUE &
            SOMATIC_CLASSIFICATION == "SOMATIC",
          "GERMLINE_HOMOZYGOUS", SOMATIC_CLASSIFICATION)) %>%
    dplyr::mutate(
      SOMATIC_CLASSIFICATION =
        dplyr::if_else(
          STATUS_LIKELY_GERMLINE_HETEROZYGOUS == TRUE &
            config[["tumor_only"]][["exclude_likely_het_germline"]] == TRUE &
            SOMATIC_CLASSIFICATION == "SOMATIC",
          "GERMLINE_HETEROZYGOUS", SOMATIC_CLASSIFICATION))

  ## set variants found in DBSNP as germline if this option is set to TRUE
  if (config[["tumor_only"]][["exclude_dbsnp_nonsomatic"]] == TRUE) {

    sample_calls <- sample_calls %>%
      dplyr::mutate(
        SOMATIC_CLASSIFICATION =
          dplyr::if_else(STATUS_DBSNP_GERMLINE == TRUE &
                           STATUS_TCGA_SOMATIC == FALSE &
                           STATUS_COSMIC == FALSE &
                           SOMATIC_CLASSIFICATION == "SOMATIC",
                         "GERMLINE_DBSNP", SOMATIC_CLASSIFICATION))

  }

  return(sample_calls)
}

#' Function that appends several tags denoting
#' evidence for somatic/germline status of variants
#'
#' @param sample_calls data frame with variants
#' @param config configuration object
#'
#' @return sample_calls
#'

assign_somatic_germline_evidence <- function(sample_calls, config) {

  invisible(
    assertthat::assert_that(
      is.data.frame(sample_calls),
      msg = paste0("Argument 'sample_calls' must be of ",
                   "type data.frame"))
  )
  ## assign STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED
  for (pop in c("EUR", "AMR", "AFR", "SAS", "EAS", "GLOBAL")) {
    sample_calls <-
      pcgrr::assign_germline_popfreq_status(
        sample_calls,
        pop = pop,
        dbquery = "1KG",
        max_tolerated_af =
          config[["tumor_only"]][[paste0("maf_onekg_", tolower(pop))]])
  }

  ## assign STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED
  for (pop in c("GLOBAL", "NFE", "AMR", "AFR",
                "SAS", "EAS", "ASJ", "FIN", "OTH")) {
    sample_calls <-
      pcgrr::assign_germline_popfreq_status(
        sample_calls,
        pop = pop,
        dbquery = "gnomAD",
        max_tolerated_af =
          config[["tumor_only"]][[paste0("maf_gnomad_", tolower(pop))]])
  }

  sample_calls <- sample_calls %>%
    max_af_gnomad() %>%
    max_af_onekg %>%
    dbsnp_germline_status() %>%
    clinvar_germline_status() %>%
    tcga_somatic_status() %>%
    cosmic_somatic_status() %>%
    hom_af_status() %>%
    pon_status() %>%
    het_af_germline_status()

  ## assign MAX_AF_1KG / MAX_AF_GNOMAD
  # gnomad_cols <- c("GLOBAL_AF_GNOMAD", "NFE_AF_GNOMAD", "AMR_AF_GNOMAD",
  #                  "AFR_AF_GNOMAD", "SAS_AF_GNOMAD",
  #                  "EAS_AF_GNOMAD", "ASJ_AF_GNOMAD",
  #                  "FIN_AF_GNOMAD", "OTH_AF_GNOMAD")
  # onekg_cols <- c("GLOBAL_AF_1KG", "AMR_AF_1KG", "AFR_AF_1KG",
  #                 "EAS_AF_1KG", "SAS_AF_1KG", "EUR_AF_1KG")
  #
  # sample_calls$MAX_AF_1KG <- 0
  # for (c in onekg_cols) {
  #   if (nrow(sample_calls[!is.na(sample_calls[, c]) &
  #                        sample_calls[, c] > sample_calls$MAX_AF_1KG,]) > 0) {
  #     sample_calls[!is.na(sample_calls[, c]) &
  #                    sample_calls[, c] > sample_calls$MAX_AF_1KG, "MAX_AF_1KG"] <-
  #       sample_calls[!is.na(sample_calls[, c]) &
  #                      sample_calls[, c] > sample_calls$MAX_AF_1KG, c]
  #   }
  #
  # }
  #
  # sample_calls$MAX_AF_GNOMAD <- 0
  # for (c in gnomad_cols) {
  #   if (nrow(sample_calls[!is.na(sample_calls[, c]) &
  #                        sample_calls[, c] > sample_calls$MAX_AF_GNOMAD, ]) > 0) {
  #     sample_calls[!is.na(sample_calls[, c]) &
  #                    sample_calls[, c] > sample_calls$MAX_AF_GNOMAD,
  #                  "MAX_AF_GNOMAD"] <-
  #       sample_calls[!is.na(sample_calls[, c]) &
  #                      sample_calls[, c] > sample_calls$MAX_AF_GNOMAD, c]
  #   }
  #
  # }
  #
  # ## assign STATUS_DBSNP_GERMLINE status to all calls recorded in
  # ## dbSNP (except relevant in a somatic setting, as defined by ClinVar/DoCM)
  # if ("DBSNPRSID" %in% colnames(sample_calls)) {
  #   sample_calls <- sample_calls %>%
  #     dplyr::mutate(
  #       STATUS_DBSNP_GERMLINE =
  #                     dplyr::if_else(!is.na(DBSNPRSID), TRUE, FALSE)) %>%
  #     dplyr::mutate(
  #       STATUS_DBSNP_GERMLINE =
  #         dplyr::if_else(STATUS_DBSNP_GERMLINE == T & !is.na(DOCM_PMID),
  #                        FALSE, STATUS_DBSNP_GERMLINE)) %>%
  #     dplyr::mutate(
  #       STATUS_DBSNP_GERMLINE =
  #         dplyr::if_else(
  #           STATUS_DBSNP_GERMLINE == T &
  #             !is.na(CLINVAR_MSID) &
  #             stringr::str_detect(CLINVAR_VARIANT_ORIGIN, "somatic"),
  #           FALSE, STATUS_DBSNP_GERMLINE))
  # }

  ## assign STATUS_CLINVAR_GERMLINE status to all calls recorded
  ## in ClinVar with a "germline" variant-of-origin
  # if ("CLINVAR_MSID" %in% colnames(sample_calls)) {
  #   sample_calls <- sample_calls %>%
  #     dplyr::mutate(
  #       STATUS_CLINVAR_GERMLINE =
  #         dplyr::if_else(
  #           !is.na(CLINVAR_MSID) &
  #             stringr::str_detect(CLINVAR_VARIANT_ORIGIN, "germline") &
  #             !stringr::str_detect(CLINVAR_VARIANT_ORIGIN, "somatic"),
  #           TRUE, FALSE))
  # }

  ## assign STATUS_LIKELY_GERMLINE_HOMOZYGOUS to all calls
  ## with 100% allelic fraction of alternative allele
  # if ("AF_TUMOR" %in% colnames(sample_calls)) {
  #   sample_calls <- sample_calls %>%
  #     dplyr::mutate(
  #       STATUS_LIKELY_GERMLINE_HOMOZYGOUS =
  #         dplyr::if_else(!is.na(AF_TUMOR) & AF_TUMOR == 1, TRUE, FALSE))
  # }


  # ## assign STATUS_COSMIC to all calls with an identifier in COSMIC
  # if ("COSMIC_MUTATION_ID" %in% colnames(sample_calls)) {
  #   sample_calls <- sample_calls %>%
  #     dplyr::mutate(
  #       STATUS_COSMIC =
  #         dplyr::if_else(!is.na(COSMIC_MUTATION_ID), TRUE, FALSE))
  # }

  ## assign STATUS_PON to all calls overlapping the
  ## user-defined panel-of-normals VCF
  # if ("PANEL_OF_NORMALS" %in% colnames(sample_calls)) {
  #   sample_calls <- sample_calls %>%
  #     dplyr::mutate(
  #       STATUS_PON =
  #         dplyr::if_else(PANEL_OF_NORMALS == TRUE, TRUE, FALSE))
  # }

  ## assign STATUS_LIKELY_GERMLINE_HETEROZYGOUS to all calls
  ## that i) have the alternative allele
  ## in the [0.40,0.60] AF range, ii) are registered in dbSNP,
  ## iii) in gnomAD (yet below the user-defined thresholds,
  ## and iv) not present in COSMIC/TCGA
  # if ("AF_TUMOR" %in% colnames(sample_calls) &
  #     "MAX_AF_GNOMAD" %in% colnames(sample_calls) &
  #     "STATUS_COSMIC" %in% colnames(sample_calls) &
  #     "STATUS_TCGA_SOMATIC" %in% colnames(sample_calls)) {
  #   sample_calls <- sample_calls %>%
  #     dplyr::mutate(
  #       STATUS_LIKELY_GERMLINE_HETEROZYGOUS =
  #         dplyr::if_else(!is.na(MAX_AF_GNOMAD) &
  #                          STATUS_DBSNP_GERMLINE == TRUE &
  #                          !is.na(AF_TUMOR) &
  #                          AF_TUMOR >= 0.40 & AF_TUMOR <= 0.60 &
  #                          STATUS_TCGA_SOMATIC == FALSE &
  #                          STATUS_COSMIC == FALSE, TRUE, FALSE))
  # }

  return(sample_calls)
}

#' Function that sets STATUS_POPFREQ_1KG_ABOVE_TOLERATED/
#' STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED to TRUE for variants
#' if any population frequency exceeds max_tolerated_af
#'
#' @param sample_calls data frame with variants
#' @param pop population code (1000 Genomes/gnomAD)
#' @param dbquery 1KG or gnomAD
#' @param max_tolerated_af max tolerated germline allele frequency
#'
#' @return sample_calls
#'
assign_germline_popfreq_status <- function(sample_calls,
                                           pop = "EUR",
                                           dbquery = "1KG",
                                           max_tolerated_af = 0.01) {


  if (dbquery == "1KG") {
    if (!("STATUS_POPFREQ_1KG_ABOVE_TOLERATED" %in% colnames(sample_calls))) {
      sample_calls$STATUS_POPFREQ_1KG_ABOVE_TOLERATED <- FALSE
    }
    col <- paste0(pop, "_AF_1KG")
    if (any(grepl(paste0("^", col, "$"), names(sample_calls)))) {

      sample_calls$max_tolerated_af <- max_tolerated_af
      if (nrow(
        sample_calls[!is.na(sample_calls[, col]) &
                     sample_calls[, col] > sample_calls$max_tolerated_af, ]) > 0){
        sample_calls[!is.na(sample_calls[, col]) &
                     sample_calls[, col] > sample_calls$max_tolerated_af,
                     "STATUS_POPFREQ_1KG_ABOVE_TOLERATED"] <- TRUE
      }
      sample_calls$max_tolerated_af <- NULL
    }
  }
  if (dbquery == "gnomAD") {
    if (!("STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED" %in% colnames(sample_calls))) {
      sample_calls$STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED <- FALSE
    }
    col <- paste0(pop, "_AF_GNOMAD")
    if (any(grepl(paste0("^", col, "$"), names(sample_calls)))) {

      sample_calls$max_tolerated_af <- max_tolerated_af

      if(nrow(
        sample_calls[!is.na(sample_calls[, col]) &
                     sample_calls[, col] > sample_calls$max_tolerated_af, ]) > 0){
        sample_calls[!is.na(sample_calls[, col]) &
                       sample_calls[, col] > sample_calls$max_tolerated_af,
                     "STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED"] <- TRUE
      }
      sample_calls$max_tolerated_af <- NULL
    }
  }

  return(sample_calls)
}


#' Function that retrieves name of VCF INFO tag and
#' population description for gnomad/1000G population
#'
#' @param population_code three-letter code
#' @param db 1KG or GNOMAD
#' @param subset NA or "non_cancer" (for GNOMAD)
#'
#' @return pop_tag_info
#'
get_population_tag <- function(population_code, db = "1KG", subset = NA) {
  pop_tag_info <-
    list("vcf_tag" = paste0(toupper(population_code), "_AF_", db),
         "pop_description" = NA)
  if (db == "GNOMAD" & subset == "non_cancer") {
    pop_tag_info <-
      list("vcf_tag" =
             paste0("NON_CANCER_AF_", toupper(population_code)),
           "pop_description" = NA)
  }

  pop_descriptions_1KG <-
    data.frame(code = "afr",
               pop_description = "African", stringsAsFactors = F) %>%
    rbind(data.frame(
      code = "amr",
      pop_description = "Admixed American", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "eur",
      pop_description = "European", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "eas",
      pop_description = "East Asian", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "sas",
      pop_description = "South Asian", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "global",
      pop_description = "global", stringsAsFactors = F))

  pop_descriptions_gnomad <-
    data.frame(code = "afr",
               pop_description = "African", stringsAsFactors = F) %>%
    rbind(data.frame(
      code = "amr",
      pop_description = "Admixed American", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "nfe",
      pop_description = "Non-Finnish European", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "fin",
      pop_description = "Finnish", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "oth",
      pop_description = "Other", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "asj",
      pop_description = "Ashkenazi Jewish", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "eas",
      pop_description = "East Asian", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "sas",
      pop_description = "South Asian", stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "global",
      pop_description = "global", stringsAsFactors = F))

  pop_descriptions_gnomad_non_cancer <-
    data.frame(code = "afr",
               pop_description = "African non-cancer subset",
               stringsAsFactors = F) %>%
    rbind(data.frame(
      code = "amr",
      pop_description = "Admixed American non-cancer subset",
      stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "nfe",
      pop_description = "Non-Finnish European non-cancer subset",
      stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "fin",
      pop_description = "Finnish non-cancer subset",
      stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "oth",
      pop_description = "Other non-cancer subset",
      stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "asj",
      pop_description = "Ashkenazi Jewish non-cancer subset",
      stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "eas",
      pop_description = "East Asian non-cancer subset",
      stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "sas",
      pop_description = "South Asian non-cancer subset",
      stringsAsFactors = F)) %>%
    rbind(data.frame(
      code = "global",
      pop_description = "Global non-cancer subset",
      stringsAsFactors = F))

  if (db == "1KG") {
    pop_entry <- dplyr::filter(pop_descriptions_1KG,
                               code == population_code)
    pop_tag_info[["pop_description"]] <- pop_entry$pop_description
  }
  if (db == "GNOMAD") {
    pop_entry <- dplyr::filter(pop_descriptions_gnomad,
                               code == population_code)
    pop_tag_info[["pop_description"]] <- pop_entry$pop_description
    if (subset == "non_cancer") {
      pop_entry <- dplyr::filter(pop_descriptions_gnomad_non_cancer,
                                 code == population_code)
      pop_tag_info[["pop_description"]] <- pop_entry$pop_description
    }

  }
  return(pop_tag_info)
}

#' Function that makes input data for an UpSet plot
#' (filtering/intersection results) for the somatic-germline
#' classification procedure
#'
#' @param calls unfiltered calls (germline + somatic)
#' @param config config
#'
#' @return upset data
#'
make_upset_plot_data <- function(calls, config) {

  columns <- c()
  if (config[["tumor_only"]][["exclude_pon"]] == TRUE) {
    columns <- c(columns, "STATUS_PON")
  }
  if (config[["tumor_only"]][["exclude_likely_hom_germline"]] == TRUE) {
    columns <- c(columns, "STATUS_LIKELY_GERMLINE_HOMOZYGOUS")
  }
  if (config[["tumor_only"]][["exclude_likely_het_germline"]] == TRUE) {
    columns <- c(columns, "STATUS_LIKELY_GERMLINE_HETEROZYGOUS")
  }
  if (config[["tumor_only"]][["exclude_dbsnp_nonsomatic"]] == TRUE) {
    columns <- c(columns, "STATUS_DBSNP_GERMLINE")
  }
  assertable::assert_colnames(
    calls, c("VAR_ID",
             "STATUS_POPFREQ_1KG_ABOVE_TOLERATED",
             "STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED",
             "STATUS_CLINVAR_GERMLINE"),
    only_colnames = F, quiet = T)
  df <- dplyr::select(calls, VAR_ID, STATUS_POPFREQ_1KG_ABOVE_TOLERATED,
                      STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED,
                      STATUS_CLINVAR_GERMLINE)
  for (c in columns) {
    if (c %in% colnames(calls)) {
      df[, c] <- calls[, c]
    }
  }

  for (v in colnames(df)) {
    if (v != "VAR_ID") {
      df[, v] <- as.integer(df[, v])
    }
  }
  df <- dplyr::rename(df, gnomAD = STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED,
                      OneKGP = STATUS_POPFREQ_1KG_ABOVE_TOLERATED,
                      ClinVar = STATUS_CLINVAR_GERMLINE)
  if ("STATUS_PON" %in% colnames(df)) {
    df <- dplyr::rename(df, Panel_Of_Normals = STATUS_PON)
  }
  if ("STATUS_LIKELY_GERMLINE_HOMOZYGOUS" %in% colnames(df)) {
    df <- dplyr::rename(df, HomAF = STATUS_LIKELY_GERMLINE_HOMOZYGOUS)
  }
  if ("STATUS_LIKELY_GERMLINE_HETEROZYGOUS" %in% colnames(df)) {
    df <- dplyr::rename(df, HetAF = STATUS_LIKELY_GERMLINE_HETEROZYGOUS)
  }
  if ("STATUS_DBSNP_GERMLINE" %in% colnames(df)) {
    df <- dplyr::rename(df, dbSNP = STATUS_DBSNP_GERMLINE)
  }
  return(df)

}
#' Function that makes an upset calls for germline-filtered variants
#' classification procedure
#'
#' @param upset_data unfiltered calls (germline + somatic)
#'
#' @return p
#'
upset_plot_tumor_only <- function(upset_data) {

  isets <- c()
  all_negative_filters <- c()
  for (c in colnames(upset_data)) {
    if (c != "VAR_ID") {
      if (length(unique(upset_data[, c] == 1)) == 1) {
        if (unique(upset_data[, c] == 1) == T) {
          isets <- c(isets, c)
        }else{
          all_negative_filters <- c(all_negative_filters, c)
        }
      }else{
        isets <- c(isets, c)
      }
    }
  }

  for (m in all_negative_filters) {
    upset_data[, m] <- NULL
  }

  p <- UpSetR::upset(upset_data, sets = isets,
                     sets.bar.color = "#56B4E9",
                     order.by = "freq", nintersects = 20,
                     text.scale = 1.5, show.numbers = T,
                     point.size = 6, color.pal = "Blues",
                     empty.intersections = "on")
  return(p)

}
