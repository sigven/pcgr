
#' Function that assigns a maximum value to a variable (MAX_AF_GNOMAD) reflecting
#' the maximum allele frequency for a given variant across gnomAD populations
#'
#' @param sample_calls data frame with sample variant calls
#'
#' @export
max_af_gnomad <- function(sample_calls) {
  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## set maximum AF from gnomAD (all populations)
  gnomad_cols <- c("gnomADg_AF",
                   "gnomADg_NFE_AF",
                   "gnomADg_AMR_AF",
                   "gnomADg_AFR_AF",
                   "gnomADg_SAS_AF",
                   "gnomADg_EAS_AF",
                   "gnomADg_ASJ_AF",
                   "gnomADg_FIN_AF",
                   "gnomADg_OTH_AF",
                   "gnomADe_AF",
                   "gnomADe_NFE_AF",
                   "gnomADe_AMR_AF",
                   "gnomADe_AFR_AF",
                   "gnomADe_SAS_AF",
                   "gnomADe_EAS_AF",
                   "gnomADe_ASJ_AF",
                   "gnomADe_FIN_AF",
                   "gnomADe_OTH_AF")
  sample_calls$MAX_AF_GNOMAD <- 0
  for (c in gnomad_cols) {
    if (c %in% colnames(sample_calls)) {
      if (NROW(
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


#' Function that assigns a logical to STATUS_CLINVAR_GERMLINE based on
#' whether a ClinVar entry of germline origin is found for a given variant
#' (for entries in a data frame)
#'
#' @param sample_calls data frame with sample calls
#'
#'
#' @export
clinvar_germline_status <- function(sample_calls) {

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## assign STATUS_CLINVAR_GERMLINE status to all calls recorded
  ## in ClinVar with a "germline" variant-of-origin
  if (("CLINVAR_MSID" %in% colnames(sample_calls)) &
      ("CLINVAR_VARIANT_ORIGIN" %in% colnames(sample_calls))) {
    sample_calls <- sample_calls |>
      dplyr::mutate(
        STATUS_CLINVAR_GERMLINE =
          dplyr::if_else(
            !is.na(.data$CLINVAR_MSID) &
              stringr::str_detect(.data$CLINVAR_VARIANT_ORIGIN, "germline") &
              !stringr::str_detect(.data$CLINVAR_VARIANT_ORIGIN, "somatic"),
            TRUE, FALSE))
  }
  return(sample_calls)
}

#' Function that assigns a logical (STATUS_DBSNP) reflecting whether
#' a variant co-incides with an entry in dbSNP (germline)
#'
#' @param sample_calls data frame with sample variant calls
#'
#' @export
dbsnp_germline_status <- function(sample_calls) {

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## assign STATUS_DBSNP_GERMLINE status to all calls recorded in
  ## dbSNP (except relevant in a somatic setting, as defined by ClinVar/DoCM)
  if ("DBSNP_RSID" %in% colnames(sample_calls) &
      "CLINVAR_MSID" %in% colnames(sample_calls) &
      "CLINVAR_VARIANT_ORIGIN" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls |>
      dplyr::mutate(
        STATUS_DBSNP =
          dplyr::if_else(!is.na(.data$DBSNP_RSID), TRUE, FALSE)) |>
      dplyr::mutate(
        STATUS_DBSNP =
          dplyr::if_else(
            .data$STATUS_DBSNP == T &
              !is.na(.data$CLINVAR_MSID) &
              stringr::str_detect(.data$CLINVAR_VARIANT_ORIGIN, "somatic"),
            FALSE,
            .data$STATUS_DBSNP))
  }
  return(sample_calls)
}

#' Function that assigns a logical (STATUS_TCGA_SOMATIC) reflecting whether
#' a variant co-incides with an entry in TCGA (somatic)
#'
#' @param sample_calls data frame with sample variant calls
#'
#' @export
tcga_somatic_status <- function(sample_calls) {

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## assign STATUS_TCGA_SOMATIC to calls with presence
  ## in any of the TCGA cohorts (TCGA_PANCANCER_COUNT > 0)
  if ("TCGA_PANCANCER_COUNT" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls |>
      dplyr::mutate(
        STATUS_TCGA =
          dplyr::if_else(
            !is.na(.data$TCGA_PANCANCER_COUNT),
            TRUE, FALSE))
  }
  return(sample_calls)

}

#' Function that assigns a logical (STATUS_COSMIC) reflecting whether
#' a variant co-incides with an entry in COSMIC (germline)
#'
#' @param sample_calls data frame with sample variant calls
#'
#' @export
cosmic_somatic_status <- function(sample_calls) {

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )

  ## assign STATUS_COSMIC to all calls with an identifier in COSMIC
  if ("COSMIC_ID" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls |>
      dplyr::mutate(
        STATUS_COSMIC =
          dplyr::if_else(
            !is.na(.data$COSMIC_ID),
            TRUE, FALSE))
  }
  return(sample_calls)

}

#' Function that assigns a logical (STATUS_LIKELY_GERMLINE_HOMOZYGOUS) reflecting whether
#' a variant is likely homozygous (germline) - based on allelic fraction (VAF_TUMOR)
#'
#' @param sample_calls data frame with sample variant calls
#'
#' @export
hom_af_status <- function(sample_calls) {


  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                          msg = paste0("Argument 'sample_calls' must be of ",
                                       "type data.frame"))
  )
  ## assign STATUS_GERMLINE_HOM to all calls
  ## with 100% allelic fraction of alternative allele
  if ("VAF_TUMOR" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls |>
      dplyr::mutate(
        STATUS_GERMLINE_HOM =
          dplyr::if_else(
            !is.na(.data$VAF_TUMOR) &
              .data$VAF_TUMOR == 1,
            TRUE,
            FALSE))
  }
  return(sample_calls)
}

#' Function that assigns a logical (STATUS_PON) reflecting whether
#' a variant is co-inciding with a variant present in a panel-of-normals database
#' (PANEL_OF_NORMALS column is TRUE)
#'
#' @param sample_calls data frame with sample variant calls
#'
#' @export
pon_status <- function(sample_calls) {

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## assign STATUS_PON to all calls overlapping the
  ## user-defined panel-of-normals VCF ("PANEL_OF_NORMALS" == T)
  if ("PANEL_OF_NORMALS" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls |>
      dplyr::mutate(
        STATUS_PON =
          dplyr::if_else(
            .data$PANEL_OF_NORMALS == TRUE,
            TRUE,
            FALSE))
  }
  return(sample_calls)
}

#' Function that assigns a logical (STATUS_LIKELY_GERMLINE_HETEROZYGOUS) reflecting whether
#' a variant is likely heterozygous (germline) - based on allelic fraction (VAF_TUMOR),
#' presence in gnomAD and dbSNP, and no presence in TCGA and COSMIC
#'
#' @param sample_calls data frame with sample variant calls
#'
#' @export
het_af_germline_status <- function(sample_calls) {

  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )
  ## assign STATUS_GERMLINE_HET to all calls
  ## that i) have the alternative allele
  ## in the [0.40,0.60] AF range, ii) are registered in dbSNP,
  ## iii) in gnomAD
  ## and iv) not present in COSMIC/TCGA
  if ("VAF_TUMOR" %in% colnames(sample_calls) &
      "MAX_AF_GNOMAD" %in% colnames(sample_calls) &
      "STATUS_COSMIC" %in% colnames(sample_calls) &
      "STATUS_TCGA" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls |>
      dplyr::mutate(
        STATUS_GERMLINE_HET =
          dplyr::if_else(
            !is.na(.data$MAX_AF_GNOMAD) &
              .data$STATUS_DBSNP == TRUE &
              !is.na(.data$VAF_TUMOR) &
              .data$VAF_TUMOR >= 0.40 & .data$VAF_TUMOR <= 0.60 &
              .data$STATUS_TCGA == FALSE &
              .data$STATUS_COSMIC == FALSE, TRUE, FALSE))
  }
  return(sample_calls)
}


#' Function that assigns a SOMATIC_CLASSIFICATION to variants
#' based on evidence found in variant set,
#' potentially limited by user-defined options
#'
#' @param sample_calls data frame with putative somatic variants
#' @param settings PCGR run/configuration settings
#'
#' @return sample_calls
#'
#' @export

assign_somatic_classification <- function(sample_calls, settings) {

  tumor_only_settings <-
    settings$conf$somatic_snv$tumor_only

  ## Assign non-somatic classification based on various evidence criteria
  ## 1) Frequency of minor allele in any of the gnomAD populations is
  ##    greater than the defined thresholds by the user (by default)
  ##
  ## User-defined filtering options:
  ##
  ## 2) Variant is recorded in ClinVar as germline
  ## 3) Variant is found in the user-defined panel-of-normals VCF
  ## 4) Evidence for a likely homozygous germline variant -
  ##      allelic fraction in tumor sample (VAF_TUMOR) is 100%
  ##      (very rare scenario for true somatic variants)
  ## 5) Evidence for a likely heterozygous germline variant must
  ##    satisfy three criteria:
  ##    i) Allelic fraction of alternative allele in tumor sample
  ##        (VAF_TUMOR) is 40-60%,
  ##    ii) Variant is present in dbSNP AND gnomAD
  ##    iii) Variant is neither in COSMIC nor TCGA
  ## 6) Variant is recorded in dbSNP (non-somatic ClinVar/COSMIC/TCGA)

  pcgrr::log4r_info("Applying variant filters on tumor-only calls - assigning somatic classification")
  sample_calls <- sample_calls |>
    dplyr::mutate(
      GERMLINE_GNOMAD =
        dplyr::if_else(
          .data$gnomADe_AF_ABOVE_TOLERATED == TRUE |
            .data$gnomADg_AF_ABOVE_TOLERATED == TRUE,
          "GERMLINE_GNOMAD",
          "")) |>
    dplyr::mutate(
      GERMLINE_CLINVAR =
        dplyr::if_else(
          .data$STATUS_CLINVAR_GERMLINE == TRUE &
            tumor_only_settings[["exclude_clinvar_germline"]] == TRUE,
          "GERMLINE_CLINVAR",
          "")) |>
    dplyr::mutate(
      GERMLINE_PON =
        dplyr::if_else(
          .data$STATUS_PON == TRUE &
            tumor_only_settings[["exclude_pon"]] == TRUE,
          "GERMLINE_PON",
          "")) |>

    dplyr::mutate(
      GERMLINE_HOM =
        dplyr::if_else(
          .data$STATUS_GERMLINE_HOM == TRUE &
            as.logical(
              tumor_only_settings[["exclude_likely_hom_germline"]]) == TRUE,
          "GERMLINE_HOM",
          "")) |>
    dplyr::mutate(
      GERMLINE_HET =
        dplyr::if_else(
          .data$STATUS_GERMLINE_HET == TRUE &
            as.logical(
              tumor_only_settings[["exclude_likely_het_germline"]]) == TRUE,
          "GERMLINE_HET",
          "")) |>
    dplyr::mutate(
      GERMLINE_DBSNP =
        dplyr::if_else(
          .data$STATUS_DBSNP == TRUE &
            .data$STATUS_TCGA == FALSE &
            .data$STATUS_COSMIC == FALSE &
            as.logical(
              tumor_only_settings[["exclude_dbsnp_nonsomatic"]]) == TRUE,
          "GERMLINE_DBSNP",
          "")) |>
    tidyr::unite("SOMATIC_CLASSIFICATION",
                 c("GERMLINE_CLINVAR",
                   "GERMLINE_DBSNP",
                   "GERMLINE_GNOMAD",
                   "GERMLINE_HET",
                   "GERMLINE_HOM",
                   "GERMLINE_PON"),
                 sep="|",
                 remove = TRUE) |>
    dplyr::mutate(
      SOMATIC_CLASSIFICATION = stringr::str_replace_all(
      .data$SOMATIC_CLASSIFICATION,
      "(\\|{1,}$)|^(\\|{1,})",""
    )) |>
    dplyr::mutate(
      SOMATIC_CLASSIFICATION = stringr::str_replace_all(
      .data$SOMATIC_CLASSIFICATION,
      "(\\|{2,})","|"
    )) |>
    dplyr::mutate(
      SOMATIC_CLASSIFICATION = dplyr::if_else(
        .data$SOMATIC_CLASSIFICATION == "",
        "SOMATIC",
        .data$SOMATIC_CLASSIFICATION
      )
    ) |>
    dplyr::select(
      -c(
        "STATUS_TCGA",
        "STATUS_COSMIC",
        "STATUS_DBSNP",
        "STATUS_GERMLINE_HET",
        "STATUS_GERMLINE_HOM",
        "STATUS_PON",
        "gnomADe_AF_ABOVE_TOLERATED",
        "gnomADg_AF_ABOVE_TOLERATED",
        "STATUS_CLINVAR_GERMLINE"))

  return(sample_calls)
}

#' Function that appends several tags denoting
#' evidence for somatic/germline status of variants
#'
#' @param sample_calls data frame with variants
#' @param settings PCGR run/configuration settings
#'
#' @return sample_calls
#'
#' @export

assign_somatic_germline_evidence <- function(
    sample_calls,
    settings = NULL) {

  invisible(
    assertthat::assert_that(
      is.data.frame(sample_calls),
      msg = paste0("Argument 'sample_calls' must be of ",
                   "type data.frame"))
  )

  tumor_only_settings <-
    settings$conf$somatic_snv$tumor_only

  ## assign STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED
  for (pop in c("GLOBAL", "NFE", "AMR", "AFR",
                "SAS", "EAS", "ASJ", "FIN", "OTH")) {
    sample_calls <-
      pcgrr::assign_germline_popfreq_status(
        sample_calls,
        pop = pop,
        dbquery = c("gnomADe","gnomADg"),
        max_tolerated_af =
          tumor_only_settings[[paste0("maf_gnomad_", tolower(pop))]])
  }

  sample_calls <- sample_calls |>
    max_af_gnomad() |>
    dbsnp_germline_status() |>
    clinvar_germline_status() |>
    tcga_somatic_status() |>
    cosmic_somatic_status() |>
    hom_af_status() |>
    pon_status() |>
    het_af_germline_status()

  return(sample_calls)
}

#' Function that sets STATUS_POPFREQ_1KG_ABOVE_TOLERATED/
#' STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED to TRUE for variants
#' if any population frequency exceeds max_tolerated_af
#'
#' @param sample_calls data frame with variants
#' @param pop population code (1000 Genomes/gnomAD)
#' @param dbquery character vector with germline db sources,
#' e.g. c("gnomADe","gnomADg")
#' @param max_tolerated_af max tolerated germline allele frequency
#'
#' @return sample_calls
#'
#' @export
assign_germline_popfreq_status <- function(sample_calls,
                                           pop = "NFE",
                                           dbquery = c("gnomADe","gnomADg"),
                                           max_tolerated_af = 0.01) {

  if(pop == "GLOBAL"){
    pop = ""
  }

  for(db in dbquery){
    if (db == "gnomADe") {
      if (!("gnomADe_AF_ABOVE_TOLERATED" %in% colnames(sample_calls))) {
        sample_calls$gnomADe_AF_ABOVE_TOLERATED <- FALSE
      }
      col <- stringr::str_replace(
        paste0(db,"_",pop, "_AF"),"__","_")
      if (any(grepl(paste0("^", col, "$"), names(sample_calls)))) {

        sample_calls$max_tolerated_af <- max_tolerated_af

        if (NROW(
          sample_calls[!is.na(sample_calls[, col]) &
                       sample_calls[, col] > sample_calls$max_tolerated_af, ]) > 0) {
          sample_calls[!is.na(sample_calls[, col]) &
                         sample_calls[, col] > sample_calls$max_tolerated_af,
                       "gnomADe_AF_ABOVE_TOLERATED"] <- TRUE
        }
        sample_calls$max_tolerated_af <- NULL
      }
    }

    if (db == "gnomADg") {
      if (!("gnomADg_AF_ABOVE_TOLERATED" %in% colnames(sample_calls))) {
        sample_calls$gnomADg_AF_ABOVE_TOLERATED <- FALSE
      }
      col <- stringr::str_replace(
        paste0(db,"_",pop, "_AF"),"__","_")
      if (any(grepl(paste0("^", col, "$"), names(sample_calls)))) {

        sample_calls$max_tolerated_af <- max_tolerated_af

        if (NROW(
          sample_calls[!is.na(sample_calls[, col]) &
                       sample_calls[, col] > sample_calls$max_tolerated_af, ]) > 0) {
          sample_calls[!is.na(sample_calls[, col]) &
                         sample_calls[, col] > sample_calls$max_tolerated_af,
                       "gnomADg_AF_ABOVE_TOLERATED"] <- TRUE
        }
        sample_calls$max_tolerated_af <- NULL
      }
    }

  }


  return(sample_calls)
}


#' Function that makes input data for an UpSet plot
#' (filtering/intersection results) for the somatic-germline
#' classification procedure
#'
#' param calls unfiltered calls (germline + somatic)
#' param config config
#'
#' return upset data
#'
#' export
#' make_upset_plot_data <- function(calls, config) {
#'
#'   columns <- c()
#'   if (config[["tumor_only"]][["exclude_pon"]] == TRUE) {
#'     columns <- c(columns, "STATUS_PON")
#'   }
#'   if (config[["tumor_only"]][["exclude_likely_hom_germline"]] == TRUE) {
#'     columns <- c(columns, "STATUS_GERMLINE_HOM")
#'   }
#'   if (config[["tumor_only"]][["exclude_likely_het_germline"]] == TRUE) {
#'     columns <- c(columns, "STATUS_GERMLINE_HET")
#'   }
#'   if (config[["tumor_only"]][["exclude_dbsnp_nonsomatic"]] == TRUE) {
#'     columns <- c(columns, "STATUS_DBSNP")
#'   }
#'   assertable::assert_colnames(
#'     calls, c("VAR_ID",
#'              "gnomADe_AF_ABOVE_TOLERATED",
#'              "STATUS_CLINVAR_GERMLINE"),
#'     only_colnames = F, quiet = T)
#'   df <- dplyr::select(calls, .data$VAR_ID,
#'                       .data$gnomADe_AF_ABOVE_TOLERATED,
#'                       .data$STATUS_CLINVAR_GERMLINE)
#'   for (c in columns) {
#'     if (c %in% colnames(calls)) {
#'       df[, c] <- calls[, c]
#'     }
#'   }
#'
#'   for (v in colnames(df)) {
#'     if (v != "VAR_ID") {
#'       df[, v] <- as.integer(df[, v])
#'     }
#'   }
#'   df <- df |>
#'     dplyr::rename(
#'       gnomAD = .data$gnomADe_AF_ABOVE_TOLERATED,
#'       ClinVar = .data$STATUS_CLINVAR_GERMLINE)
#'   if ("STATUS_PON" %in% colnames(df)) {
#'     df <- dplyr::rename(df, Panel_Of_Normals = .data$STATUS_PON)
#'   }
#'   if ("STATUS_GERMLINE_HOM" %in% colnames(df)) {
#'     df <- dplyr::rename(df, HomAF = .data$STATUS_GERMLINE_HOM)
#'   }
#'   if ("STATUS_GERMLINE_HET" %in% colnames(df)) {
#'     df <- dplyr::rename(df, HetAF = .data$STATUS_GERMLINE_HET)
#'   }
#'   if ("STATUS_DBSNP" %in% colnames(df)) {
#'     df <- dplyr::rename(df, dbSNP = .data$STATUS_DBSNP)
#'   }
#'   return(df)
#'
#' }
#'
#' #' #' Function that makes an upset calls for germline-filtered variants
#' #' classification procedure
#' #'
#' #' param upset_data unfiltered calls (germline + somatic)
#' #'
#' #' return p
#' #'
#' #' export
#' upset_plot_tumor_only <- function(upset_data) {
#'
#'   isets <- c()
#'   all_negative_filters <- c()
#'   for (c in colnames(upset_data)) {
#'     if (c != "VAR_ID") {
#'       if (length(unique(upset_data[, c] == 1)) == 1) {
#'         if (unique(upset_data[, c] == 1) == T) {
#'           isets <- c(isets, c)
#'         }else{
#'           all_negative_filters <- c(all_negative_filters, c)
#'         }
#'       }else{
#'         isets <- c(isets, c)
#'       }
#'     }
#'   }
#'
#'   for (m in all_negative_filters) {
#'     upset_data[, m] <- NULL
#'   }
#'
#'   p <- UpSetR::upset(upset_data, sets = isets,
#'                      sets.bar.color = "#56B4E9",
#'                      order.by = "freq", nintersects = 20,
#'                      text.scale = 1.5, show.numbers = T,
#'                      point.size = 6, color.pal = "Blues",
#'                      empty.intersections = "on")
#'   return(p)
#'
#' }

#' Function that generates variant filtering statistics (i.e. removal of
#' likely non-somatic/germline events) for callsets coming from
#' tumor-only sequencing
#'
#' #' param callset list object with unfiltered calls (callset$variant_unfiltered)
#' #' param settings list object with PCGR run configuration settings
#' #'
#' #' export
#' tumor_only_vfilter_stats <-
#'   function(callset = NULL,
#'            settings = NULL) {
#'
#'     ## raw, unfiltered set of (PASSED) variant calls that
#'     ## should be subject to filtering
#'     vcalls <- callset$variant_unfiltered
#'     to_settings <- settings$conf$somatic_snv$tumor_only
#'
#'     ## initiate report object with tumor-only variant filter stats
#'     to_stats <-
#'       pcgrr::init_tumor_only_content()
#'
#'     to_stats[['vfilter']][['unfiltered_n']] <-
#'       NROW(vcalls)
#'
#'     ## Assign statistics to successive filtering levels for
#'     ## different evidence criteria
#'     ## excluded germline calls found in gnomAD
#'     to_stats[["vfilter"]][["gnomad_n_remain"]] <-
#'       NROW(vcalls) -
#'       NROW(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_GNOMAD", ])
#'     # pcgrr::log4r_info(paste0("Excluding coinciding germline variants in ",
#'     #                          "gnomAD populations"))
#'     # pcgrr::log4r_info(paste0("Total sample calls remaining: ",
#'     #                          to_stats$vfilter[["gnomad_n_remain"]]))
#'
#'     ## excluded germline calls found in ClinVar
#'     to_stats$vfilter[["clinvar_n_remain"]] <-
#'       to_stats$vfilter[["gnomad_n_remain"]] -
#'       NROW(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_CLINVAR", ])
#'     # pcgrr::log4r_info(paste0("Excluding coinciding germline variants in ClinVar"))
#'     # pcgrr::log4r_info(paste0("Total sample calls remaining: ",
#'     #                          to_stats$vfilter[["clinvar_n_remain"]]))
#'
#'
#'     ## excluded germline calls found in panel of normals (if provided)
#'     to_stats$vfilter[["pon_n_remain"]] <-
#'       to_stats$vfilter[["clinvar_n_remain"]]
#'     if (as.logical(to_settings[["exclude_pon"]]) == TRUE) {
#'       to_stats$vfilter[["pon_n_remain"]] <-
#'         to_stats$vfilter[["pon_n_remain"]] -
#'         NROW(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_PON", ])
#'       # pcgrr::log4r_info(
#'       #   paste0("Excluding putative germline variants found in calls ",
#'       #          "from panel-of-normals (PON)"))
#'       # pcgrr::log4r_info(
#'       #   paste0("Total sample calls remaining: ",
#'       #          to_stats$vfilter[["pon_n_remain"]]))
#'     }
#'
#'     ## excluded germline calls found with 100% allelic fraction
#'     ## (likely homozygous germline variants)
#'     to_stats$vfilter[["hom_n_remain"]] <-
#'       to_stats$vfilter[["pon_n_remain"]]
#'     if (as.logical(to_settings[["exclude_likely_hom_germline"]]) == TRUE) {
#'       to_stats$vfilter[["hom_n_remain"]] <-
#'         to_stats$vfilter[["hom_n_remain"]] -
#'         NROW(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_HOMOZYGOUS", ])
#'       # pcgrr::log4r_info(
#'       #   paste0("Excluding likely homozygous germline variants found ",
#'       #          "as variants with 100% allelic fraction"))
#'       # pcgrr::log4r_info(paste0("Total sample calls remaining: ",
#'       #                          to_stats$vfilter[["hom_n_remain"]]))
#'     }
#'
#'     ## excluded germline calls found as likely heterozygous germline variants
#'     to_stats$vfilter[["het_n_remain"]] <-
#'       to_stats$vfilter[["hom_n_remain"]]
#'     if (as.logical(to_settings[["exclude_likely_het_germline"]]) == TRUE) {
#'       to_stats$vfilter[["het_n_remain"]] <-
#'         to_stats$vfilter[["het_n_remain"]] -
#'         NROW(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_HETEROZYGOUS", ])
#'       # pcgrr::log4r_info(paste0(
#'       #   "Excluding likely heterozygous germline variants found as variants ",
#'       #   "with 40-60% allelic fraction and recorded in gnomAD + dbSNP"))
#'       # pcgrr::log4r_info(paste0("Total sample calls remaining: ",
#'       #                          to_stats$vfilter[["het_n_remain"]]))
#'     }
#'
#'     ## excluded calls with dbSNP germline status (if set in config)
#'     to_stats$vfilter[["dbsnp_n_remain"]] <-
#'       to_stats$vfilter[["het_n_remain"]]
#'     if (as.logical(to_settings[["exclude_dbsnp_nonsomatic"]]) == TRUE) {
#'
#'       #pcgrr::log4r_info(
#'       #  paste0("Excluding non-somatically associated dbSNP variants ",
#'       #         "(dbSNP - not recorded as somatic in ClinVar",
#'       #         "and not registered in COSMIC or found in TCGA"))
#'
#'       to_stats$vfilter[["dbsnp_n_remain"]] <-
#'         to_stats$vfilter[["dbsnp_n_remain"]] -
#'         NROW(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_DBSNP", ])
#'       #pcgrr::log4r_info(paste0("Total sample calls remaining: ",
#'       #                         to_stats$vfilter[["dbsnp_n_remain"]]))
#'     }
#'
#'     ## excluded non-exonic calls (if set in config)
#'     to_stats$vfilter[["nonexonic_n_remain"]] <-
#'       to_stats$vfilter[["dbsnp_n_remain"]]
#'     if (as.logical(to_settings[["exclude_nonexonic"]]) == TRUE) {
#'
#'       #pcgrr::log4r_info(
#'       #  paste0("Excluding non-somatically associated dbSNP variants ",
#'       #         "(dbSNP - not recorded as somatic in ClinVar",
#'       #         "and not registered in COSMIC or found in TCGA"))
#'
#'       to_stats$vfilter[["nonexonic_n_remain"]] <-
#'         to_stats$vfilter[["nonexonic_n_remain"]] -
#'         NROW(vcalls[vcalls$EXONIC_STATUS == "nonexonic", ])
#'       #pcgrr::log4r_info(paste0("Total sample calls remaining: ",
#'       #                         to_stats$vfilter[["dbsnp_n_remain"]]))
#'     }
#'
#'     to_stats[["eval"]] <- TRUE
#'
#'     for (db_filter in c("gnomad", "dbsnp", "pon",
#'                         "clinvar", "hom", "het", "nonexonic")) {
#'       if (to_stats[["vfilter"]][[paste0(db_filter, "_n_remain")]] > 0 &
#'           to_stats[["vfilter"]][["unfiltered_n"]] > 0) {
#'         to_stats[["vfilter"]][[paste0(db_filter, "_frac_remain")]] <-
#'           round((as.numeric(to_stats[["vfilter"]][[paste0(db_filter,
#'                                                               "_n_remain")]]) /
#'                    to_stats[["vfilter"]][["unfiltered_n"]]) * 100, digits = 2)
#'       }
#'     }
#'     return(to_stats)
#'
#'   }


#' Function that generates a pie chart for germline filtering statistics
#' for callsets coming from tumor-only sequencing
#'
#' @param report list object with PCGR report content
#' @param plot_margin_top top margin of the plot
#' @param plot_margin_bottom bottom margin of the plot
#' @param plot_margin_left left margin of the plot
#' @param plot_margin_right right margin of the plot
#' @param font_family font family for plot text
#' @param font_size font size for plot text
#' @param pie_line_width line width for pie chart segments
#' @param opacity_filtered_categories opacity for filtered categories
#' @param hole_size_pie size of the hole in the pie chart
#'
#' @return filtering_stats list with data frame and plotly pie chart
#' @export
#'
#'
plot_filtering_stats_germline <- function(
    report = NULL,
    plot_margin_top = 50,
    plot_margin_bottom = 20,
    plot_margin_left = 20,
    plot_margin_right = 20,
    font_family = "Helvetica",
    font_size = 15,
    pie_line_width = 3,
    opacity_filtered_categories = 0.4,
    hole_size_pie = 0.4) {

  invisible(assertthat::assert_that(
    !is.null(report),
    msg = paste0("Argument 'report' must be provided")))
  invisible(assertthat::assert_that(
    is.list(report),
    msg = paste0("Argument 'report' must be of type list")))
  invisible(assertthat::assert_that(
    "content" %in% names(report),
    msg = paste0(
      "Argument 'report' must contain ",
      "a 'content' entry")))
  invisible(assertthat::assert_that(
    "snv_indel" %in% names(report$content),
    msg = paste0(
      "Argument 'report$content' must contain ",
      "a 'snv_indel' entry")))
  invisible(assertthat::assert_that(
    "callset" %in% names(report$content$snv_indel),
    msg = paste0(
      "Argument 'report$content$snv_indel' must contain ",
      "a 'callset' entry")))
  invisible(assertthat::assert_that(
    "variant_unfiltered" %in% names(report$content$snv_indel$callset),
    msg = paste0(
      "Argument 'report$content$snv_indel$callset' must contain ",
      "a 'variant_unfiltered' entry")))
  invisible(assertthat::assert_that(
    is.data.frame(report$content$snv_indel$callset$variant_unfiltered),
    msg = paste0(
      "Argument 'report$content$snv_indel$callset$variant_unfiltered' ",
      "must be of type data.frame")))

  invisible(assertable::assert_colnames(
    report$content$snv_indel$callset$variant_unfiltered,
    c("SOMATIC_CLASSIFICATION"), only_colnames = F, quiet = T))

  df <- report$content$snv_indel$callset$variant_unfiltered

  filtering_stats <- list()
  filtering_stats[['df']] <- data.frame()
  filtering_stats[['plot']] <- NULL

  if(NROW(df) > 0){

    filtering_stats[['df']] <-
      plyr::count(
        df$SOMATIC_CLASSIFICATION
      ) |>
      dplyr::arrange(dplyr::desc(freq)) |>
      dplyr::rename(FILTER = x, FREQUENCY = freq) |>
      dplyr::mutate(FILTER = stringr::str_replace_all(
        FILTER, "\\|", ", ")) |>
      dplyr::mutate(FILTER = dplyr::if_else(
        !(FILTER %in% pcgrr::germline_filter_levels) & FILTER != "SOMATIC",
        "MULTIPLE FILTERS",
        FILTER
      )) |>
      dplyr::group_by(FILTER) |>
      dplyr::reframe(FREQUENCY = sum(FREQUENCY)) |>
      dplyr::arrange(dplyr::desc(FREQUENCY)) |>
      dplyr::mutate(FILTER = factor(
        FILTER, levels = pcgrr::germline_filter_levels)) |>
      dplyr::mutate(PERCENT = scales::percent(
        FREQUENCY / sum(FREQUENCY), accuracy = 0.1))

    germline_filters <- unique(filtering_stats[['df']]$FILTER)
    germline_filters <- c("SOMATIC", setdiff(germline_filters, "SOMATIC"))
    hex_colors_germline <- head(
      pcgrr::color_palette$tier$values,
      length(germline_filters))

    if(length(germline_filters) == 2){
      plot_margin_bottom <- 100
    }

    rgba_colors <- c()
    i <- 1
    while(i <= length(germline_filters)){
      alpha <- opacity_filtered_categories
      ## full opacity for 'SOMATIC' category
      if(i == 1){
        alpha <- 1
      }
      rgba_colors <- c(
        rgba_colors,
        pcgrr::hex_to_rgba(
          hex_colors_germline[i],
          alpha = alpha))
      i <- i + 1
    }

    t <- list(
      family = font_family,
      size = font_size)

    ## pie chart for germline filtering stats
    filtering_stats[['plot']] <-
      plotly::plot_ly(
        filtering_stats[['df']],
        marker = list(
          colors = rgba_colors,
          line = list(
            color = '#FFFFFF',
            width = pie_line_width))
      ) |>
      plotly::add_pie(
        filtering_stats[['df']],
        labels =~ FILTER,
        values = ~FREQUENCY,
        textinfo = "PERCENT",
        type = 'pie',
        hole = hole_size_pie) |>
      plotly::layout(
        ## legend at the bottom
        legend = list(
          orientation = "h",
          font = t,
          xanchor = "center",
          x = 0.5),
        margin = list(
          l = plot_margin_left,  # left margin
          r = plot_margin_right,  # right margin
          b = plot_margin_bottom,  # bottom margin
          t = plot_margin_top   # top margin
        ))


  }

  return(filtering_stats)

}


#' Function that generates a pie chart for exonic/non-exonic variant
#' statistics (for callsets coming from tumor-only sequencing)
#'
#' @param report list object with PCGR report content
#' @param plot_margin_top top margin of the plot
#' @param plot_margin_bottom bottom margin of the plot
#' @param plot_margin_left left margin of the plot
#' @param plot_margin_right right margin of the plot
#' @param font_family font family for plot text
#' @param font_size font size for plot text
#' @param pie_line_width line width for pie chart segments
#' @param opacity_filtered_categories opacity for filtered categories
#' @param hole_size_pie size of the hole in the pie chart
#'
#' @return filtering_stats list with data frame and plotly pie chart
#' @export
#'
#'
plot_filtering_stats_exonic <- function(
    report = NULL,
    plot_margin_top = 50,
    plot_margin_bottom = 20,
    plot_margin_left = 20,
    plot_margin_right = 20,
    font_family = "Helvetica",
    font_size = 15,
    pie_line_width = 3,
    opacity_filtered_categories = 0.4,
    hole_size_pie = 0.4) {

  invisible(assertthat::assert_that(
    !is.null(report),
    msg = paste0("Argument 'report' must be provided")))
  invisible(assertthat::assert_that(
    is.list(report),
    msg = paste0("Argument 'report' must be of type list")))
  invisible(assertthat::assert_that(
    "content" %in% names(report),
    msg = paste0(
      "Argument 'report' must contain ",
      "a 'content' entry")))
  invisible(assertthat::assert_that(
    "snv_indel" %in% names(report$content),
    msg = paste0(
      "Argument 'report$content' must contain ",
      "a 'snv_indel' entry")))
  invisible(assertthat::assert_that(
    "callset" %in% names(report$content$snv_indel),
    msg = paste0(
      "Argument 'report$content$snv_indel' must contain ",
      "a 'callset' entry")))
  invisible(assertthat::assert_that(
    "variant_unfiltered" %in% names(report$content$snv_indel$callset),
    msg = paste0(
      "Argument 'report$content$snv_indel$callset' must contain ",
      "a 'variant_unfiltered' entry")))
  invisible(assertthat::assert_that(
    is.data.frame(report$content$snv_indel$callset$variant_unfiltered),
    msg = paste0(
      "Argument 'report$content$snv_indel$callset$variant_unfiltered' ",
      "must be of type data.frame")))

  assertable::assert_colnames(
    report$content$snv_indel$callset$variant_unfiltered,
    c("EXONIC_STATUS","SOMATIC_CLASSIFICATION"),
    only_colnames = F, quiet = T)

  df <- report$content$snv_indel$callset$variant_unfiltered |>
    dplyr::filter(.data$SOMATIC_CLASSIFICATION == "SOMATIC")

  filtering_stats <- list()
  filtering_stats[['df']] <- data.frame()
  filtering_stats[['plot']] <- NULL

  if(NROW(df) > 0){

    filtering_stats[['df']] <-
      ## Exonic status can be either of 'exonic' or 'nonexonic'
      plyr::count(
        df$EXONIC_STATUS) |>
      dplyr::arrange(dplyr::desc(freq)) |>
      dplyr::rename(FILTER = x, FREQUENCY = freq) |>
      dplyr::mutate(FILTER = toupper(paste0("SOMATIC - ", FILTER))) |>
      dplyr::group_by(FILTER) |>
      dplyr::reframe(FREQUENCY = sum(FREQUENCY)) |>
      dplyr::arrange(dplyr::desc(FREQUENCY)) |>
      dplyr::mutate(PERCENT = scales::percent(
        FREQUENCY / sum(FREQUENCY), accuracy = 0.1))

    ## if either category is zero, add a row with zero frequency
    if(!("SOMATIC - EXONIC" %in% filtering_stats[['df']]$FILTER)){
      filtering_stats[['df']] <- dplyr::bind_rows(
        filtering_stats[['df']],
        data.frame(
          FILTER = "SOMATIC - EXONIC",
          FREQUENCY = 0,
          PERCENT = "0.0%"))
    }
    if(!("SOMATIC - NONEXONIC" %in% filtering_stats[['df']]$FILTER)){
      filtering_stats[['df']] <- dplyr::bind_rows(
        filtering_stats[['df']],
        data.frame(
          FILTER = "SOMATIC - NONEXONIC",
          FREQUENCY = 0,
          PERCENT = "0.0%"))
    }

    filtering_stats$df$FILTER <- factor(
      filtering_stats$df$FILTER,
      levels = pcgrr::exonic_filter_levels)

    rgba_colors <- c(
      pcgrr::hex_to_rgba(pcgrr::color_palette$tier$values[1], alpha = 1),
      pcgrr::hex_to_rgba(
        pcgrr::color_palette$tier$values[2],
        alpha = opacity_filtered_categories
      ))

    t <- list(
      family = font_family,
      size = font_size)

    ## pie chart for germline filtering stats
    filtering_stats[['plot']] <-
      plotly::plot_ly(
        filtering_stats[['df']],
        marker = list(
          colors = rgba_colors,
          line = list(
            color = '#FFFFFF',
            width = pie_line_width))) |>
      plotly::add_pie(
        filtering_stats[['df']],
        labels =~ FILTER,
        values = ~FREQUENCY,
        textinfo = "PERCENT",
        type = 'pie',
        hole = hole_size_pie) |>
      plotly::layout(
        ## legend at the bottom
        legend = list(
          orientation = "h",
          font = t,
          xanchor = "center",
          y = -0.2,
          x = 0.5),
        margin = list(
          l = plot_margin_left,  # left margin
          r = plot_margin_right,  # right margin
          b = plot_margin_bottom,  # bottom margin
          t = plot_margin_top   # top margin
        ))


  }

  return(filtering_stats)

}

#' Function that generates a string with filtering criteria
#' for callsets coming from tumor-only sequencing
#'
#' @param conf list object with PCGR run configuration settings
#' @return string with filtering criteria
#'
#' @export
#'
get_tumor_only_filtering_criteria <- function(conf){

  invisible(
    assertthat::assert_that(
      is.list(conf),
      msg = paste0("Argument 'conf' must be of type list"))
  )
  invisible(
    assertthat::assert_that(
      "somatic_snv" %in% names(conf),
      msg = paste0("Argument 'conf' must contain a 'somatic_snv' entry"))
  )
  invisible(
    assertthat::assert_that(
      "tumor_only" %in% names(conf$somatic_snv),
      msg = paste0("Argument 'conf$somatic_snv' must contain a ",
                   "'tumor_only' entry"))
  )

  criteria <- c("gnomAD")
  if(as.logical(
    conf$somatic_snv$tumor_only$exclude_clinvar_germline) == TRUE){
    criteria <- c(criteria, "ClinVar (germline)")
  }
  if(as.logical(
    conf$somatic_snv$tumor_only$exclude_pon) == TRUE){
    criteria <- c(criteria, "Panel of Normals")
  }
  if(as.logical(
    conf$somatic_snv$tumor_only$exclude_likely_hom_germline) == TRUE){
    criteria <- c(criteria, "Likely germline (Homozygous AF)")
  }
  if(as.logical(
    conf$somatic_snv$tumor_only$exclude_likely_het_germline) == TRUE){
    criteria <- c(criteria, "Likely germline (Heterozygous AF)")
  }
  if(as.logical(
    conf$somatic_snv$tumor_only$exclude_dbsnp_nonsomatic) == TRUE){
    criteria <- c(criteria, "dbSNP (non-somatic)")
  }

  ## don't include exonic filter here

  all_criteria <- paste(criteria, collapse = ", ")
  return(all_criteria)

}
