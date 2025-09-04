#' Function that assigns a maximum value to a variable
#' (gnomAD_NC_FAF_GRPMAX) reflecting the filter allele frequency
#' (GrpMAX) for a given variant in the non-cancer
#' gnomAD subset (v3.1)
#'
#' @param sample_calls data frame with sample variant calls
#'
#' @export
grpmax_faf_nc_gnomad <- function(sample_calls) {
  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )

  if ("gnomAD_NC_FAF_GRPMAX" %in% colnames(sample_calls)) {
    sample_calls$gnomAD_NC_FAF_GRPMAX <- NULL
  }

  # Identify exome and genome AF columns
  af_cols <-
    grep("^tmp_gNC_FAF_(AFR|AMR|NFE|SAS|EAS|GLOBAL)$",
         names(sample_calls), value = TRUE)

  # Function to compute max AF per row
  get_max_allele_freq_nc <- function(row) {
    af_nc_vals <- as.numeric(row[af_cols])
    af_max <- 0
    if (all(is.na(af_nc_vals))) {
      return(0)
    } else {
      af_max <- max(af_nc_vals, na.rm = TRUE)
    }

    return(af_max)

  }

  # Apply to each row
  sample_calls$gnomAD_NC_FAF_GRPMAX <-
    apply(sample_calls, 1, get_max_allele_freq_nc)
  return(sample_calls)

}



#' Function that assigns a maximum value to a variable
#' (gnomAD_AF_POPMAX) reflecting the maximum allele frequency
#' for a given variant across gnomAD populations
#'
#' @param sample_calls data frame with sample variant calls
#'
#' @export
popmax_af_gnomad <- function(sample_calls) {
  invisible(
    assertthat::assert_that(is.data.frame(sample_calls),
                            msg = paste0("Argument 'sample_calls' must be of ",
                                         "type data.frame"))
  )

  if ("gnomAD_AF_POPMAX" %in% colnames(sample_calls)) {
    sample_calls$gnomAD_AF_POPMAX <- NULL
  }

  # Identify exome and genome AF columns
  exome_cols <-
    grep("^gnomADe_((NFE|AFR|AMR|SAS|EAS|FIN)_)?AF",
         names(sample_calls), value = TRUE)
  genome_cols <-
    grep("^gnomADg_((NFE|AFR|AMR|SAS|EAS|FIN)_)?AF",
         names(sample_calls), value = TRUE)

  # Function to compute max AF per row
  get_max_allele_freq <- function(row) {
    exome_vals <- as.numeric(row[exome_cols])
    genome_vals <- as.numeric(row[genome_cols])

    retmax <- 0
    genome_max <- 0
    exome_max <- 0
    exome_all_nas <- FALSE

    ## prefer exome AF values over genome AF values (if both are present)
    ## reason: 1) larger sample size in exome AF values
    ## if all exome AF values are NA (or zero and genome non-zero),
    ## then use genome AF values
    if (all(is.na(exome_vals))) {
      if (all(is.na(genome_vals))) {
        #return(NA_real_)
        return(0)
      } else {
        exome_all_nas <- TRUE
        genome_max <- max(genome_vals, na.rm = TRUE)
        #if (genome_max > 0) {
        retmax <- genome_max
        #}
      }
    } else {
      exome_max <- max(exome_vals, na.rm = TRUE)
      if (!all(is.na(genome_vals))) {
        genome_max <- max(genome_vals, na.rm = TRUE)
      }
      retmax <- exome_max
    }

    if(exome_all_nas == TRUE | (exome_max == 0 & genome_max > 0)){
      retmax <- genome_max
    }
    return(retmax)

  }

  # Apply to each row
  sample_calls$gnomAD_AF_POPMAX <-
    apply(sample_calls, 1, get_max_allele_freq)
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
      "POPMAX_AF_GNOMAD" %in% colnames(sample_calls) &
      "STATUS_COSMIC" %in% colnames(sample_calls) &
      "STATUS_TCGA" %in% colnames(sample_calls)) {
    sample_calls <- sample_calls |>
      dplyr::mutate(
        STATUS_GERMLINE_HET =
          dplyr::if_else(
            !is.na(.data$POPMAX_AF_GNOMAD) &
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
  ##    greater than the defined threshold by the user
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
          .data$gnomAD_AF_ABOVE_TOLERATED == TRUE,
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
        "gnomAD_AF_ABOVE_TOLERATED",
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

  sample_calls <-
    pcgrr::assign_germline_popfreq_status(
      sample_calls,
      max_af =
        tumor_only_settings[["gnomad_popmax_af_tolerated"]])
        #tumor_only_settings[[paste0("maf_gnomad_", tolower(pop))]])

  ## assign STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED
  # for (pop in c("GLOBAL", "NFE", "AMR", "AFR",
  #               "SAS", "EAS", "ASJ", "FIN", "REMAINING")) {
  #   sample_calls <-
  #     pcgrr::assign_germline_popfreq_status(
  #       sample_calls,
  #       pop = pop,
  #       dbquery = "gnomADe",
  #       max_tolerated_af =
  #         tumor_only_settings[[paste0("maf_gnomad_", tolower(pop))]])
  # }

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

#' Function that sets gnomAD_AF_ABOVE_TOLERATED to TRUE for variants
#' if any gnomAD population frequency exceeds max_tolerated_af
#'
#' @param sample_calls data frame with variants
#' @param max_af max tolerated germline allele frequency
#'
#' @return sample_calls
#'
#' @export
assign_germline_popfreq_status <- function(sample_calls,
                                           max_af = 0.01) {

  if (!("gnomAD_AF_ABOVE_TOLERATED" %in% colnames(sample_calls))) {
      sample_calls$gnomAD_AF_ABOVE_TOLERATED <- FALSE
  }

  # Identify exome and genome AF columns
  exome_cols <-
    grep("^gnomADe_((NFE|AFR|AMR|SAS|EAS|FIN)_)?AF",
         names(sample_calls), value = TRUE)
  genome_cols <-
    grep("^gnomADg_((NFE|AFR|AMR|SAS|EAS|FIN)_)?AF",
         names(sample_calls), value = TRUE)


  # Function to check max_tolerated_af per row
  af_exceeds_threshold <- function(row) {
    exome_vals <- as.numeric(row[exome_cols])
    genome_vals <- as.numeric(row[genome_cols])

    ## prefer exome AF values over genome AF values (if both are present)
    ## reason: larger sample size in exome AF values
    if (!all(is.na(exome_vals))) {
      return(any(exome_vals > max_tolerated_af, na.rm = TRUE))
    } else if (!all(is.na(genome_vals))) {
      return(any(genome_vals > max_tolerated_af, na.rm = TRUE))
    } else {
      return(FALSE)  # No AF data available
    }
  }

  # Apply to each row
  sample_calls$gnomAD_AF_ABOVE_TOLERATED <-
    apply(sample_calls, 1, af_exceeds_threshold)
  return(sample_calls)
}

  # if (dbquery == "gnomADe") {
  #   if (!("gnomADe_AF_ABOVE_TOLERATED" %in% colnames(sample_calls))) {
  #     sample_calls$gnomADe_AF_ABOVE_TOLERATED <- FALSE
  #   }
  #   col <- paste0(dbquery,"_",pop, "_AF")
  #   if (any(grepl(paste0("^", col, "$"), names(sample_calls)))) {
  #
  #     sample_calls$max_tolerated_af <- max_tolerated_af
  #
  #     if (NROW(
  #       sample_calls[!is.na(sample_calls[, col]) &
  #                    sample_calls[, col] > sample_calls$max_tolerated_af, ]) > 0) {
  #       sample_calls[!is.na(sample_calls[, col]) &
  #                      sample_calls[, col] > sample_calls$max_tolerated_af,
  #                    "gnomADe_AF_ABOVE_TOLERATED"] <- TRUE
  #     }
  #     sample_calls$max_tolerated_af <- NULL
  #   }
  # }

  #return(sample_calls)
#}


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
