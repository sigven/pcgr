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

    if (exome_all_nas == TRUE | (exome_max == 0 & genome_max > 0)) {
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
  ## user-defined panel-of-normals VCF ("PANEL_OF_NORMALS" == TRUE)
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

  ## Ensure all STATUS columns used below exist — each is only created
  ## conditionally upstream (when the prerequisite source column is present)
  status_cols <- c(
    "gnomAD_AF_ABOVE_TOLERATED",
    "STATUS_CLINVAR_GERMLINE",
    "STATUS_PON",
    "STATUS_GERMLINE_HOM",
    "STATUS_GERMLINE_HET",
    "STATUS_DBSNP",
    "STATUS_TCGA",
    "STATUS_COSMIC"
  )
  for (col in status_cols) {
    if (!(col %in% colnames(sample_calls))) {
      sample_calls[[col]] <- FALSE
    }
  }

  log4r_info("Applying variant filters on tumor-only calls - assigning somatic classification")
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
      -dplyr::any_of(c(
        "STATUS_TCGA",
        "STATUS_COSMIC",
        "STATUS_DBSNP",
        "STATUS_GERMLINE_HET",
        "STATUS_GERMLINE_HOM",
        "STATUS_PON",
        "gnomAD_AF_ABOVE_TOLERATED",
        "STATUS_CLINVAR_GERMLINE")))

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
    assign_germline_popfreq_status(
      sample_calls,
      max_tolerated_af =
        tumor_only_settings[["gnomad_popmax_af_tolerated"]])

  # for (pop in c("GLOBAL", "NFE", "AMR", "AFR",
  #               "SAS", "EAS", "ASJ", "FIN", "OTH")) {
  #   sample_calls <-
  #     assign_germline_popfreq_status_old(
  #       sample_calls,
  #       pop = pop,
  #       dbquery = c("gnomADe","gnomADg"),
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
assign_germline_popfreq_status_old <- function(sample_calls,
                                           pop = "NFE",
                                           dbquery = c("gnomADe","gnomADg"),
                                           max_tolerated_af = 0.01) {

  if (pop == "GLOBAL") {
    pop = ""
  }

  for (db in dbquery) {
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




#' Function that sets gnomAD_AF_ABOVE_TOLERATED to TRUE for variants
#' if any gnomAD population frequency exceeds max_tolerated_af
#'
#' @param sample_calls data frame with variants
#' @param max_tolerated_af max tolerated germline allele frequency
#'
#' @return sample_calls
#'
#' @export
assign_germline_popfreq_status <- function(sample_calls,
                                           max_tolerated_af = 0.01) {

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
    hole_size_pie = 0.4,
    color_palette = pcgrr::color_palette) {

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
    c("SOMATIC_CLASSIFICATION"), only_colnames = FALSE, quiet = TRUE))

  df <- report$content$snv_indel$callset$variant_unfiltered

  filtering_stats <- list()
  filtering_stats[['df']] <- data.frame()
  filtering_stats[['plot']] <- NULL

  if (NROW(df) > 0) {

    filtering_stats[['df']] <-
      plyr::count(
        df$SOMATIC_CLASSIFICATION
      ) |>
      dplyr::arrange(dplyr::desc(.data$freq)) |>
      dplyr::rename(FILTER = .data$x, FREQUENCY = .data$freq) |>
      dplyr::mutate(FILTER = stringr::str_replace_all(
        .data$FILTER, "\\|", ", ")) |>
      dplyr::mutate(FILTER = dplyr::if_else(
        !(.data$FILTER %in% pcgrr::germline_filter_levels) &
          .data$FILTER != "SOMATIC",
        "MULTIPLE FILTERS",
        .data$FILTER
      )) |>
      dplyr::group_by(.data$FILTER) |>
      dplyr::reframe(FREQUENCY = sum(.data$FREQUENCY)) |>
      dplyr::arrange(dplyr::desc(.data$FREQUENCY)) |>
      dplyr::mutate(FILTER = factor(
        .data$FILTER, levels = pcgrr::germline_filter_levels)) |>
      dplyr::mutate(PERCENT = scales::percent(
        .data$FREQUENCY / sum(.data$FREQUENCY), accuracy = 0.1))

    germline_filters <- unique(filtering_stats[['df']]$FILTER)
    germline_filters <- c("SOMATIC", setdiff(germline_filters, "SOMATIC"))
    hex_colors_germline <- utils::head(
      color_palette$multi$values,
      length(germline_filters))

    if (length(germline_filters) == 2) {
      plot_margin_bottom <- 100
    }

    rgba_colors <- c()
    i <- 1
    while(i <= length(germline_filters)) {
      alpha <- opacity_filtered_categories
      ## full opacity for 'SOMATIC' category
      if (i == 1) {
        alpha <- 1
      }
      rgba_colors <- c(
        rgba_colors,
        hex_to_rgba(
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
    hole_size_pie = 0.4,
    color_palette = pcgrr::color_palette) {

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
    only_colnames = FALSE, quiet = TRUE)

  df <- report$content$snv_indel$callset$variant_unfiltered |>
    dplyr::filter(.data$SOMATIC_CLASSIFICATION == "SOMATIC")

  filtering_stats <- list()
  filtering_stats[['df']] <- data.frame()
  filtering_stats[['plot']] <- NULL

  if (NROW(df) > 0) {

    filtering_stats[['df']] <-
      ## Exonic status can be either of 'exonic' or 'nonexonic'
      plyr::count(
        df$EXONIC_STATUS) |>
      dplyr::arrange(dplyr::desc(.data$freq)) |>
      dplyr::rename(FILTER = .data$x,
                    FREQUENCY = .data$freq) |>
      dplyr::mutate(FILTER = toupper(
        paste0("SOMATIC - ", .data$FILTER))) |>
      dplyr::group_by(.data$FILTER) |>
      dplyr::reframe(FREQUENCY = sum(.data$FREQUENCY)) |>
      dplyr::arrange(dplyr::desc(.data$FREQUENCY)) |>
      dplyr::mutate(PERCENT = scales::percent(
        .data$FREQUENCY / sum(.data$FREQUENCY), accuracy = 0.1))

    ## if either category is zero, add a row with zero frequency
    if (!("SOMATIC - EXONIC" %in% filtering_stats[['df']]$FILTER)) {
      filtering_stats[['df']] <- dplyr::bind_rows(
        filtering_stats[['df']],
        data.frame(
          FILTER = "SOMATIC - EXONIC",
          FREQUENCY = 0,
          PERCENT = "0.0%"))
    }
    if (!("SOMATIC - NONEXONIC" %in% filtering_stats[['df']]$FILTER)) {
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
      hex_to_rgba(
        color_palette$multi$values[1], alpha = 1),
      hex_to_rgba(
        color_palette$multi$values[2],
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
get_tumor_only_filtering_criteria <- function(conf) {

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
  if (as.logical(
    conf$somatic_snv$tumor_only$exclude_clinvar_germline) == TRUE) {
    criteria <- c(criteria, "ClinVar (germline)")
  }
  if (as.logical(
    conf$somatic_snv$tumor_only$exclude_pon) == TRUE) {
    criteria <- c(criteria, "Panel of Normals")
  }
  if (as.logical(
    conf$somatic_snv$tumor_only$exclude_likely_hom_germline) == TRUE) {
    criteria <- c(criteria, "Likely germline (Homozygous AF)")
  }
  if (as.logical(
    conf$somatic_snv$tumor_only$exclude_likely_het_germline) == TRUE) {
    criteria <- c(criteria, "Likely germline (Heterozygous AF)")
  }
  if (as.logical(
    conf$somatic_snv$tumor_only$exclude_dbsnp_nonsomatic) == TRUE) {
    criteria <- c(criteria, "dbSNP (non-somatic)")
  }

  ## don't include exonic filter here

  all_criteria <- paste(criteria, collapse = ", ")
  return(all_criteria)

}
