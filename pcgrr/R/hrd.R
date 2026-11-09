#' Function that initiates report element with HRD results
#'
#' The HRD report element holds the result of an external HRD assay/tool
#' ('external', \code{--input_hrd}). The genomic scar score estimated by PCGR
#' from the copy number segments (\code{--estimate_hrd}) is kept in the
#' settings (\code{settings$conf$somatic_cna}); 'internal_eval' indicates
#' whether it is available.
#'
#' @export
init_hrd_content <- function() {
  rep <- list()
  rep[["eval"]] <- FALSE
  rep[["internal_eval"]] <- FALSE
  rep[["external"]] <- init_external_result()
  return(rep)
}

#' Broad categories of HRD scores (HRD_SCORE_TYPE of an external HRD result)
#'
#' @param score_type HRD score type, one of 'genomic_scar',
#' 'mutational_signature', 'composite' or 'other'
#' @return label for display (NA for an unknown score type)
#'
#' @export
hrd_score_type_label <- function(score_type = NULL) {
  labels <- c(
    genomic_scar = "Genomic scar score",
    mutational_signature = "Mutational signature-based classifier",
    composite = "Composite model (genomic scars and mutational signatures)",
    other = "Other")
  if (is.null(score_type) || length(score_type) != 1 || is.na(score_type) ||
      !(score_type %in% names(labels))) {
    return(NA_character_)
  }
  unname(labels[score_type])
}

#' Function that loads the normalised result of an external HRD assay/tool
#'
#' Reads the (single-row) tab-separated file with the HRD result, as written
#' by the PCGR Python workflow (columns SAMPLE_ID, HRD_STATUS,
#' HRD_STATUS_BASIS, HRD_CALLER, CALLER_VERSION, HRD_SCORE_TYPE,
#' HRD_SCORE_DESCRIPTION, HRD_SCORE, HRD_SCORE_THRESHOLD_USED, HRD_LOH,
#' HRD_LST, HRD_TAI).
#'
#' @param fname name of the file with the normalised external HRD result
#' @return data frame with one row
#'
#' @export
load_external_hrd <- function(fname = NULL) {
  load_external_result(
    fname,
    required_cols = c("HRD_STATUS", "HRD_STATUS_BASIS", "HRD_CALLER",
                      "HRD_SCORE_TYPE"),
    num_cols = c("HRD_SCORE", "HRD_SCORE_THRESHOLD_USED"),
    int_cols = c("HRD_LOH", "HRD_LST", "HRD_TAI"),
    label = "HRD")
}

#' Function that generates report data for an external HRD result
#'
#' @param settings PCGR run/configuration settings (the normalised external
#' HRD result is given by \code{settings$molecular_data$fname_hrd_tsv})
#' @return list with elements 'eval', 'result' (one-row data frame) and
#' 'threshold_user_defined' (logical: the HRD status was derived from the
#' score with a user-defined threshold, \code{--input_hrd_threshold})
#'
#' @export
generate_report_data_hrd_external <- function(settings = NULL) {
  generate_report_data_external(settings, type = "hrd")
}

#' Function that checks whether the tumor is HR-deficient according to an
#' external HRD result
#'
#' @param hrd_content HRD report element
#' @return logical
#'
#' @export
hrd_status_is_deficient <- function(hrd_content = NULL) {
  isTRUE(hrd_content[["external"]][["eval"]]) &&
    identical(hrd_content[["external"]][["result"]][["HRD_STATUS"]],
              "HR-deficient")
}

#' Function that prepares display strings for an external HRD result
#'
#' @param hrd_content HRD report element
#' @return list with elements: 'status', 'caller', 'caller_label' (caller,
#' with the threshold on the score if the status is derived), 'score_type',
#' 'score_type_label', 'score_title', 'score_value', 'components' (named
#' vector with the genomic scar components given, possibly empty),
#' 'basis_note' (how the status was obtained) and 'internal_note' (NA unless
#' the PCGR genomic scar score is also available)
#'
#' @export
hrd_external_display <- function(hrd_content = NULL) {
  res <- hrd_content[["external"]][["result"]]
  caller <- res$HRD_CALLER
  if (!is.na(res$CALLER_VERSION)) {
    caller <- paste0(caller, " ", res$CALLER_VERSION)
  }
  score_type_label <- hrd_score_type_label(res$HRD_SCORE_TYPE)

  score_title <- "HRD score"
  score_kind <- tolower(score_type_label)
  if (!is.na(res$HRD_SCORE_DESCRIPTION)) {
    score_kind <- paste0(score_kind, ": ", res$HRD_SCORE_DESCRIPTION)
  }
  score_value <- "Not available"
  if (!is.na(res$HRD_SCORE)) {
    score_value <- as.character(res$HRD_SCORE)
  }

  components <- c(
    "HRD-LOH" = res$HRD_LOH, "LST" = res$HRD_LST, "TAI" = res$HRD_TAI)
  components <- components[!is.na(components)]

  caller_label <- caller
  threshold_text <- NA_character_
  if (!is.na(res$HRD_SCORE_THRESHOLD_USED)) {
    threshold_text <- as.character(res$HRD_SCORE_THRESHOLD_USED)
  }
  if (identical(res$HRD_STATUS_BASIS, "derived")) {
    basis_note <- paste0(
      "The HRD status was derived by PCGR from the score reported by ",
      caller, ", with a user-defined threshold: HR-deficient if the score is ",
      "at least ", threshold_text, ".")
    caller_label <- paste0(
      caller, " (HR-deficient if score >= ", threshold_text,
      ", user-defined threshold)")
  } else {
    basis_note <- paste0(
      "The HRD status was reported by the HRD assay/tool (", caller, ")")
    if (!is.na(threshold_text)) {
      basis_note <- paste0(
        basis_note, ", with a threshold on the score of ", threshold_text)
    }
    basis_note <- paste0(basis_note, ".")
  }
  basis_note <- paste0(
    basis_note, " The kind of score (", score_kind,
    ") was stated in the input file. Thresholds for HR deficiency are ",
    "specific to the assay/tool, and often to the tumor type - PCGR does not ",
    "assess whether the threshold used is appropriate.")

  ## the PCGR score is a genomic scar score - it is comparable only in kind
  ## with an external genomic scar score, never in scale or threshold
  internal_note <- NA_character_
  if (isTRUE(hrd_content[["internal_eval"]])) {
    if (identical(res$HRD_SCORE_TYPE, "genomic_scar")) {
      internal_note <- paste0(
        "PCGR also estimated a genomic scar score from the copy number ",
        "segments (scarHRD algorithm). Both are genomic scar scores, but they ",
        "may differ in definition, scale and calibration (e.g. sequencing ",
        "assay, segmentation, normalisation), so the threshold of the external ",
        "assay/tool does not apply to the PCGR score.")
    } else {
      internal_note <- paste0(
        "PCGR also estimated a genomic scar score from the copy number ",
        "segments (scarHRD algorithm). The external result is of a different ",
        "kind (", tolower(score_type_label), "), and the two scores are not ",
        "directly comparable.")
    }
  }

  list(
    status = res$HRD_STATUS,
    caller = caller,
    caller_label = caller_label,
    score_type = res$HRD_SCORE_TYPE,
    score_type_label = score_type_label,
    score_title = score_title,
    score_value = score_value,
    components = components,
    basis_note = basis_note,
    internal_note = internal_note)
}
