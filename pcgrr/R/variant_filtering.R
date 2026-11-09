#' Summary of SNV/InDel variant filtering, for the variant filtering
#' sections of the HTML and PDF reports
#'
#' Combines the depth/allelic support filtering (all assay modes, see
#' filter_read_support()) and - for tumor-only input - the classification of
#' variants as likely germline/somatic (SOMATIC_CLASSIFICATION) and the
#' exclusion of non-exonic variants.
#'
#' @param report PCGR report object
#'
#' @return list with elements
#' \itemize{
#'   \item tumor_only - logical
#'   \item n_input - number of variants before filtering
#'   \item n_reported - number of variants retained (reported)
#'   \item read_support - data frame with the depth/allelic support filtering
#'   steps (STEP, THRESHOLD, APPLIED, NOTE, N_REMOVED, N_REMAINING), zero rows
#'   if no thresholds were set
#'   \item read_support_removed - actionable (Tier 1/2) variants and mutation
#'   hotspots removed by depth/allelic support filtering
#'   \item germline - (tumor-only) data frame with CLASSIFICATION, N, SHARE
#'   \item exonic - (tumor-only, non-exonic variants excluded) data frame with
#'   CLASSIFICATION, N, SHARE among somatic variants
#'   \item germline_removed - (tumor-only) actionable (Tier 1/2) variants
#'   filtered as likely germline
#'   \item tumor_only_settings - tumor-only filtering settings
#' }
#' @keywords internal
filtering_summary <- function(report) {
  callset <- report$content$snv_indel$callset
  conf <- report$settings$conf
  filt <- callset$filtering
  rs <- filt$read_support
  tumor_only <- isTRUE(as.logical(conf$assay_properties$vcf_tumor_only))

  s <- list(
    tumor_only = tumor_only,
    n_input = if (!is.null(rs$n_input)) rs$n_input else NA_integer_,
    n_reported = NROW(callset$variant),
    read_support = if (NROW(rs$steps) > 0) rs$steps else data.frame(),
    read_support_removed = if (NROW(rs$actionable_removed) > 0)
      as.data.frame(rs$actionable_removed) else data.frame(),
    germline = data.frame(),
    exonic = data.frame(),
    germline_removed = if (NROW(filt$actionable_germline) > 0)
      as.data.frame(filt$actionable_germline) else data.frame(),
    tumor_only_settings = conf$somatic_snv$tumor_only)

  unf <- callset$variant_unfiltered
  if (tumor_only && NROW(unf) > 0 &&
      "SOMATIC_CLASSIFICATION" %in% colnames(unf)) {
    cl <- stringr::str_replace_all(unf$SOMATIC_CLASSIFICATION, "\\|", ", ")
    cl <- dplyr::if_else(
      !cl %in% pcgrr::germline_filter_levels, "MULTIPLE FILTERS", cl)
    s$germline <- share_table_filtering(cl, pcgrr::germline_filter_levels)

    if (isTRUE(as.logical(s$tumor_only_settings$exclude_nonexonic)) &&
        "EXONIC_STATUS" %in% colnames(unf)) {
      som <- unf[unf$SOMATIC_CLASSIFICATION == "SOMATIC", , drop = FALSE]
      if (NROW(som) > 0) {
        s$exonic <- share_table_filtering(
          toupper(paste0("SOMATIC - ", som$EXONIC_STATUS)),
          pcgrr::exonic_filter_levels)
      }
    }
  }
  s
}

#' Counts and shares of variant classifications (filtering summary)
#'
#' Filtered classifications come first, the retained classification
#' (see retained_filter_classes()) last.
#'
#' @param x character vector of classifications
#' @param levels ordering of classifications (others last)
#' @keywords internal
share_table_filtering <- function(x, levels) {
  tab <- table(x)
  d <- data.frame(CLASSIFICATION = names(tab), N = as.integer(tab))
  d <- d[order(d$CLASSIFICATION %in% retained_filter_classes(),
               match(d$CLASSIFICATION, levels), -d$N), , drop = FALSE]
  d$SHARE <- d$N / sum(d$N)
  rownames(d) <- NULL
  d
}

#' Labels of actionable variants/mutation hotspots removed by filtering
#'
#' @param d data frame with removed variants (SYMBOL, ALTERATION/HGVSP,
#' ACTIONABILITY_TIER, MUTATION_HOTSPOT, SOMATIC_CLASSIFICATION, VAF_TUMOR,
#' DP_TUMOR - all optional except SYMBOL)
#' @param markdown bold gene symbols (markdown)
#'
#' @return character vector, e.g. "**BRAF** p.V600E (Tier 1; AF 0.04, DP 12)"
#' @keywords internal
filtered_variant_labels <- function(d, markdown = TRUE) {
  if (NROW(d) == 0) return(character())
  col <- function(x) if (x %in% colnames(d)) d[[x]] else rep(NA, NROW(d))
  alt <- html_to_text(col("ALTERATION"))
  alt <- dplyr::coalesce(
    dplyr::if_else(!is.na(alt) & nzchar(alt), alt, NA_character_),
    as.character(col("HGVSP")), as.character(col("CONSEQUENCE")), "")
  sym <- dplyr::coalesce(as.character(col("SYMBOL")), "")
  tier <- suppressWarnings(as.integer(col("ACTIONABILITY_TIER")))
  hotspot <- !is.na(col("MUTATION_HOTSPOT")) &
    nzchar(as.character(dplyr::coalesce(col("MUTATION_HOTSPOT"), "")))
  af <- suppressWarnings(as.numeric(col("VAF_TUMOR")))
  dp <- suppressWarnings(as.numeric(col("DP_TUMOR")))
  cls <- as.character(col("SOMATIC_CLASSIFICATION"))

  details <- vapply(seq_len(NROW(d)), function(i) {
    a <- c(if (!is.na(tier[i]) && tier[i] <= 2) paste0("Tier ", tier[i]),
           if (hotspot[i]) "hotspot")
    b <- c(if (!is.na(af[i])) paste0("AF ", sprintf("%.2f", af[i])),
           if (!is.na(dp[i])) paste0("DP ", dp[i]),
           if (!is.na(cls[i])) stringr::str_replace_all(cls[i], "\\|", ", "))
    paste(c(paste(a, collapse = ", "), paste(b, collapse = ", "))[
      c(length(a) > 0, length(b) > 0)], collapse = "; ")
  }, character(1))
  sym_fmt <- if (markdown) paste0("**", md_escape_hl(sym), "**") else sym
  alt_fmt <- if (markdown) md_escape_hl(alt) else alt
  trimws(paste0(sym_fmt, " ", alt_fmt,
                ifelse(nzchar(details), paste0(" (", details, ")"), "")))
}

#' Threshold of a depth/allelic support filtering step, as displayed
#' (">= 10" -> "≥ 10")
#' @keywords internal
fmt_threshold_filtering <- function(x) {
  x <- stringr::str_replace(x, "^>=", "≥")
  stringr::str_replace(x, "^<=", "≤")
}

#' Classifications of variants retained by the tumor-only germline filtering
#' and exonic status filtering
#' @keywords internal
retained_filter_classes <- function() {
  c("SOMATIC", "SOMATIC - EXONIC")
}

#' Display labels/descriptions of the tumor-only germline filter categories
#'
#' @param settings tumor-only settings (conf$somatic_snv$tumor_only), used for
#' the gnomAD allele frequency threshold
#' @keywords internal
germline_filter_descriptions <- function(settings = NULL) {
  gnomad_af <- settings$gnomad_popmax_af_tolerated
  c("SOMATIC" = "Somatic (retained)",
    "GERMLINE_GNOMAD" = if (is.null(gnomad_af))
      "Germline - in gnomAD (excluded)" else
        paste0("Germline - in gnomAD with allele frequency > ", gnomad_af,
               " (excluded)"),
    "GERMLINE_CLINVAR" = "Germline - ClinVar (excluded)",
    "GERMLINE_DBSNP" = "Germline - in dbSNP, not in TCGA/COSMIC (excluded)",
    "GERMLINE_HET" = "Germline - likely heterozygous, VAF 0.4-0.6 (excluded)",
    "GERMLINE_HOM" = "Germline - likely homozygous, VAF = 1 (excluded)",
    "GERMLINE_PON" = "Germline - panel-of-normals (excluded)",
    "MULTIPLE FILTERS" = "Germline - multiple filters (excluded)",
    "SOMATIC - EXONIC" = "Somatic - exonic (retained)",
    "SOMATIC - NONEXONIC" = "Somatic - non-exonic (excluded)")
}

#' Tumor-only germline filtering settings, as a markdown bullet list
#'
#' @param settings tumor-only settings (conf$somatic_snv$tumor_only)
#' @keywords internal
tumor_only_settings_markdown <- function(settings) {
  yn <- function(x) if (isTRUE(as.logical(x))) "yes" else "no"
  paste0(
    "- gnomAD (*GERMLINE_GNOMAD*): maximum population allele frequency ",
    "tolerated for somatic variants **", settings$gnomad_popmax_af_tolerated,
    "**\n",
    "- ClinVar germline variants (*GERMLINE_CLINVAR*): **",
    yn(settings$exclude_clinvar_germline), "**\n",
    "- dbSNP variants not found in TCGA/COSMIC (*GERMLINE_DBSNP*): **",
    yn(settings$exclude_dbsnp_nonsomatic), "**\n",
    "- Likely heterozygous germline, VAF 0.4-0.6 (*GERMLINE_HET*): **",
    yn(settings$exclude_likely_het_germline), "**\n",
    "- Likely homozygous germline, VAF = 1 (*GERMLINE_HOM*): **",
    yn(settings$exclude_likely_hom_germline), "**\n",
    "- Panel-of-normals (*GERMLINE_PON*): **",
    yn(settings$exclude_pon), "**\n",
    "- Non-exonic variants excluded: **",
    yn(settings$exclude_nonexonic), "**\n")
}

#' Warning (markdown) on actionable variants/mutation hotspots removed by
#' variant filtering, for the Key findings and variant filtering sections
#'
#' @param s filtering summary, see filtering_summary()
#' @param max_show maximum number of variants listed per filtering type
#'
#' @return character vector (one element per filtering type, empty if none)
#' @keywords internal
filtering_warnings <- function(s, max_show = 5) {
  w <- character()
  n_rs <- NROW(s$read_support_removed)
  if (n_rs > 0) {
    w <- c(w, paste0(
      "n = ", n_rs, " actionable variant", if (n_rs > 1) "s" else "",
      " (Tier 1/2) or mutation hotspot", if (n_rs > 1) "s" else "",
      " removed by depth/allelic support filtering: ",
      fmt_list_hl(filtered_variant_labels(s$read_support_removed), max_show)))
  }
  n_gl <- NROW(s$germline_removed)
  if (n_gl > 0) {
    w <- c(w, paste0(
      "n = ", n_gl, " actionable variant", if (n_gl > 1) "s" else "",
      " (Tier 1/2) filtered as likely germline (tumor-only): ",
      fmt_list_hl(filtered_variant_labels(s$germline_removed), max_show)))
  }
  w
}
