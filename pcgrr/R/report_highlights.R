## Key findings ("highlights") of a PCGR report, as a short templated text -
## the most critical findings (actionable alterations) first, supported by
## other findings (oncogenic SNV/InDel variants, composite genomic
## biomarkers, copy number variants, mutational signatures, germline and
## transcriptome findings). Markdown output, for use in both the PDF and HTML reports.

#' Escape markdown emphasis characters in plain text (e.g. 'p.R1450*')
#'
#' @param x character vector
#' @keywords internal
md_escape_hl <- function(x) {
  stringr::str_replace_all(x, "([\\\\*_`\\[\\]<>])", "\\\\\\1")
}

#' Deterministic ordering of gene-level records: cancer association ranks
#' (TISSUE_ASSOC_RANK, then GLOBAL_ASSOC_RANK; highest first - as in the
#' variant tables of the report), then optional extra keys, then SYMBOL
#'
#' @param d data frame with column SYMBOL (rank columns optional)
#' @param first optional leading sort keys (numeric, larger = first)
#' @param then optional sort keys applied after the ranks (larger = first)
#' @return row order (integer)
#' @keywords internal
rank_order_hl <- function(d, first = NULL, then = NULL) {
  num <- function(x) if (is.null(x)) rep(0, NROW(d)) else
    dplyr::coalesce(suppressWarnings(as.numeric(x)), 0)
  keys <- c(lapply(first, function(k) -num(k)),
            list(-num(d$TISSUE_ASSOC_RANK), -num(d$GLOBAL_ASSOC_RANK)),
            lapply(then, function(k) -num(k)),
            list(as.character(d$SYMBOL)))
  do.call(order, unname(keys))
}

#' Percentile as ordinal number ('99th'), rounded and kept within 1-99
#'
#' @param x percentile(s) (0-100)
#' @keywords internal
ordinal_hl <- function(x) {
  n <- as.integer(pmin(pmax(round(x), 1), 99))
  suffix <- ifelse(n %% 100 %in% 11:13, "th",
                   c("th", "st", "nd", "rd", rep("th", 6))[n %% 10 + 1])
  paste0(n, suffix)
}

#' Compact list of items ('A, B, C +k more')
#'
#' @param x items (character)
#' @param max_show maximum number of items shown
#' @keywords internal
fmt_list_hl <- function(x, max_show = 5) {
  x <- unique(x[!is.na(x) & x != ""])
  n <- length(x)
  if (n == 0) return("")
  s <- paste(utils::head(x, max_show), collapse = ", ")
  if (n > max_show) s <- paste0(s, " (+", n - max_show, " more)")
  s
}

#' Somatic alterations with an actionability tier (SNVs/InDels, CNAs, RNA
#' fusions) as display labels
#'
#' @param report PCGR report object
#' @param tier actionability tiers
#' @param clnsig type of biomarker evidence the tiers refer to -
#' 'therapeutic_sensitivity' (default, the tier of each variant record), or
#' e.g. 'prognostic_poor', 'prognostic_better', 'diagnostic_positive' (tiers
#' from the corresponding bm_evidence classification of each callset)
#'
#' @return data frame with LABEL, SYMBOL, ALTERATION and ACTIONABILITY_TIER
#' @keywords internal
tiered_alterations_hl <- function(report, tier = c(1, 2),
                                  clnsig = "therapeutic_sensitivity") {
  content <- report$content
  ## variant records of a callset, with ACTIONABILITY_TIER according to clnsig
  ## (fusions are matched on VAR_ID/VARIANT_CLASS only, as their ENTREZGENE
  ## differs between variant and evidence records - see
  ## prep_actble_display_tbl())
  variant_tiers <- function(category) {
    callset <- content[[category]]$callset
    v <- callset$variant
    if (clnsig == "therapeutic_sensitivity" || NROW(v) == 0) return(v)
    cls <- callset$bm_evidence[[clnsig]]$classification
    if (NROW(cls) == 0) return(v[0, , drop = FALSE])
    join_cols <- c("VAR_ID", "VARIANT_CLASS",
                   if (category != "fusion") "ENTREZGENE")
    if (!all(join_cols %in% colnames(cls)) ||
        !all(join_cols %in% colnames(v))) return(v[0, , drop = FALSE])
    dplyr::inner_join(
      dplyr::select(v, -dplyr::any_of("ACTIONABILITY_TIER")),
      dplyr::distinct(dplyr::select(
        cls, dplyr::all_of(c(join_cols, "ACTIONABILITY_TIER")))),
      by = join_cols, relationship = "many-to-many")
  }
  recs <- list()
  snv <- variant_tiers("snv_indel")
  if (isTRUE(content$snv_indel$eval) && NROW(snv) > 0 &&
      all(c("SYMBOL", "ALTERATION", "ACTIONABILITY_TIER") %in% colnames(snv))) {
    recs$snv <- data.frame(
      SYMBOL = snv$SYMBOL, ALTERATION = html_to_text(snv$ALTERATION),
      TISSUE_ASSOC_RANK = snv$TISSUE_ASSOC_RANK,
      GLOBAL_ASSOC_RANK = snv$GLOBAL_ASSOC_RANK,
      ACTIONABILITY_TIER = suppressWarnings(as.integer(snv$ACTIONABILITY_TIER)))
  }
  cna <- variant_tiers("cna")
  if (isTRUE(content$cna$eval) && NROW(cna) > 0 &&
      all(c("SYMBOL", "VARIANT_CLASS", "ACTIONABILITY_TIER") %in% colnames(cna))) {
    recs$cna <- data.frame(
      SYMBOL = cna$SYMBOL,
      TISSUE_ASSOC_RANK = cna$TISSUE_ASSOC_RANK,
      GLOBAL_ASSOC_RANK = cna$GLOBAL_ASSOC_RANK,
      ALTERATION = dplyr::case_when(
        cna$VARIANT_CLASS %in% c("homdel", "hemdel", "hetdel") ~ "loss",
        .default = as.character(cna$VARIANT_CLASS)),
      ACTIONABILITY_TIER = suppressWarnings(as.integer(cna$ACTIONABILITY_TIER)))
  }
  fus <- variant_tiers("fusion")
  if (isTRUE(content$fusion$eval) && NROW(fus) > 0) {
    f <- fusion_tier_records(fus)
    if (NROW(f) > 0) recs$fusion <- f[, c("SYMBOL", "ALTERATION",
                                          "ACTIONABILITY_TIER")]
  }
  d <- dplyr::bind_rows(recs)
  if (NROW(d) == 0) {
    return(data.frame(LABEL = character(), SYMBOL = character(),
                      ALTERATION = character(), ACTIONABILITY_TIER = integer()))
  }
  d <- d[!is.na(d$ACTIONABILITY_TIER) & d$ACTIONABILITY_TIER %in% tier, ,
         drop = FALSE]
  d <- dplyr::distinct(d)
  ## fusion records carry no cancer association ranks - absent when only
  ## fusions are tiered
  for (rank_col in c("TISSUE_ASSOC_RANK", "GLOBAL_ASSOC_RANK")) {
    if (!rank_col %in% colnames(d)) d[[rank_col]] <- rep(NA_real_, NROW(d))
  }
  ## tier (1 before 2), then cancer association ranks, then symbol
  d <- d[do.call(order, list(
    d$ACTIONABILITY_TIER,
    -dplyr::coalesce(as.numeric(d$TISSUE_ASSOC_RANK), 0),
    -dplyr::coalesce(as.numeric(d$GLOBAL_ASSOC_RANK), 0),
    d$SYMBOL, d$ALTERATION)), , drop = FALSE]
  ## (paste0() on zero rows would return a single string)
  d$LABEL <- if (NROW(d) > 0) paste0("**", md_escape_hl(d$SYMBOL), "** ",
                                     md_escape_hl(d$ALTERATION)) else
    character()
  d
}

#' Key findings of a PCGR report (highlights) as markdown
#'
#' The most critical findings come first: somatic alterations of strong or
#' potential clinical significance (AMP/ASCO/CAP Tier 1/2). These are
#' followed by supporting findings, each included only when informative:
#' prognostic and diagnostic biomarkers (Tier 1/2),
#' oncogenic/likely oncogenic SNV/InDel variants, composite genomic
#' biomarkers (TMB, MSI), copy number variants (amplified proto-oncogenes,
#' deleted tumor suppressors, genome doubling, HRD score), mutational
#' signatures, kataegis, germline findings (CPSR), and transcriptome findings
#' (RNA fusions, expression outliers).
#'
#' @param report PCGR report object
#' @param max_show maximum number of items listed per finding
#' @param links finding titles as links to the corresponding report sections
#' (PDF report - section labels as in the pcgr_quarto_report/pdf templates,
#' e.g. 'pcgr-cna'); each finding is only included when its section is part
#' of the report
#'
#' @return list with 'lead' (markdown on actionable alterations), 'pills'
#' (positive composite biomarkers), 'pill_targets' (section label per pill),
#' 'warnings' (markdown on actionable
#' variants removed by variant filtering, "" if none) and 'supporting'
#' (markdown bullet list, "" if none)
#' @export
report_highlights <- function(report, max_show = 5, links = FALSE) {
  content <- report$content
  conf <- report$settings$conf
  ## markdown link to a section of the PDF report (labelled heading)
  sec_link <- function(text, id) {
    if (isTRUE(links)) paste0("[", text, "](#", id, ")") else text
  }

  ## --- Lead: actionable alterations ---
  lead <- c()
  actionable <- tiered_alterations_hl(report, tier = c(1, 2))
  ## one list item per tier
  tier_item <- function(t, desc) {
    lbl <- actionable$LABEL[actionable$ACTIONABILITY_TIER == t]
    if (length(lbl) == 0) return(NULL)
    paste0("- **Tier ", t, "** (", desc, "): ", fmt_list_hl(lbl, max_show))
  }
  if (NROW(actionable) > 0) {
    lead <- c(lead, paste0(
      sec_link("**Actionable alterations**", "pcgr-actionability"),
      " (AMP/ASCO/CAP)\n\n",
      paste(c(tier_item(1, "strong clinical significance"),
              tier_item(2, "potential clinical significance")),
            collapse = "\n")))
  } else if (isTRUE(content$snv_indel$eval) || isTRUE(content$cna$eval) ||
             isTRUE(content$fusion$eval)) {
    lead <- c(lead, paste0(
      sec_link("**No actionable alterations**", "pcgr-actionability"),
      " of strong or potential clinical ",
      "significance (AMP/ASCO/CAP Tier 1/2) were identified."))
  }

  ## --- Supporting findings ---
  supporting <- c()

  ## prognostic and diagnostic biomarkers (AMP/ASCO/CAP Tier 1/2) - one
  ## finding per evidence type, listing alterations with their tier; only
  ## included when any such evidence is found
  bm_items <- function(types) {
    items <- c()
    for (desc in names(types)) {
      d <- tiered_alterations_hl(report, tier = c(1, 2), clnsig = types[[desc]])
      if (NROW(d) == 0) next
      d <- d[!duplicated(d$LABEL), , drop = FALSE]
      items <- c(items, paste0(desc, ": ", fmt_list_hl(
        paste0(d$LABEL, " (Tier ", d$ACTIONABILITY_TIER, ")"), max_show)))
    }
    items
  }
  progn <- bm_items(c("poor outcome" = "prognostic_poor",
                      "better outcome" = "prognostic_better"))
  if (length(progn) > 0) {
    supporting <- c(supporting, paste0(
      sec_link("**Prognostic biomarkers:**", "pcgr-diagn-progn"), " ",
      paste(progn, collapse = "; ")))
  }
  diagn <- bm_items(c("supporting diagnosis" = "diagnostic_positive"))
  if (length(diagn) > 0) {
    supporting <- c(supporting, paste0(
      sec_link("**Diagnostic biomarkers:**", "pcgr-diagn-progn"), " ",
      paste(diagn, collapse = "; ")))
  }

  ## SNV/InDel variants - number of oncogenic / likely oncogenic variants
  ## (all, including the actionable ones - their number is stated), with the
  ## first max_show of those not already listed as actionable (oncogenicity
  ## score, then cancer association ranks, then symbol)
  snv <- content$snv_indel$callset$variant
  if (isTRUE(content$snv_indel$eval) && NROW(snv) > 0 &&
      all(c("SYMBOL", "ALTERATION", "ONCOGENICITY") %in% colnames(snv))) {
    onc_class <- tolower(stringr::str_replace_all(snv$ONCOGENICITY, "_", " "))
    onc <- snv[onc_class %in% c("oncogenic", "likely oncogenic"), ,
               drop = FALSE]
    if (NROW(onc) > 0) {
      onc <- onc[rank_order_hl(onc, first = list(onc$ONCOGENICITY_SCORE)), ,
                 drop = FALSE]
      onc_alt <- html_to_text(onc$ALTERATION)
      key <- paste(onc$SYMBOL, onc_alt)
      is_actionable <- key %in% paste(actionable$SYMBOL, actionable$ALTERATION)
      labels <- paste0("**", md_escape_hl(onc$SYMBOL), "** ",
                       md_escape_hl(onc_alt))
      n_total <- length(unique(key))
      n_actionable <- length(unique(key[is_actionable]))
      other <- unique(labels[!is_actionable])
      supporting <- c(supporting, paste0(
        sec_link("**SNV/InDel variants:**", "pcgr-oncogenic-variants"),
        " **n = ", n_total,
        "** oncogenic/likely oncogenic",
        if (n_actionable > 0) paste0(
          " (", if (n_actionable == n_total) "all" else n_actionable,
          " also actionable, listed above)"),
        if (length(other) > 0) paste0(" - ", fmt_list_hl(other, max_show))))
    }
  }

  ## composite genomic biomarkers - shown as pills (below the actionable
  ## alterations) only when positive (TMB-high, MSI-high); none otherwise
  pills <- character()
  pill_targets <- character()
  if (isTRUE(content$tmb$eval) && isTRUE(content$tmb$tmb_high)) {
    est <- content$tmb$sample_estimate
    tmb_display <- conf$somatic_snv$tmb$tmb_display
    if (NROW(est) > 0 && "TMB_MEASURE" %in% colnames(est)) {
      est <- est[est$TMB_MEASURE == paste0("TMB_", tmb_display), , drop = FALSE]
    }
    pills <- c(pills, if (NROW(est) > 0) paste0(
      "TMB-high (", round(est$TMB_ESTIMATE[1], 1), " mut/Mb)") else "TMB-high")
    pill_targets <- c(pill_targets, "pcgr-tmb")
  }
  ## MSI-high only for high-confidence predictions
  if (isTRUE(content$msi$eval) && !isTRUE(content$msi$missing_data) &&
      msi_status_is_high(content$msi, high_confidence_only = TRUE)) {
    pills <- c(pills, "MSI-high")
    pill_targets <- c(pill_targets, "pcgr-msi")
  }

  ## copy number features
  if (isTRUE(content$cna$eval)) {
    cna <- content$cna$callset$variant
    cn <- c()
    if (NROW(cna) > 0 && all(c("SYMBOL", "VARIANT_CLASS") %in% colnames(cna))) {
      ranked <- function(d) unique(d$SYMBOL[rank_order_hl(d)])
      amp_onc <- ranked(cna[cna$VARIANT_CLASS == "amplification" &
                              cna$ONCOGENE %in% TRUE, , drop = FALSE])
      del_tsg <- ranked(cna[cna$VARIANT_CLASS %in% c("homdel", "hemdel") &
                              cna$TUMOR_SUPPRESSOR %in% TRUE, , drop = FALSE])
      if (length(unique(amp_onc)) > 0) {
        cn <- c(cn, paste0("amplified proto-oncogenes: ",
                           fmt_list_hl(amp_onc, max_show)))
      }
      if (length(unique(del_tsg)) > 0) {
        cn <- c(cn, paste0("deep deletions of tumor suppressors: ",
                           fmt_list_hl(del_tsg, max_show)))
      }
    }
    cna_conf <- conf$somatic_cna
    if (isTRUE(as.logical(cna_conf$genome_doubled))) {
      cn <- c(cn, "whole-genome doubling")
    }
    if (isTRUE(as.logical(cna_conf$estimate_hrd)) && !is.null(cna_conf$hrd_sum)) {
      cn <- c(cn, paste0("HRD score: ", cna_conf$hrd_sum))
    }
    if (length(cn) > 0) {
      supporting <- c(supporting, paste0(
        sec_link("**Copy number variants:**", "pcgr-cna"), " ",
        paste(cn, collapse = "; ")))
    }
  }

  ## mutational signatures - the (at most two) proposed aetiology groups
  ## with the largest contributions, each contributing at least 10%
  ## ('Other' - a mix of unrelated signatures - is not listed)
  msig <- content$mutational_signatures
  if (isTRUE(msig$eval) && !isTRUE(msig$missing_data)) {
    pg <- msig$result$contributions$per_group
    if (NROW(pg) > 0 && all(c("group", "prop_group") %in% colnames(pg))) {
      pg <- pg[!pg$group %in% "Other" & pg$prop_group >= 0.1, , drop = FALSE]
      pg <- utils::head(pg[order(-pg$prop_group, as.character(pg$group)), , drop = FALSE], 2)
      if (NROW(pg) > 0) {
        groups <- ifelse(pg$group == "Unknown", "unknown-aetiology",
                         md_escape_hl(pg$group))
        supporting <- c(supporting, paste0(
          sec_link("**Mutational signatures:**", "pcgr-signatures"),
          " largest contributions from ",
          paste0(groups, " (", round(100 * pg$prop_group), "%)",
                 collapse = " and "), " signatures"))
      }
    }
  }
  if (isTRUE(content$kataegis$eval) && NROW(content$kataegis$events) > 0) {
    supporting <- c(supporting, paste0(
      sec_link("**Kataegis:**", "pcgr-kataegis"), " ",
      NROW(content$kataegis$events),
      " potential event(s) of localized hypermutation"))
  }

  ## germline findings (CPSR)
  if (isTRUE(content$germline_classified$eval)) {
    gv <- content$germline_classified$callset$variant_display
    plp <- if (NROW(gv) > 0 && "CLINICAL_SIGNIFICANCE" %in% colnames(gv))
      gv[gv$CLINICAL_SIGNIFICANCE %in% c("Pathogenic", "Likely Pathogenic"), ,
         drop = FALSE] else data.frame()
    ## pathogenic before likely pathogenic, then cancer association ranks
    if (NROW(plp) > 0) {
      plp <- plp[rank_order_hl(
        plp, first = list(plp$CLINICAL_SIGNIFICANCE == "Pathogenic")), ,
        drop = FALSE]
    }
    supporting <- c(supporting, paste0(
      sec_link("**Germline findings (CPSR):**", "pcgr-germline"), " ",
      if (NROW(plp) > 0) paste0(
        "Pathogenic/likely pathogenic variants in cancer predisposition genes: ",
        fmt_list_hl(paste0("**", md_escape_hl(plp$SYMBOL), "** ",
                           md_escape_hl(html_to_text(plp$ALTERATION))),
                    max_show)) else
          "No pathogenic/likely pathogenic variants in cancer predisposition genes"))
  }

  ## transcriptome findings - fusions with proto-oncogene partners and
  ## expression outliers in cancer genes
  rna <- c()
  if (isTRUE(content$fusion$eval)) {
    fd <- content$fusion$callset$variant_display
    if (NROW(fd) > 0) {
      ## fusions ordered by number of split reads (highest first), then name
      fd <- fd[order(-dplyr::coalesce(as.numeric(fd$SPLIT_READS), 0),
                     as.character(fd$FUSION_GENE)), , drop = FALSE]
      onc_fus <- html_to_text(fd$FUSION_GENE2[
        fd$ONCOGENE_5P %in% TRUE | fd$ONCOGENE_3P %in% TRUE])
      rna <- c(rna, paste0(
        NROW(fd), " RNA fusion(s)",
        if (length(onc_fus) > 0) paste0(
          ", with proto-oncogene partners: ", fmt_list_hl(onc_fus, max_show))))
    }
  }
  if (isTRUE(content$expression$eval) && NROW(content$expression$outliers) > 0) {
    out <- content$expression$outliers
    ## cancer association ranks first, then |Z-score|
    top <- function(d) {
      d <- d[rank_order_hl(d, then = list(abs(d$Z_SCORE))), , drop = FALSE]
      paste0(d$SYMBOL, " (", ordinal_hl(d$PERCENTILE), " pctile)")
    }
    high <- out[out$EXPR_OUTLIER %in% "Increased expression" &
                  out$ONCOGENE %in% TRUE, , drop = FALSE]
    low <- out[out$EXPR_OUTLIER %in% "Reduced expression" &
                 out$TUMOR_SUPPRESSOR %in% TRUE, , drop = FALSE]
    if (NROW(high) > 0) {
      rna <- c(rna, paste0("overexpressed proto-oncogenes: ",
                           fmt_list_hl(top(high), max_show)))
    }
    if (NROW(low) > 0) {
      rna <- c(rna, paste0("underexpressed tumor suppressors: ",
                           fmt_list_hl(top(low), max_show)))
    }
  }
  if (length(rna) > 0) {
    supporting <- c(supporting, paste0(
      sec_link("**Transcriptome findings:**", "pcgr-rna"), " ",
      paste(rna, collapse = "; ")))
  }

  ## actionable variants/mutation hotspots removed by variant filtering
  ## (depth/allelic support; tumor-only germline filtering)
  warnings <- character()
  if (isTRUE(content$snv_indel$eval)) {
    warnings <- filtering_warnings(filtering_summary(report), max_show)
  }

  list(
    lead = paste(lead, collapse = "\n\n"),
    pills = pills,
    pill_targets = pill_targets,
    warnings = if (length(warnings) > 0) paste0(
      "**Filtered actionable variants** (not considered in this report - ",
      "see ", sec_link("*Variant filtering*", "pcgr-variant-filtering"),
      ")\n\n",
      paste0("- ", warnings, collapse = "\n")) else "",
    supporting = if (length(supporting) > 0)
      paste0("- ", supporting, collapse = "\n") else "")
}
