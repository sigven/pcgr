## PDF report (Quarto/Typst) - table builders (tinytable objects)

#' Variant table (SNVs/InDels, or CNAs) for the PDF report - plain text
#' columns, actionability tier and oncogenicity as colored pills
#'
#' @param variants data frame with variants
#' @param color main report color (hex)
#' @param clnsig 'therapeutic_sensitivity' or 'therapeutic_resistance'
#' (tier colors)
#' @param variant_class show the variant class column
#' @param oncokb show the OncoKB oncogenicity column (ONCOGENICITY_OKB, as
#' colored pills), when OncoKB annotations are present
#'
#' @return tinytable object
#' @keywords internal
pdf_table_variants <- function(variants,
                              color = pcgrr::color_palette$report_color$values[1],
                              clnsig = "therapeutic_sensitivity",
                              variant_class = TRUE,
                              oncokb = FALSE) {
  cols <- c(
    "Gene" = "SYMBOL",
    "Alteration" = "ALTERATION",
    "Consequence" = "CONSEQUENCE",
    "Variant class" = "VARIANT_CLASS",
    "Copy number" = "CN_TOTAL",
    "Cytoband" = "CYTOBAND",
    "Tier" = "ACTIONABILITY_TIER",
    "Oncogenicity" = "ONCOGENICITY",
    "Genotype" = "GENOTYPE",
    "Clinical significance" = "CLINICAL_SIGNIFICANCE",
    "VAF (tumor)" = "VAF_TUMOR",
    "Depth (tumor)" = "DP_TUMOR")
  tbl <- variants |>
    dplyr::select(dplyr::any_of(cols)) |>
    dplyr::mutate(
      dplyr::across(dplyr::everything(),
                    ~ dplyr::coalesce(html_to_text(.x), "")))
  ## gene symbols linked to Open Targets (Ensembl gene id)
  j_raw <- c()
  if ("Gene" %in% colnames(tbl)) {
    tbl[["Gene"]] <- pdf_gene_link(tbl[["Gene"]], variants[["ENSEMBL_GENE_ID"]])
    j_raw <- c(j_raw, which(colnames(tbl) == "Gene"))
  }
  ## variant class not needed where implicit (e.g. SNV/InDel-only tables)
  if (!variant_class) tbl[["Variant class"]] <- NULL
  ## OncoKB oncogenicity (e.g. CNA table), only when annotated
  okb <- if ("ONCOGENICITY_OKB" %in% colnames(variants))
    dplyr::coalesce(as.character(variants$ONCOGENICITY_OKB), "") else NULL
  if (oncokb && !is.null(okb) && any(okb != "" & okb != ".")) {
    okb[okb == "."] <- ""
    tbl[["OncoKB"]] <- pdf_pill_class(type = "oncogenicity", okb)
    j_raw <- c(j_raw, which(colnames(tbl) == "OncoKB"))
  }
  ## consequence terms (e.g. 'splice_region_variant&stop_gained') with
  ## break opportunities, to wrap within narrow columns
  if ("Consequence" %in% colnames(tbl)) {
    tbl[["Consequence"]] <- tbl[["Consequence"]] |>
      stringr::str_replace_all("&", ", ") |>
      stringr::str_replace_all("_", " ")
  }
  if ("VAF (tumor)" %in% colnames(tbl)) {
    tbl[["VAF (tumor)"]] <- pdf_vaf_bar(tbl[["VAF (tumor)"]])
    j_raw <- c(j_raw, which(colnames(tbl) == "VAF (tumor)"))
  }
  if ("Tier" %in% colnames(tbl)) {
    tbl[["Tier"]] <- pdf_pill_class(type = "tier", 
      suppressWarnings(as.integer(tbl[["Tier"]])), clnsig = clnsig)
    j_raw <- c(j_raw, which(colnames(tbl) == "Tier"))
  }
  if ("Oncogenicity" %in% colnames(tbl)) {
    ## oncogenicity pill, followed by the cancer hotspot symbol for variants
    ## located in a mutation hotspot (as in the HTML report)
    ## (pill and symbol in a grid, vertically centered - also when the pill
    ## text wraps onto two lines)
    onc_pill <- pdf_pill_class(type = "oncogenicity", tbl[["Oncogenicity"]])
    tbl[["Oncogenicity"]] <- ifelse(
      is_hotspot_variant(variants) & onc_pill != "",
      paste0('#box(grid(columns: 2, column-gutter: 0.3em, align: horizon, [',
             onc_pill, '], image("typst/logo_cancerhotspots.png", height: 1em)))'),
      onc_pill)
    j_raw <- c(j_raw, which(colnames(tbl) == "Oncogenicity"))
  }
  tt <- tinytable::tt(tbl, width = rep(1, ncol(tbl)) / ncol(tbl))
  j_esc <- setdiff(seq_len(ncol(tbl)), j_raw)
  if (length(j_esc) > 0) {
    tt <- tinytable::format_tt(tt, j = j_esc, escape = TRUE)
  }
  tt |>
    pdf_style_table(color = color) |>
    tinytable::style_tt(j = 1, bold = TRUE)
}

#' Biomarker evidence table (therapeutic sensitivity/resistance) for the PDF
#' report - one row per matched alteration, with tier, matched therapies
#' (tumor types), and evidence sources (cf. render_actble_bm_table() of the
#' HTML report)
#'
#' @param rctbl_recs list with 'main' data frame, from
#' prep_actble_display_tbl()
#' @param clnsig 'therapeutic_sensitivity' or 'therapeutic_resistance'
#' @param color main report color (hex)
#'
#' @return tinytable object, or NULL if there are no evidence records
#' @keywords internal
pdf_table_bm_evidence <- function(
    rctbl_recs,
    clnsig = "therapeutic_sensitivity",
    color = pcgrr::color_palette$report_color$values[1]) {
  recs <- rctbl_recs[["main"]]
  if (NROW(recs) == 0) {
    return(NULL)
  }
  ## evidence sources as logos (as in the HTML report), see
  ## pdf_legend("bm_sources")
  sources <- vapply(strsplit(dplyr::coalesce(recs$BM_SOURCES, ""), "\\|"),
                    function(s) paste(pdf_logo(s, "1.4em"), collapse = " "),
                    character(1))
  ## matched therapies: '<b>therapy</b> (tumor types) | ...' - therapy in
  ## bold, one match per line
  therapy <- vapply(recs$THERAPY_MATCH, function(x) {
    x <- stringr::str_replace(x, "^ - ", "")
    m <- strsplit(x, " \\| ")[[1]]
    m <- vapply(m, function(y) {
      parts <- stringr::str_match(y, "(?s)^<b>(.*?)</b>(.*)$")
      if (is.na(parts[1, 1])) return(typst_escape(html_to_text(y)))
      paste0("*", typst_escape(parts[1, 2]), "* ",
             typst_escape(html_to_text(parts[1, 3])))
    }, character(1))
    paste(m, collapse = " #linebreak() ")
  }, character(1), USE.NAMES = FALSE)
  ## evidence levels of the matched evidence items (nested records)
  levels <- rep("", NROW(recs))
  nested <- rctbl_recs[["nested"]]
  if (NROW(nested) > 0 && "BM_EVIDENCE_LEVEL" %in% colnames(nested)) {
    lvl <- nested |>
      dplyr::mutate(LEVEL = stringr::str_replace(
        as.character(.data$BM_EVIDENCE_LEVEL), ":.*$", "")) |>
      dplyr::filter(!is.na(.data$LEVEL) & .data$LEVEL != "") |>
      dplyr::group_by(.data$VAR_ID, .data$ENTREZGENE,
                      .data$ACTIONABILITY_TIER) |>
      dplyr::summarise(LEVELS = list(sort_evidence_levels(unique(.data$LEVEL))),
                       .groups = "drop")
    idx <- match(
      paste(recs$VAR_ID, recs$ENTREZGENE, recs$ACTIONABILITY_TIER),
      paste(lvl$VAR_ID, lvl$ENTREZGENE, lvl$ACTIONABILITY_TIER))
    levels <- vapply(idx, function(i) {
      if (is.na(i)) return("")
      paste(pdf_pill_class(type = "elevel", lvl$LEVELS[[i]]), collapse = " ")
    }, character(1))
  }
  tbl <- data.frame(
    Alteration = typst_escape(html_to_text(recs$SAMPLE_ALTERATION)),
    Tier = pdf_pill_class(type = "tier", recs$ACTIONABILITY_TIER, clnsig = clnsig),
    `Therapy (tumor types)` = therapy,
    `Evidence level` = levels,
    Sources = sources,
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.2, 0.07, 0.48, 0.12, 0.13)) |>
    pdf_style_table(color = color) |>
    tinytable::style_tt(j = 1, bold = TRUE)
}

#' Biomarker evidence tables (therapeutic sensitivity: Tier 1-3, predictive
#' evidence; therapeutic resistance) of a variant category for the PDF
#' report - same evidence selection as the HTML report
#'
#' @param report PCGR report object
#' @param variant_category 'snv_indel', 'cna', 'fusion', 'msi' or 'tmb'
#' @param tier actionability tiers to include (default: 1-3)
#' @param color main report color (hex)
#'
#' @return list with 'sensitivity' and 'resistance' tinytable objects (NULL
#' when there is no evidence)
#' @keywords internal
pdf_bm_evidence <- function(
    report, variant_category = "snv_indel", tier = c(1, 2, 3),
    color = pcgrr::color_palette$report_color$values[1]) {
  callset <- report$content[[variant_category]]$callset
  out <- list(sensitivity = NULL, resistance = NULL)
  if (is.null(callset) || is.null(callset$bm_evidence)) {
    return(out)
  }
  if (any(callset$variant_display$ACTIONABILITY_TIER %in% tier)) {
    out$sensitivity <- pdf_table_bm_evidence(
      pcgrr::prep_actble_display_tbl(
        rep = report, tier = tier,
        etype_for_tiering = c("predictive"),
        clnsig = "therapeutic_sensitivity",
        variant_category = variant_category),
      clnsig = "therapeutic_sensitivity", color = color)
  }
  res_class <- callset$bm_evidence$therapeutic_resistance$classification
  if (NROW(res_class) > 0 &&
      any(res_class$ACTIONABILITY_TIER %in% tier)) {
    out$resistance <- pdf_table_bm_evidence(
      pcgrr::prep_actble_display_tbl(
        rep = report, tier = tier,
        clnsig = "therapeutic_resistance",
        variant_category = variant_category),
      clnsig = "therapeutic_resistance", color = color)
  }
  out
}

#' Biomarker evidence of Tier 1/2 aberrations (SNVs/InDels, CNAs, RNA
#' fusions) for the actionability section of the PDF report - therapeutic
#' sensitivity and resistance tables per variant category
#'
#' @param report PCGR report object
#' @param tier actionability tiers to include
#' @param color main report color (hex)
#'
#' @return markdown with raw Typst tables (character), to be printed with
#' cat() in a chunk with output 'asis'
#' @keywords internal
pdf_section_bm_evidence <- function(
    report, tier = c(1, 2),
    color = pcgrr::color_palette$report_color$values[1]) {
  categories <- c(snv_indel = "SNVs/InDels",
                  cna = "Copy number aberrations",
                  fusion = "RNA fusions")
  raw_tt <- function(x) {
    paste0("\n```{=typst}\n", tinytable::save_tt(x, output = "typst"),
           "\n```\n\n")
  }
  out <- c()
  for (cat_id in names(categories)) {
    if (!isTRUE(report$content[[cat_id]]$eval)) next
    tbls <- pdf_bm_evidence(
      report, variant_category = cat_id, tier = tier, color = color)
    if (is.null(tbls$sensitivity) && is.null(tbls$resistance)) next
    out <- c(out, paste0("### ", categories[[cat_id]], "\n\n"))
    if (!is.null(tbls$sensitivity)) {
      out <- c(out, "**Therapeutic sensitivity**\n\n", raw_tt(tbls$sensitivity))
    }
    if (!is.null(tbls$resistance)) {
      out <- c(out, "**Therapeutic resistance**\n\n", raw_tt(tbls$resistance))
    }
  }
  if (length(out) == 0) {
    return(pdf_note(kind = "empty", paste0(
      "Based on available data in ", bm_db_label(report),
      ", no therapeutic biomarker evidence (Tier ",
      paste(tier, collapse = "/"), ") was found for somatic aberrations ",
      "in this tumor sample.")))
  }
  paste(out, collapse = "")
}

#' Table of mutational signature contributions (signature refitting) for
#' the PDF report - signature (colored as in the contribution plot),
#' contribution with bootstrap confidence interval, and aetiology
#'
#' @param per_signature data frame with signature contributions (content
#' mutational_signatures - result - contributions - per_signature)
#' @param color main report color (hex)
#'
#' @return tinytable object
#' @keywords internal
pdf_table_signatures <- function(
    per_signature, color = pcgrr::color_palette$report_color$values[1]) {
  d <- per_signature
  ## first sentence of the COSMIC aetiology description
  aetiology <- decode_unicode_escapes(dplyr::coalesce(d$aetiology, ""))
  aetiology <- stringr::str_replace(aetiology, "^(.*?[.])(\\s.*)?$", "\\1")
  tbl <- data.frame(
    Signature = pdf_pill(
      d$signature_id, dplyr::coalesce(d$col, pcgrr::color_palette$none)),
    Contribution = typst_escape(d$contribution),
    `95% CI` = typst_escape(paste0(
      round(d$prop_signature_ci_lower * 100, 2), " - ",
      round(d$prop_signature_ci_upper * 100, 2), "%")),
    `Aetiology group` = typst_escape(d$group),
    `Proposed aetiology (COSMIC)` = typst_escape(aetiology),
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.12, 0.12, 0.14, 0.17, 0.45)) |>
    pdf_style_table(color = color) |>
    tinytable::style_tt(j = 2, bold = TRUE)
}

#' Table of the reference signatures most similar to the mutational spectrum
#' of the sample (cosine similarity) for the PDF report - signatures
#' previously attributed to the tumor site in the main report color
#'
#' @param similarity data frame with SIGNATURE_ID, SIMILARITY,
#' AETIOLOGY_KEYWORD and SITE_SPECIFIC (content mutational_signatures -
#' result - signature_similarity)
#' @param n maximum number of signatures shown
#' @param color main report color (hex)
#'
#' @return tinytable object
#' @keywords internal
pdf_table_signature_similarity <- function(
    similarity, n = 10,
    color = pcgrr::color_palette$report_color$values[1]) {
  d <- utils::head(
    similarity[order(-similarity$SIMILARITY), , drop = FALSE], n)
  tbl <- data.frame(
    Signature = pdf_pill(
      d$SIGNATURE_ID,
      ifelse(d$SITE_SPECIFIC %in% "YES", color, pcgrr::color_palette$none)),
    `Cosine similarity` = sprintf("%.2f", d$SIMILARITY),
    `Aetiology` = typst_escape(dplyr::coalesce(d$AETIOLOGY_KEYWORD, "")),
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.2, 0.25, 0.55)) |>
    pdf_style_table(color = color) |>
    tinytable::style_tt(j = 2, bold = TRUE)
}

#' RNA fusion table for the PDF report - fusion events ranked with
#' rank_fusions(), with recurrence in the tumor site of the sample and
#' overall recurrence (Mitelman DB), and actionability tier as colored pills
#'
#' @param fusions data frame with fusions (fusion callset - variant_display)
#' @param site primary tumor site of the sample
#' @param max_rows maximum number of fusions shown (after ranking)
#' @param color main report color (hex)
#'
#' @return tinytable object
#' @keywords internal
pdf_table_fusions <- function(
    fusions, site = "Any", max_rows = 25,
    color = pcgrr::color_palette$report_color$values[1]) {
  d <- utils::head(rank_fusions(fusions, site), max_rows)
  ## partner genes as Open Targets links - the Ensembl gene id is taken
  ## from the (HTML) link of the partner in the fusion callset
  partner <- function(gene, oncogene) {
    gene <- as.character(gene)
    paste0(pdf_gene_link(html_to_text(gene),
                         stringr::str_extract(gene, "ENSG[0-9]+")),
           ifelse(oncogene %in% TRUE, " (oncogene)", ""))
  }
  tier <- suppressWarnings(as.integer(d$ACTIONABILITY_TIER))
  tbl <- data.frame(
    Fusion = typst_escape(html_to_text(d$FUSION_GENE2)),
    `Split reads` = as.character(d$SPLIT_READS),
    `5' partner` = partner(d$FUSION_GENE_5P, d$ONCOGENE_5P),
    `3' partner` = partner(d$FUSION_GENE_3P, d$ONCOGENE_3P),
    check.names = FALSE)
  widths <- c(0.22, 0.09, 0.2, 0.2)
  ## recurrence in the tumor site of the sample (Mitelman DB)
  if (!site %in% c("Any", "")) {
    tbl[[paste0("Mitelman DB - ", typst_escape(site), " (n)")]] <-
      as.character(d$SITE_RECURRENCE)
    widths <- c(0.2, 0.08, 0.18, 0.18, 0.13)
  }
  tbl[["Mitelman DB - all (n)"]] <- as.character(dplyr::coalesce(
    as.numeric(d$MITDB_NUM_EVIDENCE), 0))
  tbl[["Tier"]] <- pdf_pill_class(type = "tier", ifelse(tier %in% c(1, 2, 3), tier, NA))
  widths <- c(widths, 0.09, 0.07)
  tinytable::tt(tbl, width = widths / sum(widths)) |>
    pdf_style_table(color = color) |>
    tinytable::style_tt(j = 1, bold = TRUE)
}

#' Gene expression outliers in cancer genes for the PDF report - increased
#' expression of proto-oncogenes, or reduced expression of tumor suppressor
#' genes; percentile as colored pills (colors as in the HTML report), ranked
#' by absolute Z-score
#'
#' @param outliers data frame with expression outliers (content - expression
#' - outliers)
#' @param direction 'Increased expression' or 'Reduced expression'
#' @param color main report color (hex)
#'
#' @return tinytable object, or NULL when there are no outliers
#' @keywords internal
pdf_table_expression_outliers <- function(
    outliers, direction = "Increased expression",
    color = pcgrr::color_palette$report_color$values[1]) {
  ## overexpressed oncogenes / underexpressed tumor suppressors
  role_col <- if (direction == "Increased expression") "ONCOGENE" else
    "TUMOR_SUPPRESSOR"
  d <- outliers[outliers$EXPR_OUTLIER %in% direction &
                  outliers[[role_col]] %in% TRUE, , drop = FALSE]
  if (NROW(d) == 0) {
    return(NULL)
  }
  d <- d[order(-abs(d$Z_SCORE)), , drop = FALSE]
  pal <- if (direction == "Increased expression")
    pcgrr::color_palette$expression_outlier_high else
      pcgrr::color_palette$expression_outlier_low
  bin <- pmin(findInterval(d$PERCENTILE, pal$breaks) + 1, length(pal$values))
  role <- dplyr::case_when(
    d$ONCOGENE %in% TRUE & d$TUMOR_SUPPRESSOR %in% TRUE ~
      "Oncogene / tumor suppressor",
    d$ONCOGENE %in% TRUE ~ "Oncogene",
    TRUE ~ "Tumor suppressor")
  tbl <- data.frame(
    Gene = pdf_gene_link(d$SYMBOL, d[["ENSEMBL_GENE_ID"]]),
    `Gene name` = typst_escape(dplyr::coalesce(html_to_text(d$GENENAME), "")),
    Role = role,
    `log2(TPM)` = sprintf("%.2f", d$EXPR_LEVEL),
    Percentile = pdf_pill(sprintf("%.1f", d$PERCENTILE), pal$values[bin]),
    `Z-score` = sprintf("%.2f", d$Z_SCORE),
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.12, 0.38, 0.2, 0.1, 0.1, 0.1)) |>
    pdf_style_table(color = color) |>
    tinytable::style_tt(j = 1, bold = TRUE)
}

#' Most similar external samples (expression correlation) of a reference
#' source for the PDF report - samples from the primary site of the query
#' sample highlighted in the main report color
#'
#' @param similarity data frame with similarity results (content -
#' expression - similarity_analysis - source)
#' @param site primary tumor site of the sample
#' @param n maximum number of samples shown
#' @param color main report color (hex)
#'
#' @return tinytable object, or NULL when there are no results
#' @keywords internal
pdf_table_expression_similarity <- function(
    similarity, site = "Any", n = 10,
    color = pcgrr::color_palette$report_color$values[1]) {
  if (NROW(similarity) == 0 || !"CORR" %in% colnames(similarity)) {
    return(NULL)
  }
  d <- utils::head(similarity[order(-similarity$CORR), , drop = FALSE], n)
  col_or_empty <- function(x) {
    if (x %in% colnames(d)) dplyr::coalesce(as.character(d[[x]]), "") else
      rep("", NROW(d))
  }
  primary_site <- col_or_empty("EXT_PRIMARY_SITE")
  tbl <- data.frame(
    ## sample name linked to the source database (e.g. DepMap), as in the
    ## HTML report
    Sample = pdf_link(
      col_or_empty("EXT_SAMPLE_NAME"),
      ifelse(grepl("^https?://", col_or_empty("EXT_LINK")),
             col_or_empty("EXT_LINK"), NA_character_)),
    Correlation = sprintf("%.3f", d$CORR),
    Diagnosis = typst_escape(col_or_empty("EXT_PRIMARY_DIAGNOSIS")),
    `Primary site` = pdf_pill(
      primary_site,
      ifelse(primary_site == site, color, pcgrr::color_palette$none)),
    `Sample type` = typst_escape(format_sample_type(col_or_empty("EXT_SAMPLE_TYPE"))),
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.2, 0.12, 0.3, 0.2, 0.18)) |>
    pdf_style_table(color = color) |>
    tinytable::style_tt(j = 2, bold = TRUE)
}

#' Germline variants in cancer predisposition genes (CPSR classification)
#' for the PDF report - genotype, clinical significance (colored pills),
#' ACMG evidence codes and mode of inheritance
#'
#' @param variants data frame with classified germline variants (content -
#' germline_classified - callset - variant_display)
#' @param low_depth_threshold variants with a control (germline) sequencing
#' depth below this value are flagged (only when depth is known)
#' @param color main report color (hex)
#'
#' @return tinytable object
#' @keywords internal
pdf_table_germline <- function(
    variants, low_depth_threshold = 10,
    color = pcgrr::color_palette$report_color$values[1]) {
  d <- variants
  col_or_empty <- function(x) {
    if (x %in% colnames(d)) dplyr::coalesce(html_to_text(d[[x]]), "") else
      rep("", NROW(d))
  }
  gt_raw <- tolower(col_or_empty("GENOTYPE"))
  genotype <- dplyr::case_when(
    gt_raw %in% c("het", "heterozygous") ~ "Heterozygous",
    gt_raw %in% c("hom", "homalt", "hom_alt", "homozygous") ~ "Homozygous",
    gt_raw %in% c("hemi", "hemizygous") ~ "Hemizygous",
    .default = col_or_empty("GENOTYPE"))
  ## genotype as colored pills - colors as in CPSR (cpsr::color_palette -
  ## genotypes): heterozygous, homozygous/hemizygous
  gt_fill <- dplyr::case_when(
    genotype == "Heterozygous" ~ "#007a74",
    genotype %in% c("Homozygous", "Hemizygous") ~ "#004c47",
    .default = "#999999")
  genotype_md <- pdf_pill(genotype, gt_fill)
  ## DP_CONTROL is -1 when CPSR was run without depth information
  if ("DP_CONTROL" %in% colnames(d)) {
    dp <- suppressWarnings(as.numeric(d$DP_CONTROL))
    low <- !is.na(dp) & dp >= 0 & dp < low_depth_threshold
    genotype_md[low] <- paste0(
      genotype_md[low], " #linebreak() #text(size: 0.9em, fill: pcgr-muted)[",
      "low depth (", dp[low], "X)]")
  }
  consequence <- col_or_empty("CONSEQUENCE") |>
    stringr::str_replace_all("&", ", ") |>
    stringr::str_replace_all("_", " ")
  tbl <- data.frame(
    Gene = pdf_gene_link(col_or_empty("SYMBOL"), d[["ENSEMBL_GENE_ID"]]),
    Alteration = typst_escape(col_or_empty("ALTERATION")),
    Consequence = typst_escape(consequence),
    Genotype = genotype_md,
    `Clinical significance` = pdf_pill_class(type = "pathogenicity", 
      col_or_empty("CLINICAL_SIGNIFICANCE")),
    `ACMG evidence` = typst_escape(
      stringr::str_replace_all(col_or_empty("ACMG_CODE"), "\\|", ", ")),
    check.names = FALSE)
  widths <- c(0.1, 0.15, 0.17, 0.15, 0.16, 0.27)
  if ("CPG_MOI" %in% colnames(d)) {
    tbl$Inheritance <- typst_escape(col_or_empty("CPG_MOI"))
    widths <- c(0.09, 0.14, 0.16, 0.14, 0.15, 0.22, 0.1)
  }
  tinytable::tt(tbl, width = widths) |>
    pdf_style_table(color = color) |>
    tinytable::style_tt(j = 1, bold = TRUE)
}

#' Kataegis events (localized hypermutation) for the PDF report - region,
#' chromosome arm, length, number of mutations, proportion of C>X mutations
#' and confidence
#'
#' @param events data frame with kataegis events (content - kataegis -
#' events)
#' @param color main report color (hex)
#'
#' @return tinytable object, or NULL when there are no events
#' @keywords internal
pdf_table_kataegis <- function(
    events, color = pcgrr::color_palette$report_color$values[1]) {
  if (NROW(events) == 0 ||
      !all(c("chrom", "start", "end") %in% colnames(events))) {
    return(NULL)
  }
  d <- as.data.frame(events)
  col_or_empty <- function(x) {
    if (x %in% colnames(d)) dplyr::coalesce(as.character(d[[x]]), "") else
      rep("", NROW(d))
  }
  weight <- suppressWarnings(as.numeric(col_or_empty("weight.C>X")))
  tbl <- data.frame(
    Region = paste0(d$chrom, ":", format(d$start, big.mark = ",", trim = TRUE),
                    "-", format(d$end, big.mark = ",", trim = TRUE)),
    `Chromosome arm` = col_or_empty("chrom.arm"),
    `Length (bp)` = format(suppressWarnings(as.numeric(col_or_empty("length"))),
                           big.mark = ",", trim = TRUE),
    `Mutations (n)` = col_or_empty("number.mut"),
    `Weight (C>X)` = ifelse(is.na(weight), "", sprintf("%.3f", weight)),
    check.names = FALSE)
  widths <- c(0.3, 0.16, 0.16, 0.14, 0.14)
  if ("confidence" %in% colnames(d)) {
    tbl$Confidence <- col_or_empty("confidence")
    widths <- c(0.28, 0.14, 0.15, 0.13, 0.13, 0.13)
  }
  tinytable::tt(tbl, width = widths) |>
    tinytable::format_tt(escape = TRUE) |>
    pdf_style_table(color = color) |>
    tinytable::style_tt(j = 1, bold = TRUE)
}

#' Table of the depth/allelic support filtering steps of SNVs/InDels for the
#' PDF report (variant filtering section)
#'
#' @param steps data frame with STEP, THRESHOLD, APPLIED, NOTE, N_REMOVED and
#' N_REMAINING (see filter_read_support())
#' @param color main report color (hex)
#'
#' @return tinytable object
#' @keywords internal
pdf_table_read_support <- function(
    steps, color = pcgrr::color_palette$report_color$values[1]) {
  applied <- as.logical(steps$APPLIED)
  tbl <- data.frame(
    Step = typst_escape(steps$STEP),
    Threshold = typst_escape(fmt_threshold_filtering(steps$THRESHOLD)),
    Removed = ifelse(
      applied, format(steps$N_REMOVED, big.mark = ",", trim = TRUE),
      paste0("#text(fill: luma(110), style: \"italic\")[",
             typst_escape(steps$NOTE), "]")),
    Remaining = format(steps$N_REMAINING, big.mark = ",", trim = TRUE),
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.36, 0.14, 0.32, 0.18)) |>
    pdf_style_table(color = color) |>
    tinytable::style_tt(j = 4, bold = TRUE)
}

#' Table of variant classifications (tumor-only germline filtering, exonic
#' status) with counts and shares (bars) for the PDF report
#'
#' @param d data frame with CLASSIFICATION, N and SHARE (see
#' filtering_summary())
#' @param retained classifications retained (in bold)
#' @param settings tumor-only settings (conf$somatic_snv$tumor_only)
#' @param color main report color (hex)
#'
#' @return tinytable object
#' @keywords internal
pdf_table_filter_shares <- function(
    d, retained = retained_filter_classes(), settings = NULL,
    color = pcgrr::color_palette$report_color$values[1]) {
  desc <- germline_filter_descriptions(settings)
  lbl <- ifelse(d$CLASSIFICATION %in% names(desc),
                desc[d$CLASSIFICATION], d$CLASSIFICATION)
  keep <- d$CLASSIFICATION %in% retained
  bar_col <- ifelse(keep, color, "#9e9e9e")
  share <- paste0(
    '#box(grid(columns: 2, column-gutter: 0.4em, align: horizon, ',
    'box(width: 2.6cm, height: 0.6em, fill: luma(228), radius: 1pt, ',
    'clip: true, box(width: ', round(100 * d$SHARE, 1),
    '%, height: 100%, fill: rgb("', bar_col, '"))), [',
    sprintf("%.1f", 100 * d$SHARE), '%]))')
  tbl <- data.frame(
    Classification = ifelse(keep, paste0("*", typst_escape(lbl), "*"),
                            typst_escape(lbl)),
    Filter = ifelse(keep, "", paste0("#text(size: 7.5pt, fill: luma(90))[",
                                     typst_escape(d$CLASSIFICATION), "]")),
    Variants = format(d$N, big.mark = ",", trim = TRUE),
    Share = share,
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.38, 0.24, 0.12, 0.26)) |>
    pdf_style_table(color = color)
}
