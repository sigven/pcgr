## Helpers for the PDF report (Quarto/Typst output, tinytable tables) -
## layout and components mirror the CPSR PDF report, colors follow the PCGR
## report theme (color_palette$report_color)

#' Main color of the PCGR report (tumor-control vs. tumor-only)
#'
#' @param report PCGR report object
#' @param color_palette PCGR color palette
#'
#' @return hex color (character)
#' @keywords internal
pdf_report_color <- function(report, color_palette = pcgrr::color_palette) {
  tumor_only <- isTRUE(as.logical(
    report$settings$conf$assay_properties$vcf_tumor_only))
  color_palette$report_color$values[ifelse(tumor_only, 2, 1)]
}

#' Style the header row of a PDF report table
#'
#' Dark header (main report color, darkened) with white, bold text, white
#' separator lines between the header columns, and light zebra striping of
#' the body rows
#'
#' @param tt tinytable object
#' @param color main report color (hex), see pdf_report_color()
#'
#' @return styled tinytable object
#' @keywords internal
style_tt_header_pdf <- function(tt, color = pcgrr::color_palette$report_color$values[1]) {
  header_col <- grDevices::colorRampPalette(c(color, "black"))(10)[4]
  stripe_col <- grDevices::colorRampPalette(c(color, "white"))(20)[19]
  tt <- tinytable::style_tt(
    tt, i = 0, background = header_col, color = "white", bold = TRUE)
  n_cols <- ncol(tt@data)
  n_rows <- nrow(tt@data)
  if (n_rows >= 2) {
    tt <- tinytable::style_tt(
      tt, i = seq(2, n_rows, by = 2), background = stripe_col)
  }
  if (n_cols > 1) {
    tt <- tinytable::style_tt(
      tt, i = 0, j = seq_len(n_cols - 1), line = "r",
      line_color = "white", line_width = 0.15)
  }
  tt
}

#' Escape Typst markup characters in plain text
#'
#' @param x character vector
#' @keywords internal
typst_escape <- function(x) {
  stringr::str_replace_all(x, "([\\\\#*_\\[\\]<>@$`])", "\\\\\\1")
}

#' Character vector as Typst string literal content (escapes backslashes
#' and double quotes), for values passed to Typst functions
#'
#' @param x character vector
#' @keywords internal
typst_str <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- stringr::str_replace_all(x, stringr::fixed("\\"), "\\\\")
  stringr::str_replace_all(x, stringr::fixed('"'), '\\"')
}

#' Colored pill (rounded box with white, bold text) in a PDF report table
#' cell - raw Typst markup; the column must not be escaped by format_tt()
#'
#' @param label pill text(s)
#' @param fill background color(s) (hex)
#'
#' @return character vector with Typst markup ("" for empty labels)
#' @keywords internal
pill_pdf <- function(label, fill) {
  esc <- typst_escape(label)
  ifelse(
    is.na(label) | label == "", "",
    paste0('#box(fill: rgb("', fill, '"), radius: 3pt, ',
           'inset: (x: 4pt, y: 2.5pt))[#text(fill: white, weight: "bold", ',
           'size: 0.92em)[', esc, ']]'))
}

#' Actionability tier as colored pills in a PDF report table - colors as in
#' the HTML report (color_palette$tier_sensitivity / tier_resistance)
#'
#' @param tiers actionability tier per table row (1, 2, 3)
#' @param clnsig 'therapeutic_sensitivity' or 'therapeutic_resistance'
#'
#' @return character vector with Typst markup (see pill_pdf)
#' @keywords internal
tier_pill_pdf <- function(tiers, clnsig = "therapeutic_sensitivity") {
  pal <- if (clnsig == "therapeutic_resistance") {
    pcgrr::color_palette$tier_resistance
  } else pcgrr::color_palette$tier_sensitivity
  idx <- match(as.character(tiers), as.character(pal$levels))
  pill_pdf(ifelse(is.na(tiers), NA, paste0("T", tiers)),
           ifelse(is.na(idx), "#999999", pal$values[idx]))
}

#' Oncogenicity classification as colored pills in a PDF report table -
#' colors as in the HTML report (color_palette$oncogenicity)
#'
#' @param oncogenicity classification per table row (e.g. 'Oncogenic',
#' 'Likely_Oncogenic', 'VUS')
#'
#' @return character vector with Typst markup (see pill_pdf); unknown
#' categories in grey
#' @keywords internal
oncogenicity_pill_pdf <- function(oncogenicity) {
  pal <- pcgrr::color_palette$oncogenicity
  lbl <- stringr::str_replace_all(oncogenicity, "_", " ")
  idx <- match(tolower(lbl), tolower(pal$levels))
  pill_pdf(lbl, ifelse(is.na(idx), "#999999", pal$values[idx]))
}

#' Text as a link in a PDF report table cell (raw Typst markup; the column
#' must not be escaped by format_tt())
#'
#' @param label link text(s)
#' @param url URL(s); NA for plain (escaped) text
#'
#' @return character vector with Typst markup
#' @keywords internal
link_pdf <- function(label, url) {
  esc <- typst_escape(label)
  ifelse(is.na(url) | is.na(label) | label == "", esc,
         paste0('#link("', url, '")[', esc, ']'))
}

#' Gene symbols as links to the Open Targets Platform target page (by
#' Ensembl gene id) in a PDF report table cell (raw Typst markup; the column
#' must not be escaped by format_tt()). Symbols without a valid Ensembl gene
#' id are shown as plain (escaped) text.
#'
#' @param symbol gene symbol(s)
#' @param ensembl_gene_id Ensembl gene id(s) (e.g. ENSG00000133703)
#'
#' @return character vector with Typst markup
#' @keywords internal
gene_link_pdf <- function(symbol, ensembl_gene_id = NULL) {
  symbol <- dplyr::coalesce(as.character(symbol), "")
  if (is.null(ensembl_gene_id)) ensembl_gene_id <- rep(NA_character_, length(symbol))
  id <- stringr::str_extract(as.character(ensembl_gene_id), "ENSG[0-9]+")
  link_pdf(symbol, ifelse(
    is.na(id), NA_character_,
    paste0("https://platform.opentargets.org/target/", id)))
}

#' Strip HTML markup from annotation values of the display tables (e.g.
#' hyperlinked gene symbols) - plain text for the PDF report
#'
#' @param x character vector
#' @keywords internal
strip_html_pdf <- function(x) {
  x <- stringr::str_replace_all(as.character(x), "<[^>]+>", "")
  stringr::str_squish(x)
}

#' Convert the HTML elements of a documentation note (markdown with inline
#' HTML, as used in the HTML report) to markdown, for the PDF report
#'
#' @param x documentation note (character)
#'
#' @return documentation note with links (<a>) and <b>/<i> as markdown, and
#' line breaks (<br>) removed
#' @keywords internal
doc_note_pdf <- function(x) {
  x |>
    ## HTML links as markdown links
    stringr::str_replace_all(
      "(?s)<a\\s+href=['\"]([^'\"]+)['\"][^>]*>(.*?)</a>", "[\\2](\\1)") |>
    stringr::str_replace_all("(?s)<b[^>]*>(.*?)</b>", "**\\1**") |>
    stringr::str_replace_all("(?s)<i>(.*?)</i>", "*\\1*") |>
    stringr::str_replace_all("<br\\s*/?>", "")
}

#' Message box for a PDF report section without findings
#'
#' @param text message (plain text; Typst markup characters are escaped)
#' @param title short title shown in bold before the message
#'
#' @return raw Typst block (character), to be printed with cat() in a chunk
#' with output 'asis'
#' @keywords internal
no_findings_pdf <- function(text, title = "No findings") {
  paste0(
    "\n```{=typst}\n",
    "#block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt, ",
    "fill: luma(245), stroke: 0.6pt + luma(200), above: 0.8em, below: 1em)[",
    "#text(weight: \"bold\")[", typst_escape(title), "] #h(0.4em) ",
    "#text(fill: luma(70))[", typst_escape(text), "]]\n",
    "```\n\n")
}

#' Placeholder box for PDF report content that is not yet implemented
#' (report skeleton)
#'
#' @param text description of the content to be added (plain text)
#'
#' @return raw Typst block (character), to be printed with cat() in a chunk
#' with output 'asis'
#' @keywords internal
placeholder_pdf <- function(text) {
  paste0(
    "\n```{=typst}\n#pcgr-placeholder[", typst_escape(text), "]\n```\n\n")
}

#' Gene symbols (alphabetical) as a compact string ("A · B", "A · B +k")
#'
#' @param symbols gene symbols
#' @param max_show maximum number of symbols shown
#' @keywords internal
fmt_gene_symbols_pdf <- function(symbols, max_show = 3) {
  g <- sort(unique(symbols[!is.na(symbols) & symbols != ""]))
  n <- length(g)
  val <- paste(utils::head(g, max_show), collapse = " · ")
  if (n > max_show) val <- paste0(val, " +", n - max_show)
  val
}

#' Variant allelic fraction as a horizontal bar with the value (Typst
#' markup), as in the HTML report tables
#'
#' @param vaf allelic fractions (numeric or character, 0-1)
#' @param width bar width (Typst length)
#'
#' @return character vector with Typst markup ("" for missing values)
#' @keywords internal
vaf_bar_pdf <- function(vaf, width = "1.3cm") {
  v <- suppressWarnings(as.numeric(vaf))
  ifelse(
    is.na(v), "",
    paste0(
      '#box(grid(columns: 2, column-gutter: 0.35em, align: horizon, ',
      'box(width: ', width, ', height: 0.6em, fill: luma(228), radius: 1pt, ',
      'clip: true, box(width: ', round(100 * pmin(pmax(v, 0), 1), 1),
      '%, height: 100%, fill: rgb("#5D5165"))), [', sprintf("%.2f", v), ']))'))
}

#' Variants located in a cancer mutation hotspot (cancerhotspots.org), as
#' flagged in the HTML report (MUTATION_HOTSPOT_CANCERTYPE)
#'
#' @param variants data frame with variants
#'
#' @return logical vector (one element per row; FALSE if the hotspot column
#' is missing)
#' @keywords internal
is_hotspot_pdf <- function(variants) {
  h <- variants[["MUTATION_HOTSPOT_CANCERTYPE"]]
  if (is.null(h)) return(rep(FALSE, NROW(variants)))
  !is.na(h) & nzchar(trimws(as.character(h)))
}

#' Cancer hotspot symbol (cancerhotspots.org logo, as in the HTML report) as
#' Typst markup - image copied with the templates to the rendering
#' directory (typst/logo_cancerhotspots.png)
#'
#' @param height symbol height (Typst length)
#' @keywords internal
hotspot_symbol_pdf <- function(height = "1em") {
  paste0('#box(baseline: 20%, image("typst/logo_cancerhotspots.png", height: ',
         height, '))')
}

#' Legend of the cancer hotspot symbol (raw Typst block), printed below a
#' variant table when at least one variant is located in a hotspot
#'
#' @param variants data frame with the variants listed in the table
#'
#' @return character (raw Typst block; "" if no variant is in a hotspot)
#' @keywords internal
hotspot_legend_pdf <- function(variants) {
  if (!any(is_hotspot_pdf(variants))) return("")
  paste0("\n```{=typst}\n#block(above: 0.4em, below: 0.8em)[",
         "#set text(size: 8.5pt, fill: pcgr-muted)\n",
         hotspot_symbol_pdf("1.2em"),
         " Variant located in a cancer mutation hotspot (",
         "#link(\"https://www.cancerhotspots.org\")[cancerhotspots.org])",
         "]\n```\n\n")
}

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
variant_table_pdf <- function(variants,
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
                    ~ dplyr::coalesce(strip_html_pdf(.x), "")))
  ## gene symbols linked to Open Targets (Ensembl gene id)
  j_raw <- c()
  if ("Gene" %in% colnames(tbl)) {
    tbl[["Gene"]] <- gene_link_pdf(tbl[["Gene"]], variants[["ENSEMBL_GENE_ID"]])
    j_raw <- c(j_raw, which(colnames(tbl) == "Gene"))
  }
  ## variant class not needed where implicit (e.g. SNV/InDel-only tables)
  if (!variant_class) tbl[["Variant class"]] <- NULL
  ## OncoKB oncogenicity (e.g. CNA table), only when annotated
  okb <- if ("ONCOGENICITY_OKB" %in% colnames(variants))
    dplyr::coalesce(as.character(variants$ONCOGENICITY_OKB), "") else NULL
  if (oncokb && !is.null(okb) && any(okb != "" & okb != ".")) {
    okb[okb == "."] <- ""
    tbl[["OncoKB"]] <- oncogenicity_pill_pdf(okb)
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
    tbl[["VAF (tumor)"]] <- vaf_bar_pdf(tbl[["VAF (tumor)"]])
    j_raw <- c(j_raw, which(colnames(tbl) == "VAF (tumor)"))
  }
  if ("Tier" %in% colnames(tbl)) {
    tbl[["Tier"]] <- tier_pill_pdf(
      suppressWarnings(as.integer(tbl[["Tier"]])), clnsig = clnsig)
    j_raw <- c(j_raw, which(colnames(tbl) == "Tier"))
  }
  if ("Oncogenicity" %in% colnames(tbl)) {
    ## oncogenicity pill, followed by the cancer hotspot symbol for variants
    ## located in a mutation hotspot (as in the HTML report)
    ## (pill and symbol in a grid, vertically centered - also when the pill
    ## text wraps onto two lines)
    onc_pill <- oncogenicity_pill_pdf(tbl[["Oncogenicity"]])
    tbl[["Oncogenicity"]] <- ifelse(
      is_hotspot_pdf(variants) & onc_pill != "",
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
    style_tt_header_pdf(color = color) |>
    tinytable::style_tt(j = 1, bold = TRUE)
}

#' Sample metadata of the PDF report cover page (pre-formatted strings)
#'
#' @param report PCGR report object
#' @param color_palette PCGR color palette
#'
#' @return named list (cover metadata)
#' @keywords internal
pdf_cover_values <- function(report, color_palette = pcgrr::color_palette) {
  conf <- report$settings$conf
  na_str <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x %in% c("NA", "None")) {
      "N/A"
    } else as.character(x)
  }

  ## value with its source (e.g. '0.75 (provided)', '2.1 (estimated)')
  with_source <- function(x, src) {
    if (na_str(x) == "N/A") return("N/A")
    if (na_str(src) == "N/A") na_str(x) else paste0(x, " (", src, ")")
  }

  genome_assembly <- dplyr::case_when(
    tolower(report$settings$genome_assembly) == "grch38" ~ "GRCh38",
    tolower(report$settings$genome_assembly) == "grch37" ~ "GRCh37",
    .default = as.character(report$settings$genome_assembly))
  v <- list(
    sample_id = report$settings$sample_id,
    genome_assembly = genome_assembly,
    ## assay type, mode and genome assembly (e.g. 'WES · Tumor-control · GRCh38')
    assay = paste0(
      na_str(conf$assay_properties$type), " · ",
      if (isTRUE(as.logical(conf$assay_properties$vcf_tumor_only)))
        "Tumor-only" else "Tumor-control", " · ", genome_assembly),
    tumor_site = na_str(conf$sample_properties$site),
    purity_ploidy = paste0(
      with_source(conf$sample_properties$tumor_purity,
                  conf$sample_properties$tumor_purity_source), " / ",
      with_source(conf$sample_properties$tumor_ploidy,
                  conf$sample_properties$tumor_ploidy_source)),
    sex_age = paste0(
      na_str(conf$subject_properties$sex), " / ",
      if (na_str(conf$subject_properties$age) == "N/A") "N/A" else
        paste(conf$subject_properties$age, "years")),
    pcgr_version = na_str(report$settings$software$pcgr_version))

  v
}

#' PDF report cover page: sample metadata (pcgr-cover(), defined in
#' typst/_preamble.typ), key findings (report_highlights()) and contents
#'
#' @param report PCGR report object
#'
#' @return raw Typst block (character), to be printed with cat() in a chunk
#' with output 'asis'
#' @keywords internal
pdf_cover_typst <- function(report) {
  v <- pdf_cover_values(report)
  args <- paste0("  ", names(v), ': "', vapply(v, typst_str, character(1)),
                 '",', collapse = "\n")
  ## key findings (markdown, see report_highlights()) in a framed box
  hl <- report_highlights(report)
  paste0(
    "```{=typst}\n#pcgr-cover(\n", args, "\n)\n",
    "#pcgr-highlights[\n#pcgr-hl-lead[\n```\n\n",
    hl$lead, "\n\n",
    "```{=typst}\n",
    ## positive composite biomarkers (TMB-high, MSI-high) as pills
    if (length(hl$pills) > 0) paste0(
      "#pcgr-hl-pills((", paste0('"', typst_str(hl$pills), '",', collapse = " "),
      "))\n"),
    "]\n",
    if (nzchar(hl$supporting)) paste0(
      "#pcgr-hl-support[\n```\n\n", hl$supporting, "\n\n",
      "```{=typst}\n]\n"),
    "]\n#pcgr-cover-end()\n```\n")
}

#' Row of KPI boxes (pcgr-kpi-box(), defined in typst/_preamble.typ) in a
#' PDF report section - as the bslib value boxes of the HTML report
#'
#' @param boxes list of KPI boxes, each a list with 'title', 'value', and
#' optionally 'sub' and 'color' (hex; default: color_palette$bg_dark)
#'
#' @return raw Typst block (character), to be printed with cat() in a chunk
#' with output 'asis'
#' @keywords internal
kpi_row_pdf <- function(boxes) {
  calls <- vapply(boxes, function(b) {
    paste0(
      'pcgr-kpi-box("', typst_str(b$title), '", "', typst_str(b$value), '"',
      if (!is.null(b$sub)) paste0(', sub: "', typst_str(b$sub), '"'),
      ', color: rgb("', if (is.null(b$color)) pcgrr::color_palette$bg_dark
      else b$color, '"))')
  }, character(1))
  paste0(
    "\n```{=typst}\n#block(above: 0.6em, below: 1em, grid(\n",
    "  columns: (", paste(rep("1fr", length(boxes)), collapse = ", "), "),\n",
    "  column-gutter: 7pt,\n  ", paste(calls, collapse = ",\n  "), "))\n```\n\n")
}

#' Note box (callout) in a PDF report section - as the callouts of the HTML
#' report
#'
#' @param text note (markdown; rendered by Quarto)
#' @param title short title shown in bold above the note
#' @param kind 'note' (light, main report color) or 'warning' (amber)
#'
#' @return markdown with raw Typst block delimiters (character), to be
#' printed with cat() in a chunk with output 'asis'
#' @keywords internal
note_box_pdf <- function(text, title = "Note", kind = "note") {
  stroke <- if (kind == "warning") pcgrr::color_palette$warning else
    "pcgr-main"
  stroke <- if (kind == "warning") paste0('rgb("', stroke, '")') else stroke
  fill <- if (kind == "warning") 'rgb("#fdf5e6")' else "pcgr-stripe"
  paste0(
    "\n```{=typst}\n#block(width: 100%, inset: (x: 10pt, y: 8pt), ",
    "radius: (right: 3pt), fill: ", fill, ", ",
    "stroke: (left: 3pt + ", stroke, "), above: 0.8em, below: 1em)[\n",
    "#text(weight: \"bold\")[", typst_escape(title), "]\n\n",
    "#set text(size: 9pt)\n```\n\n",
    text,
    "\n\n```{=typst}\n]\n```\n\n")
}

#' Sort evidence levels in the order of color_palette$elevel (A-E, then
#' OncoKB levels)
#'
#' @param levels evidence level codes
#' @keywords internal
sort_elevels_pdf <- function(levels) {
  levels[order(match(levels, names(pcgrr::color_palette$elevel)), levels)]
}

#' Actionable RNA fusions as records for the PDF report tables (gene pair
#' as 'SYMBOL', as for SNVs/InDels and CNAs)
#'
#' @param fusions data frame with fusions (fusion callset 'variant')
#'
#' @return data frame with SYMBOL, ALTERATION, VARIANT_CLASS and
#' ACTIONABILITY_TIER (integer)
#' @keywords internal
fusion_tier_records_pdf <- function(fusions) {
  if (NROW(fusions) == 0 ||
      !all(c("FUSION_GENE2", "ACTIONABILITY_TIER") %in% colnames(fusions))) {
    return(data.frame())
  }
  data.frame(
    SYMBOL = as.character(fusions$FUSION_GENE2),
    ALTERATION = "fusion",
    VARIANT_CLASS = as.character(fusions$VARIANT_CLASS),
    ACTIONABILITY_TIER = suppressWarnings(
      as.integer(fusions$ACTIONABILITY_TIER)))
}

#' Biomarker evidence sources (CIViC, CGI, OncoKB) as logos (Typst markup)
#' - images copied with the templates to the rendering directory
#' (typst/logo_*.png); unknown sources as plain text
#'
#' @param sources source identifiers (e.g. 'civic', 'cgi', 'oncokb')
#' @param height logo height (Typst length)
#' @keywords internal
bm_source_logo_pdf <- function(sources, height = "1.4em") {
  ids <- tolower(sources)
  ifelse(
    ids %in% c("civic", "cgi", "oncokb"),
    paste0('#box(baseline: 20%, image("typst/logo_', ids, '.png", height: ',
           height, '))'),
    typst_escape(sources))
}

#' Legend of the biomarker evidence source logos (raw Typst block)
#'
#' @param sources source identifiers shown in the legend
#' @keywords internal
bm_source_legend_pdf <- function(sources = c("civic", "cgi", "oncokb")) {
  names <- c(civic = "CIViC", cgi = "CGI", oncokb = "OncoKB")
  items <- paste0(bm_source_logo_pdf(sources, height = "1.2em"), " ",
                  names[tolower(sources)])
  paste0("\n```{=typst}\n#block(above: 0.4em, below: 0.8em)[",
         "#set text(size: 8.5pt, fill: pcgr-muted)\n",
         "Evidence sources: ", paste(items, collapse = " #h(0.8em) "),
         "]\n```\n\n")
}

#' Evidence levels (CIViC/CGI: A-E, OncoKB: 1-4, R1/R2) as colored pills -
#' colors as in the HTML report (color_palette$elevel)
#'
#' @param levels evidence level codes (e.g. 'A', '3B')
#'
#' @return character vector with Typst markup ("" for empty levels)
#' @keywords internal
elevel_pill_pdf <- function(levels) {
  pal <- pcgrr::color_palette$elevel
  vapply(levels, function(l) {
    if (is.na(l) || l == "") return("")
    col <- if (l %in% names(pal)) pal[[l]] else list(bg = "#999999", fg = "#ffffff")
    paste0('#box(fill: rgb("', col$bg, '"), radius: 3pt, ',
           'inset: (x: 3.5pt, y: 2.5pt))[#text(fill: rgb("', col$fg,
           '"), weight: "bold", size: 0.92em)[', typst_escape(l), ']]')
  }, character(1), USE.NAMES = FALSE)
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
bm_evidence_table_pdf <- function(
    rctbl_recs,
    clnsig = "therapeutic_sensitivity",
    color = pcgrr::color_palette$report_color$values[1]) {
  recs <- rctbl_recs[["main"]]
  if (NROW(recs) == 0) {
    return(NULL)
  }
  ## evidence sources as logos (as in the HTML report), see
  ## bm_source_legend_pdf()
  sources <- vapply(strsplit(dplyr::coalesce(recs$BM_SOURCES, ""), "\\|"),
                    function(s) paste(bm_source_logo_pdf(s), collapse = " "),
                    character(1))
  ## matched therapies: '<b>therapy</b> (tumor types) | ...' - therapy in
  ## bold, one match per line
  therapy <- vapply(recs$THERAPY_MATCH, function(x) {
    x <- stringr::str_replace(x, "^ - ", "")
    m <- strsplit(x, " \\| ")[[1]]
    m <- vapply(m, function(y) {
      parts <- stringr::str_match(y, "(?s)^<b>(.*?)</b>(.*)$")
      if (is.na(parts[1, 1])) return(typst_escape(strip_html_pdf(y)))
      paste0("*", typst_escape(parts[1, 2]), "* ",
             typst_escape(strip_html_pdf(parts[1, 3])))
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
      dplyr::summarise(LEVELS = list(sort_elevels_pdf(unique(.data$LEVEL))),
                       .groups = "drop")
    idx <- match(
      paste(recs$VAR_ID, recs$ENTREZGENE, recs$ACTIONABILITY_TIER),
      paste(lvl$VAR_ID, lvl$ENTREZGENE, lvl$ACTIONABILITY_TIER))
    levels <- vapply(idx, function(i) {
      if (is.na(i)) return("")
      paste(elevel_pill_pdf(lvl$LEVELS[[i]]), collapse = " ")
    }, character(1))
  }
  tbl <- data.frame(
    Alteration = typst_escape(strip_html_pdf(recs$SAMPLE_ALTERATION)),
    Tier = tier_pill_pdf(recs$ACTIONABILITY_TIER, clnsig = clnsig),
    `Therapy (tumor types)` = therapy,
    `Evidence level` = levels,
    Sources = sources,
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.2, 0.07, 0.48, 0.12, 0.13)) |>
    style_tt_header_pdf(color = color) |>
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
bm_evidence_tables_pdf <- function(
    report, variant_category = "snv_indel", tier = c(1, 2, 3),
    color = pcgrr::color_palette$report_color$values[1]) {
  callset <- report$content[[variant_category]]$callset
  out <- list(sensitivity = NULL, resistance = NULL)
  if (is.null(callset) || is.null(callset$bm_evidence)) {
    return(out)
  }
  if (any(callset$variant_display$ACTIONABILITY_TIER %in% tier)) {
    out$sensitivity <- bm_evidence_table_pdf(
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
    out$resistance <- bm_evidence_table_pdf(
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
bm_evidence_section_pdf <- function(
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
    tbls <- bm_evidence_tables_pdf(
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
    return(no_findings_pdf(paste0(
      "Based on available data in ", bm_db_label_pdf(report),
      ", no therapeutic biomarker evidence (Tier ",
      paste(tier, collapse = "/"), ") was found for somatic aberrations ",
      "in this tumor sample.")))
  }
  paste(out, collapse = "")
}

#' Biomarker databases used for actionability classification (as named in
#' the HTML report)
#'
#' @param report PCGR report object
#'
#' @return character, e.g. 'CIViC and CGI'
#' @keywords internal
bm_db_label_pdf <- function(report) {
  conf <- report$settings$conf
  token <- conf$oncokb$api_token
  oncokb_active <- !is.null(token) && !is.na(token) && token != "None"
  if (oncokb_active && isTRUE(as.logical(conf$oncokb$exclusive))) {
    "OncoKB"
  } else if (oncokb_active) {
    "CIViC, CGI, and OncoKB"
  } else {
    "CIViC and CGI"
  }
}

#' Font settings for (static) ggplot2 figures of the PDF report - generic
#' sans-serif font family and smaller text than the HTML report plots
#' (which set 'Helvetica', not available on all systems)
#'
#' @param p ggplot2 object
#' @param base_size base text size (pt)
#'
#' @return ggplot2 object
#' @keywords internal
pdf_plot_theme <- function(p, base_size = 9) {
  p + ggplot2::theme(
    text = ggplot2::element_text(family = "sans", size = base_size),
    axis.text.x = ggplot2::element_text(family = "sans", size = base_size,
                                        face = "plain"),
    axis.text.y = ggplot2::element_text(family = "sans", size = base_size),
    axis.title.x = ggplot2::element_text(family = "sans", size = base_size),
    legend.text = ggplot2::element_text(family = "sans", size = base_size),
    plot.margin = ggplot2::margin(4, 6, 4, 4))
}

#' Donut chart of variant statistics (static, ggplot2) for the PDF report -
#' counterpart of plotly_pie_chart() in the HTML report
#'
#' @param df_variant_stats data frame with 'N', 'Pct' and the category
#' column, from stats_type_snv_indel()
#' @param category category column (e.g. 'CODING_STATUS')
#' @param title chart title
#' @param center_label text in the center of the donut (default: total N)
#' @param show_counts show counts (N) with the legend labels
#' @param color_palette PCGR color palette
#'
#' @return ggplot2 object
#' @keywords internal
pie_chart_pdf <- function(df_variant_stats,
                          category = "CODING_STATUS",
                          title = NULL,
                          center_label = NULL,
                          show_counts = TRUE,
                          color_palette = pcgrr::color_palette) {
  d <- df_variant_stats
  d$LABEL <- stringr::str_replace_all(as.character(d[[category]]), "_", " ")
  d <- d[order(-d$N), , drop = FALSE]
  ## legend labels truncated, so that they fit the fixed legend margin
  d$LABEL <- stringr::str_trunc(d$LABEL, 30)
  if (show_counts) d$LABEL <- paste0(d$LABEL, " (", d$N, ")")
  d$LABEL <- factor(d$LABEL, levels = d$LABEL)
  if (is.null(center_label)) center_label <- paste0("N = ", sum(d$N))
  d$ymax <- cumsum(d$N)
  d$ymin <- d$ymax - d$N
  d$ymid <- (d$ymax + d$ymin) / 2
  d$PCT_LABEL <- ifelse(d$Pct >= 6, paste0(round(d$Pct), "%"), "")
  cols <- utils::head(color_palette$multi$values, NROW(d))

  ggplot2::ggplot(d) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = 2, xmax = 3.4, ymin = .data$ymin,
                   ymax = .data$ymax, fill = .data$LABEL),
      colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(x = 2.7, y = .data$ymid, label = .data$PCT_LABEL),
      colour = "white", fontface = "bold", size = 2.7) +
    ggplot2::annotate(
      "text", x = 0.6, y = 0, label = center_label,
      fontface = "bold", size = 3.2, colour = color_palette$bg_dark) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(c(0.6, 3.4)) +
    ggplot2::scale_fill_manual(values = cols) +
    ## a title line is always reserved, so that all donuts (with the same
    ## figure size) are drawn at the same size
    ggplot2::labs(title = if (is.null(title)) " " else title) +
    ggplot2::theme_void(base_size = 7) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 8, colour = color_palette$bg_dark,
        margin = ggplot2::margin(b = 2)),
      ## legend placed in a fixed right margin (not next to the panel), so
      ## that the donut size does not depend on the legend width
      legend.position = "inside",
      legend.position.inside = c(1.04, 0.5),
      legend.justification.inside = c(0, 0.5),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 6.5),
      legend.key.size = grid::unit(0.28, "cm"),
      plot.margin = ggplot2::margin(2, 125, 2, 2))
}

#' Donut charts of SNV/InDel variant statistics (coding status, variant
#' type, consequence types) for the PDF report
#'
#' @param report PCGR report object
#'
#' @return list of ggplot2 objects (coding consequence chart only when
#' coding variants are present)
#' @keywords internal
varstats_plots_pdf <- function(report) {
  variant_stats <- pcgrr::stats_type_snv_indel(
    var_df = report$content$snv_indel$callset$variant_display,
    pct_other_limit = 3.5)
  plots <- list(
    coding = pie_chart_pdf(
      variant_stats[["coding"]], "CODING_STATUS", "Coding effect"),
    variant_class = pie_chart_pdf(
      variant_stats[["variant_class"]], "VARIANT_CLASS", "Variant type"),
    consequence = pie_chart_pdf(
      variant_stats[["consequence"]], "CONSEQUENCE",
      "Consequence type - all variants"))
  if (NROW(variant_stats[["consequence_coding"]]) > 0) {
    plots[["consequence_coding"]] <- pie_chart_pdf(
      variant_stats[["consequence_coding"]], "CONSEQUENCE",
      "Consequence type - coding variants")
  }
  plots
}

#' Decode escaped Unicode code points ('<U+2019>') in text - written to the
#' report data when signature annotations are read in a non-UTF-8 locale
#'
#' @param x character vector
#' @keywords internal
decode_unicode_escapes_pdf <- function(x) {
  stringr::str_replace_all(
    x, "<U\\+([0-9A-Fa-f]{4,6})>",
    function(m) intToUtf8(strtoi(stringr::str_sub(m, 4, -2), 16L)))
}

#' Signature ID as a colored pill (Typst markup)
#'
#' @param ids signature identifiers (e.g. 'SBS13')
#' @param fill background color(s) (hex)
#' @keywords internal
signature_pill_pdf <- function(ids, fill) {
  pill_pdf(ids, fill)
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
signature_table_pdf <- function(
    per_signature, color = pcgrr::color_palette$report_color$values[1]) {
  d <- per_signature
  ## first sentence of the COSMIC aetiology description
  aetiology <- decode_unicode_escapes_pdf(dplyr::coalesce(d$aetiology, ""))
  aetiology <- stringr::str_replace(aetiology, "^(.*?[.])(\\s.*)?$", "\\1")
  tbl <- data.frame(
    Signature = signature_pill_pdf(
      d$signature_id, dplyr::coalesce(d$col, pcgrr::color_palette$none)),
    Contribution = typst_escape(d$contribution),
    `95% CI` = typst_escape(paste0(
      round(d$prop_signature_ci_lower * 100, 2), " - ",
      round(d$prop_signature_ci_upper * 100, 2), "%")),
    `Aetiology group` = typst_escape(d$group),
    `Proposed aetiology (COSMIC)` = typst_escape(aetiology),
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.12, 0.12, 0.14, 0.17, 0.45)) |>
    style_tt_header_pdf(color = color) |>
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
signature_similarity_table_pdf <- function(
    similarity, n = 10,
    color = pcgrr::color_palette$report_color$values[1]) {
  d <- utils::head(
    similarity[order(-similarity$SIMILARITY), , drop = FALSE], n)
  tbl <- data.frame(
    Signature = signature_pill_pdf(
      d$SIGNATURE_ID,
      ifelse(d$SITE_SPECIFIC %in% "YES", color, pcgrr::color_palette$none)),
    `Cosine similarity` = sprintf("%.2f", d$SIMILARITY),
    `Aetiology` = typst_escape(dplyr::coalesce(d$AETIOLOGY_KEYWORD, "")),
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.2, 0.25, 0.55)) |>
    style_tt_header_pdf(color = color) |>
    tinytable::style_tt(j = 2, bold = TRUE)
}

#' Copy number plot settings (tumor ploidy and effective amplification /
#' gain / deletion thresholds) for the PDF report - as used for the copy
#' number profile plots of the HTML report
#'
#' @param report PCGR report object
#'
#' @return list with tumor_ploidy, amp_threshold_effective, threshold_mode,
#' gain_threshold_effective and del_threshold_effective
#' @keywords internal
cna_plot_settings_pdf <- function(report) {
  cna_conf <- report$settings$conf$somatic_cna
  ploidy <- report$settings$conf$sample_properties$tumor_ploidy
  ploidy <- if (is.null(ploidy) || is.na(ploidy) || ploidy == "NA") 2 else
    as.numeric(ploidy)
  amp <- cna_conf$amp_threshold_effective
  if (is.null(amp) || is.na(amp)) amp <- cna_conf$amp_threshold_absolute
  mode <- as.character(cna_conf$threshold_mode)
  relative <- identical(mode, "relative")
  list(
    tumor_ploidy = ploidy,
    amp_threshold_effective = as.numeric(amp),
    threshold_mode = mode,
    gain_threshold_effective = if (relative)
      ploidy * as.numeric(cna_conf$gain_threshold_relative) else
        as.numeric(cna_conf$gain_threshold_absolute),
    del_threshold_effective = if (relative)
      ploidy * as.numeric(cna_conf$del_threshold_relative) else
        as.numeric(cna_conf$del_threshold_absolute))
}

#' Number of Mitelman DB records of a fusion in a given tumor site, parsed
#' from the annotated Mitelman evidence ('<b>Site</b> (<i>n = 3</i>) - ...')
#'
#' @param evidence Mitelman evidence per fusion (MITDB_EVIDENCE)
#' @param site primary tumor site (e.g. 'Colon/Rectum')
#'
#' @return integer vector (0 when the site is not recorded)
#' @keywords internal
mitdb_site_count_pdf <- function(evidence, site) {
  if (is.null(site) || site %in% c("Any", "")) {
    return(rep(0L, length(evidence)))
  }
  site_re <- stringr::str_replace_all(site, "([.|()\\^{}+$*?\\[\\]\\\\])",
                                      "\\\\\\1")
  m <- stringr::str_match(
    dplyr::coalesce(evidence, ""),
    paste0("<b>", site_re, "</b>\\s*\\(<i>n\\s*=\\s*([0-9]+)</i>\\)"))
  dplyr::coalesce(as.integer(m[, 2]), 0L)
}

#' Rank RNA fusions for the PDF report: actionability tier (Tier 1 first, no
#' tier last), then recurrence in the tumor site of the sample (Mitelman
#' DB), then the number of split reads
#'
#' @param fusions data frame with fusions (fusion callset - variant_display)
#' @param site primary tumor site of the sample
#'
#' @return ordered data frame, with column SITE_RECURRENCE
#' @keywords internal
rank_fusions_pdf <- function(fusions, site = "Any") {
  d <- fusions
  if (NROW(d) == 0) return(d)
  d$SITE_RECURRENCE <- mitdb_site_count_pdf(d$MITDB_EVIDENCE, site)
  tier <- suppressWarnings(as.integer(d$ACTIONABILITY_TIER))
  d[order(dplyr::coalesce(tier, 99L), -d$SITE_RECURRENCE,
          -dplyr::coalesce(as.numeric(d$SPLIT_READS), 0)), , drop = FALSE]
}

#' RNA fusion table for the PDF report - fusion events ranked with
#' rank_fusions_pdf(), with recurrence in the tumor site of the sample and
#' overall recurrence (Mitelman DB), and actionability tier as colored pills
#'
#' @param fusions data frame with fusions (fusion callset - variant_display)
#' @param site primary tumor site of the sample
#' @param max_rows maximum number of fusions shown (after ranking)
#' @param color main report color (hex)
#'
#' @return tinytable object
#' @keywords internal
fusion_table_pdf <- function(
    fusions, site = "Any", max_rows = 25,
    color = pcgrr::color_palette$report_color$values[1]) {
  d <- utils::head(rank_fusions_pdf(fusions, site), max_rows)
  partner <- function(gene, oncogene) {
    paste0(typst_escape(strip_html_pdf(gene)),
           ifelse(oncogene %in% TRUE, " (oncogene)", ""))
  }
  tier <- suppressWarnings(as.integer(d$ACTIONABILITY_TIER))
  tbl <- data.frame(
    Fusion = typst_escape(strip_html_pdf(d$FUSION_GENE2)),
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
  tbl[["Tier"]] <- tier_pill_pdf(ifelse(tier %in% c(1, 2, 3), tier, NA))
  widths <- c(widths, 0.09, 0.07)
  tinytable::tt(tbl, width = widths / sum(widths)) |>
    style_tt_header_pdf(color = color) |>
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
expression_outlier_table_pdf <- function(
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
    Gene = gene_link_pdf(d$SYMBOL, d[["ENSEMBL_GENE_ID"]]),
    `Gene name` = typst_escape(dplyr::coalesce(strip_html_pdf(d$GENENAME), "")),
    Role = role,
    `log2(TPM)` = sprintf("%.2f", d$EXPR_LEVEL),
    Percentile = pill_pdf(sprintf("%.1f", d$PERCENTILE), pal$values[bin]),
    `Z-score` = sprintf("%.2f", d$Z_SCORE),
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.12, 0.38, 0.2, 0.1, 0.1, 0.1)) |>
    style_tt_header_pdf(color = color) |>
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
expression_similarity_table_pdf <- function(
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
    Sample = link_pdf(
      col_or_empty("EXT_SAMPLE_NAME"),
      ifelse(grepl("^https?://", col_or_empty("EXT_LINK")),
             col_or_empty("EXT_LINK"), NA_character_)),
    Correlation = sprintf("%.3f", d$CORR),
    Diagnosis = typst_escape(col_or_empty("EXT_PRIMARY_DIAGNOSIS")),
    `Primary site` = pill_pdf(
      primary_site,
      ifelse(primary_site == site, color, pcgrr::color_palette$none)),
    `Sample type` = typst_escape(format_sample_type(col_or_empty("EXT_SAMPLE_TYPE"))),
    check.names = FALSE)
  tinytable::tt(tbl, width = c(0.2, 0.12, 0.3, 0.2, 0.18)) |>
    style_tt_header_pdf(color = color) |>
    tinytable::style_tt(j = 2, bold = TRUE)
}

#' Germline variant classification as colored pills in a PDF report table -
#' pathogenicity colors as in the HTML report (color_palette$pathogenicity)
#'
#' @param classifications classification per table row (e.g. 'Pathogenic',
#' 'Likely Pathogenic', 'VUS')
#'
#' @return character vector with Typst markup (see pill_pdf); unknown
#' categories in grey
#' @keywords internal
pathogenicity_pill_pdf <- function(classifications) {
  pal <- pcgrr::color_palette$pathogenicity
  idx <- match(tolower(classifications), tolower(pal$levels))
  pill_pdf(classifications, ifelse(is.na(idx), "#999999", pal$values[idx]))
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
germline_table_pdf <- function(
    variants, low_depth_threshold = 10,
    color = pcgrr::color_palette$report_color$values[1]) {
  d <- variants
  col_or_empty <- function(x) {
    if (x %in% colnames(d)) dplyr::coalesce(strip_html_pdf(d[[x]]), "") else
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
  genotype_md <- pill_pdf(genotype, gt_fill)
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
    Gene = gene_link_pdf(col_or_empty("SYMBOL"), d[["ENSEMBL_GENE_ID"]]),
    Alteration = typst_escape(col_or_empty("ALTERATION")),
    Consequence = typst_escape(consequence),
    Genotype = genotype_md,
    `Clinical significance` = pathogenicity_pill_pdf(
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
    style_tt_header_pdf(color = color) |>
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
kataegis_table_pdf <- function(
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
    style_tt_header_pdf(color = color) |>
    tinytable::style_tt(j = 1, bold = TRUE)
}
