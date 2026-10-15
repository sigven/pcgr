## PDF report (Quarto/Typst) - inline markup and block elements: Typst
## escaping, table header styling, pills, links, bars, logos, legends,
## note boxes and KPI box rows. Layout mirrors the CPSR PDF report;
## colors follow the PCGR report theme (color_palette$report_color)

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
pdf_style_table <- function(tt, color = pcgrr::color_palette$report_color$values[1]) {
  header_col <- grDevices::colorRampPalette(c(color, "black"))(10)[4]
  stripe_col <- grDevices::colorRampPalette(c(color, "white"))(20)[19]
  ## cell contents centered vertically (rows with pills, logos or bars are
  ## taller than plain text rows)
  tt <- tinytable::style_tt(tt, align = "l", alignv = "m")
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

#' Colored pill (rounded box with bold text) in a PDF report table cell -
#' raw Typst markup; the column must not be escaped by format_tt()
#'
#' @param label pill text(s)
#' @param fill background color(s) (hex)
#' @param text_color text color(s) (hex)
#'
#' @return character vector with Typst markup ("" for empty labels)
#' @keywords internal
pdf_pill <- function(label, fill, text_color = "#ffffff") {
  esc <- typst_escape(label)
  ifelse(
    is.na(label) | label == "", "",
    paste0('#box(fill: rgb("', fill, '"), radius: 3pt, ',
           'inset: (x: 4pt, y: 2.5pt))[#text(fill: rgb("', text_color,
           '"), weight: "bold", size: 0.92em)[', esc, ']]'))
}

#' Classification values as colored pills in a PDF report table - colors
#' as in the HTML report
#'
#' @param x values per table row: actionability tier (1-3) for 'tier',
#' classification for 'oncogenicity' (e.g. 'Likely_Oncogenic') and
#' 'pathogenicity' (e.g. 'Likely Pathogenic'), evidence level code for
#' 'elevel' (e.g. 'A', '3B')
#' @param type 'tier' (color_palette$tier_sensitivity / tier_resistance),
#' 'oncogenicity', 'pathogenicity' or 'elevel' (color_palette$elevel)
#' @param clnsig for type 'tier': 'therapeutic_sensitivity' or
#' 'therapeutic_resistance'
#'
#' @return character vector with Typst markup (see pdf_pill); unknown values
#' in grey
#' @keywords internal
pdf_pill_class <- function(
    x, type = c("tier", "oncogenicity", "pathogenicity", "elevel"),
    clnsig = "therapeutic_sensitivity") {
  type <- match.arg(type)
  grey <- "#999999"
  if (type == "tier") {
    pal <- if (clnsig == "therapeutic_resistance") {
      pcgrr::color_palette$tier_resistance
    } else pcgrr::color_palette$tier_sensitivity
    idx <- match(as.character(x), as.character(pal$levels))
    return(pdf_pill(ifelse(is.na(x), NA, paste0("T", x)),
                    ifelse(is.na(idx), grey, pal$values[idx])))
  }
  if (type == "elevel") {
    ## background and text color per level (dark text on light levels)
    pal <- pcgrr::color_palette$elevel
    x <- as.character(x)
    col <- lapply(x, function(l) {
      if (!is.na(l) && l %in% names(pal)) pal[[l]] else
        list(bg = grey, fg = "#ffffff")
    })
    return(pdf_pill(x, vapply(col, `[[`, "", "bg"),
                    vapply(col, `[[`, "", "fg")))
  }
  pal <- pcgrr::color_palette[[type]]
  lbl <- if (type == "oncogenicity") stringr::str_replace_all(x, "_", " ") else x
  idx <- match(tolower(lbl), tolower(pal$levels))
  pdf_pill(lbl, ifelse(is.na(idx), grey, pal$values[idx]))
}

#' Text as a link in a PDF report table cell (raw Typst markup; the column
#' must not be escaped by format_tt())
#'
#' @param label link text(s)
#' @param url URL(s); NA for plain (escaped) text
#'
#' @return character vector with Typst markup
#' @keywords internal
pdf_link <- function(label, url) {
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
pdf_gene_link <- function(symbol, ensembl_gene_id = NULL) {
  symbol <- dplyr::coalesce(as.character(symbol), "")
  if (is.null(ensembl_gene_id)) ensembl_gene_id <- rep(NA_character_, length(symbol))
  id <- stringr::str_extract(as.character(ensembl_gene_id), "ENSG[0-9]+")
  pdf_link(symbol, ifelse(
    is.na(id), NA_character_,
    paste0("https://platform.opentargets.org/target/", id)))
}

#' Variant allelic fraction as a horizontal bar with the value (Typst
#' markup), as in the HTML report tables
#'
#' @param vaf allelic fractions (numeric or character, 0-1)
#' @param width bar width (Typst length)
#'
#' @return character vector with Typst markup ("" for missing values)
#' @keywords internal
pdf_vaf_bar <- function(vaf, width = "1.3cm") {
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
is_hotspot_variant <- function(variants) {
  h <- variants[["MUTATION_HOTSPOT_CANCERTYPE"]]
  if (is.null(h)) return(rep(FALSE, NROW(variants)))
  !is.na(h) & nzchar(trimws(as.character(h)))
}

#' Logo of an annotation source (CIViC, CGI, OncoKB, cancerhotspots.org) as
#' inline Typst markup - images copied with the templates to the rendering
#' directory (typst/logo_*.png); unknown sources as plain text
#'
#' @param id source identifier(s) ('civic', 'cgi', 'oncokb',
#' 'cancerhotspots')
#' @param height logo height (Typst length)
#' @keywords internal
pdf_logo <- function(id, height = "1.2em") {
  ids <- tolower(id)
  ifelse(
    ids %in% c("civic", "cgi", "oncokb", "cancerhotspots"),
    paste0('#box(baseline: 20%, image("typst/logo_', ids, '.png", height: ',
           height, '))'),
    typst_escape(id))
}

#' Legend line printed below a PDF report table (raw Typst block): the
#' cancer hotspot symbol ('hotspot' - only when at least one of the variants
#' is located in a hotspot), or the biomarker evidence source logos
#' ('bm_sources')
#'
#' @param type 'hotspot' or 'bm_sources'
#' @param variants data frame with the variants listed in the table (for
#' 'hotspot')
#'
#' @return character (raw Typst block; "" when there is nothing to explain)
#' @keywords internal
pdf_legend <- function(type = c("hotspot", "bm_sources"), variants = NULL) {
  type <- match.arg(type)
  if (type == "hotspot") {
    if (!any(is_hotspot_variant(variants))) return("")
    body <- paste0(
      pdf_logo("cancerhotspots", "1.2em"),
      " Variant located in a cancer mutation hotspot (",
      "#link(\"https://www.cancerhotspots.org\")[cancerhotspots.org])")
  } else {
    sources <- c(civic = "CIViC", cgi = "CGI", oncokb = "OncoKB")
    body <- paste0(
      "Evidence sources: ",
      paste(pdf_logo(names(sources), "1.2em"), sources, collapse = " #h(0.8em) "))
  }
  paste0("\n```{=typst}\n#block(above: 0.4em, below: 0.8em)[",
         "#set text(size: 8.5pt, fill: pcgr-muted)\n", body, "]\n```\n\n")
}

#' Note box in a PDF report section (raw Typst block)
#'
#' @param text note text - markdown for 'note' and 'warning' (rendered by
#' Quarto), plain text for 'empty' and 'placeholder' (Typst markup
#' characters are escaped)
#' @param title short title in bold (default: 'Note', or 'No findings' for
#' 'empty'; not used for 'placeholder')
#' @param kind 'note' (callout, main report color), 'warning' (callout,
#' amber), 'empty' (grey box for a section without findings), or
#' 'placeholder' (content not yet implemented)
#'
#' @return character, to be printed with cat() in a chunk with output 'asis'
#' @keywords internal
pdf_note <- function(text, title = NULL,
                     kind = c("note", "warning", "empty", "placeholder")) {
  kind <- match.arg(kind)
  if (kind == "placeholder") {
    return(paste0(
      "\n```{=typst}\n#pcgr-placeholder[", typst_escape(text), "]\n```\n\n"))
  }
  if (kind == "empty") {
    if (is.null(title)) title <- "No findings"
    return(paste0(
      "\n```{=typst}\n",
      "#block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt, ",
      "fill: luma(245), stroke: 0.6pt + luma(200), above: 0.8em, below: 1em)[",
      "#text(weight: \"bold\")[", typst_escape(title), "] #h(0.4em) ",
      "#text(fill: luma(70))[", typst_escape(text), "]]\n",
      "```\n\n"))
  }
  if (is.null(title)) title <- "Note"
  stroke <- if (kind == "warning")
    paste0('rgb("', pcgrr::color_palette$warning, '")') else "pcgr-main"
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

#' Row of KPI boxes (pcgr-kpi-box(), defined in typst/_preamble.typ) in a
#' PDF report section - as the bslib value boxes of the HTML report
#'
#' @param boxes list of KPI boxes, each a list with 'title', 'value', and
#' optionally 'sub' and 'color' (hex; default: color_palette$bg_dark)
#'
#' @return raw Typst block (character), to be printed with cat() in a chunk
#' with output 'asis'
#' @keywords internal
pdf_kpi_row <- function(boxes) {
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
