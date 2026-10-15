## PDF report (Quarto/Typst) - report color and cover page (sample
## metadata, key findings, contents)

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
