#' Get documentation string for loss-of-function annotation
#'
#' @return A documentation string
#'
#' @export
#'
lof_doc_note <- function() {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "lof.md", package = "pcgrr")
  lof_template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  seqontology_url <- "http://www.sequenceontology.org/miso/current_svn/term/"

  return(glue::glue(lof_template, seqontology_url = seqontology_url))

}

#' Get documentation string for AMP/ASCO/CAP clinical actionability tiers
#'
#' @return A documentation string
#'
#' @export
actionability_doc_note <- function() {

  cgi_url <- "https://www.cancergenomeinterpreter.org/2021/biomarkers"
  civic_url <- "https://civicdb.org"
  civic_docs_url <- "https://civic.readthedocs.io/en/latest/"
  oncokb_url <- "https://oncokb.org"
  oncokb_levels_url <- "https://www.oncokb.org/therapeutic-levels"

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "actionability.md", package = "pcgrr")
  actionability_template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  return(
    glue::glue(
      actionability_template,
      civic_url = civic_url,
      cgi_url = cgi_url,
      civic_docs_url = civic_docs_url,
      oncokb_url = oncokb_url,
      oncokb_levels_url = oncokb_levels_url
    )
  )

}

#' Get documentation string for oncogenicity annotation
#'
#' @return A documentation string
#'
#' @export
#'
oncogenicity_doc_note <- function() {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "oncogenicity.md", package = "pcgrr")
  oncogenicity_template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  #clingen_url <- "https://clinicalgenome.org/"
  #vicc_url <- "https://vicc.org/"

  return(
    glue::glue(
      oncogenicity_template
    )
  )

}

#' Get documentation string for tumor mutational burden (TMB) annotation
#'
#' @return A documentation string
#'
#' @export
#'
tmb_doc_note <- function() {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "tmb.md", package = "pcgrr")
  tmb_template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  pubmed_url <- "https://pubmed.ncbi.nlm.nih.gov/"

  return(
    glue::glue(
      tmb_template, pubmed_url = pubmed_url
    )
  )

}

#' Get documentation string for mutational signatures analysis
#'
#' @return A documentation string
#'
#' @export
#'
mutational_signatures_doc_note <- function() {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "mutational_signatures.md", package = "pcgrr")
  mutsig_template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  pubmed_url <- "https://pubmed.ncbi.nlm.nih.gov/"

  return(
    glue::glue(
      mutsig_template, pubmed_url = pubmed_url
    )
  )

}

#' Get documentation string for MSI status prediction
#'
#' @return A documentation string
#'
#' @export
#'
msi_doc_note <- function() {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "msi.md", package = "pcgrr")
  msi_template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  pubmed_url <- "https://pubmed.ncbi.nlm.nih.gov/"
  wikipedia_url <- "https://en.wikipedia.org/wiki/"

  return(
    glue::glue(
      msi_template, pubmed_url = pubmed_url, wikipedia_url = wikipedia_url
    )
  )

}

#' Get documentation string for RNA expression analysis (outlier, similarity)
#'
#' @return A documentation string
#'
#' @export
#'
expression_doc_note <- function() {

  doc_md_file <- system.file(
    "templates", "doc_notes_md", "expression.md", package = "pcgrr")
  expr_template <- paste0(readLines(doc_md_file, warn = FALSE), collapse = "\n")

  return(
    glue::glue(
      expr_template
    )
  )

}

