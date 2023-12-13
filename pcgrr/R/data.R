#' List of URLS and variant identifiers for variant/gene/protein domain databases
#'
#'
#' @format A data.frame with 6 rows and 5 columns that indicates URL's for various variant/gene databases
#' and how to use PCGR annotation columns to generate variant links
#' \itemize{
#'   \item \emph{name} - Name encoding for variant/gene database
#'   \item \emph{group_by_var} - Which column should be used for grouping
#'   \item \emph{url_prefix} - URL prefix
#'   \item \emph{link_key_var} - Which column to be used as the key value in link
#'   \item \emph{link_display_var} - Which column to be used as the display variable in link
#' }
#'
"variant_db_url"

#' List of URLs for a range of variant effect prediction algorithms
#'
#'
#' @format A data.frame with 21 rows and 3 columns that indicates URL's for
#' variant effect prediction algorithms
#' \itemize{
#'   \item \emph{algorithm} - Name encoding for effect prediction algorithm
#'   \item \emph{url} - URL
#'   \item \emph{display_name} - Display name for use in reporting
#' }
#'
"effect_prediction_algos"

#' Regular expression of terms indicative of cancer-related phenotypes and syndromes
#'
#' @format A long regular expression of cancer-related phenotype terms
#'
"cancer_phenotypes_regex"

#' Color encodings for report elements of PCGR/CPSR
#'
#' @format A list object with different report elements that are color-coded in
#' PCGR/CPSR reports. Each list element have two vectors: 'levels' and 'values'.
#' Currently, the following list elements are included:
#' \itemize{
#'   \item \emph{pathogenicity} - Colors for five-level pathogenicity levels (CPSR)
#'   \item \emph{clinical_evidence} - Colors for strength of evidence of cancer-variant associations (A-E)
#'   \item \emph{tier} - Colors for tier levels for variant prioritization (PCGR)
#'   \item \emph{report_color} - Colors for PCGR assay mode (tumor-control vs. tumor-only)
#'   \item \emph{warning} - Color for warning (low confidence in PCGR analysis output)
#'   \item \emph{success} - Color for success (no evident uncertainty in PCGR analysis output)
#' }
#'
"color_palette"

