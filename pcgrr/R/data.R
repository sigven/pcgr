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

#' VCF encodings of heterozygous variant states
#'
#' @format A vector of potential VCF encodings of heterozygous variant states
#'
"heterozygous_states"


#' VCF encodings of homozygous variant states
#'
#' @format A vector of potential VCF encodings of homozygous variant states
#'
"homozygous_states"

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

#' Scores and documentation of ACMG evidence criteria used for variant classification
#' in CPSR
#'
#' @format A list object with three elements: 'score2tier', 'evidence_codes',
#' 'pathogenic_range_gnomad'
#'
#' @format \bold{score2tier} - A data frame with 5 rows and two columns that indicate
#' current score thresholds for variant classification in CPSR:
#' \itemize{
#'   \item \emph{CPSR_CLASSIFICATION} - variant classification level (P, LP, VUS etc)
#'   \item \emph{CPSR_PATHOGENICITY_SCORE} - indication of CPSR "score bucket" for a given
#'   classification (HTML string)
#' }
#'
#' #' @format \bold{evidence_codes} - A data frame with 34 rows and 7 columns that document
#' all ACMG evidence criteria that are used for variant classification in CPSR:
#' \itemize{
#'   \item \emph{cpsr_evidence_code} - code for evidence criterion ('ACMG_BA1_AD' etc)
#'   \item \emph{category} - type of evidence feature ('clinpop','funcvarpop','funcvar','funccomp')
#'   \item \emph{pathogenicity_pole} - whether the given evidence support a benign variant character ('B'), or
#'   pathogenic character ('P')
#'   \item \emph{category_long} - long version of 'category' column
#'   \item \emph{description} - Verbose description for the given evidence criterion
#'   \item \emph{sherloc_code} - Corresponding code identifier in SherLoc (Nykamp et al., GiM, 2017)
#'   \item \emph{path_score} - Score associated with the given evidence criterion (negative or positive)
#' }
#'
"cpsr_acmg"
