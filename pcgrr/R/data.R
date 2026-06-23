#' List of URLS and variant identifiers for variant/gene/protein domain databases
#'
#'
#' @format A data.frame with 6 rows and 5 columns that indicates URL's
#' for various variant/gene databases and how to use PCGR annotation
#' columns to generate variant links
#' \itemize{
#'   \item \emph{name} - Name encoding for variant/gene database
#'   \item \emph{group_by_var} - Which column should be used for grouping
#'   \item \emph{url_prefix} - URL prefix
#'   \item \emph{link_key_var} - Which column to be used as the key value in link
#'   \item \emph{link_display_var} - Which column to be used as the display variable in link
#' }
#'
"variant_db_url"


#' Oncogenicity criteria (ClinGen/CGC/VICC)
#'
"oncogenicity_criteria"

#' Fixed data types/categories used for biomarker evidence, e.g. 'types','levels' etc.
#'
"bm_evidence"

#' List of coltype definitions for input files to pcgrr
#' (e.g. VCF-converted TSV, CNA TVS etc.)
#'
"data_coltype_defs"

#' List of COSMIC reference mutational signatures (SBS, v3.4)
#'
#' @format A list with two matrix objects ('all' and 'no_artefacts').
#' One matrix contains the COSMIC reference mutational signatures without signature
#' artefacts ('no_artefacts', number of columns = 68), while the other contains
#' all signatures, including artefacts ('all', number of columns = 86). Each
#' matrix has 96 rows, one for each of the 96 possible trinucleotide contexts.
#'
"cosmic_sbs_signatures"

#' Data frame with all TCGA cohorts
#'
#' @format A data.frame with 33 rows and 2 columns that indicates TCGA cohorts
"tcga_cohorts"


#' Data frame with immune cell types
#'
#' @format A data.frame with 11 rows and 2 columns that indicates immune
#' cell types used in immune contexture analysis by quanTIseq
#'
"immune_celltypes"

#' Data frame with germline filtering criteria
#'
#' @format A character vector listing all germline filtering criteria
#' applied on input callsets (SNVs/InDels) in tumor-only mode
#'
"germline_filter_levels"

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

#' TSV columns
"tsv_cols"

#' HTML Table display columns
"table_display_cols"

#' Exonic filter levels
"exonic_filter_levels"

#' Germline filter levels
"germline_filter_levels"

#' OncoKB base API URL
"oncokb_base_api_url"

#' Fixed data types and levels for biomarker evidence items
#'
#' @format A named list with three character vectors:
#' \itemize{
#'   \item \emph{types} - Biomarker evidence types (e.g. "predictive", "prognostic")
#'   \item \emph{levels} - Evidence strength levels (e.g. "A: Validated", "B: Clinical evidence")
#'   \item \emph{clinical_significance} - Clinical significance categories (e.g. "Sensitivity/Response")
#' }
#'
"biomarker_evidence"

#' Named list of biomarker categories used for variant classification
#'
#' @format A named list with six elements, each representing a biomarker category
#' (e.g. "therapeutic_sensitivity", "prognostic_poor"). Each element contains:
#' \itemize{
#'   \item \emph{etype} - Evidence type (e.g. "predictive", "prognostic", "diagnostic")
#'   \item \emph{clnsig} - Clinical significance label(s) for the category
#' }
#'
"bm_categories"

#' Character vector of tumor site names used in PCGR/CPSR reports
#'
#' @format A character vector of 31 tumor site names (e.g. "Any", "Adrenal Gland",
#' "Biliary Tract"). Used to constrain tumor site selection in configuration.
#'
"tumor_sites"


#' Character vector with OncoKB annotations coming from the
#' MafAnnotator / FusionAnnotator / CnaAnnnotator tools in the PCGR Python workflow.
#' These annotations are used for variant classification and reporting in PCGR.
#'
#' @format A character vector with 8 different OncoKB annotations, including:
#' "ONCOGENICITY_OKB", "MUTATION_EFFECT_OKB" etc
#'
"oncokb_annotations"
