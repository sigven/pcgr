# Lazy-loaded package data objects (stored as .rda files in data/) are
# available at runtime but invisible to R CMD check's static analyser.
# Declaring them here suppresses "no visible binding for global variable"
# notes without affecting runtime behaviour.
utils::globalVariables(c(
  "bm_categories",
  "bm_evidence",
  "color_palette",
  "cosmic_sbs_signatures",
  "data_coltype_defs",
  "effect_prediction_algos",
  "exonic_filter_levels",
  "germline_filter_levels",
  "immune_celltypes",
  "table_display_cols",
  "tcga_cohorts",
  "tsv_cols",
  "tumor_sites",
  "variant_db_url"
))
