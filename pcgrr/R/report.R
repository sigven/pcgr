
#' Function that initiates PCGR/CPSR report object
#'
#' @param yaml_fname Yaml file with configuration settings and paths
#' to annotated molecular datasets
#' @param report_mode Type of report ('PCGR' or 'CPSR')
#' @export
init_report <- function(yaml_fname = NULL,
                        report_mode = "PCGR") {

  # invisible(assertthat::assert_that(
  #   !is.null(type), msg = "Argument 'type' must be of type character"))
  invisible(assertthat::assert_that(
    is.character(report_mode),
    msg = "Argument 'report_mode' must be of type character"))
  invisible(assertthat::assert_that(
    report_mode == "PCGR" | report_mode == "CPSR",
    msg = "Argument 'report_mode' must be of either 'PCGR' or 'CPSR'"))

  report <- list()
  for (elem in c("settings", "content","ref_data")) {
    report[[elem]] <- list()
  }

  yaml_data <- list()
  ref_data <- list()
  if (!is.null(yaml_fname)) {
    yaml_data <- pcgrr::load_yaml(
      yaml_fname,
      report_mode = report_mode
    )

    if (!is.null(yaml_data)) {
      report[["settings"]] <- yaml_data$settings
      report[["ref_data"]] <- yaml_data$ref_data
    }
  }

  if (report_mode == "CPSR") {
    report[["content"]][["snv_indel"]] <-
      init_germline_content()

    if (!is.null(report[["settings"]][["conf"]][["variant_classification"]][["pop_gnomad"]])) {
      if (report[["settings"]][["conf"]][["variant_classification"]][["pop_gnomad"]] != "") {
        population <-
          report[["settings"]][["conf"]][["variant_classification"]][["pop_gnomad"]]
        vcf_tag_AF <- "gnomADe_non_cancer_AF"
        vcf_tag_AN <- "gnomADe_non_cancer_AN"
        vcf_tag_AC <- "gnomADe_non_cancer_AC"
        vcf_tag_NHOMALT <- "gnomADe_non_cancer_NHOMALT"
        if (population != "global") {
          vcf_tag_AF <-
            paste0("gnomADe_non_cancer_",toupper(population),"_AF")
          vcf_tag_AN <-
            paste0("gnomADe_non_cancer_",toupper(population),"_AN")
          vcf_tag_AC <-
            paste0("gnomADe_non_cancer_",toupper(population),"_AC")
          vcf_tag_NHOMALT <-
            paste0("gnomADe_non_cancer_",toupper(population),"_NHOMALT")
        }
        pop_desc_df <-
          report$ref_data$vcf_infotags[
            report$ref_data$vcf_infotags$tag == vcf_tag_AF,]
        if (NROW(pop_desc_df) == 1) {
          population_description <- pop_desc_df$description
          report[["settings"]][["conf"]][["variant_classification"]][["vcftag_gnomad_AF"]] <-
            vcf_tag_AF
          report[["settings"]][["conf"]][["variant_classification"]][["vcftag_gnomad_AN"]] <-
            vcf_tag_AN
          report[["settings"]][["conf"]][["variant_classification"]][["vcftag_gnomad_AC"]] <-
            vcf_tag_AC
          report[["settings"]][["conf"]][["variant_classification"]][["vcftag_gnomad_NHOMALT"]] <-
            vcf_tag_NHOMALT
          report[["settings"]][["conf"]][["variant_classification"]][["popdesc_gnomad"]] <-
            population_description
        }
      }
    }

  }else{

    for (a_elem in c("snv_indel",
                     "cna",
                     "sample_properties",
                     "assay_properties",
                     "mutational_signatures",
                     "tmb",
                     "msi",
                     "germline_classified",
                     "rainfall",
                     "kataegis",
                     "expression")){
                     #"clinicaltrials")) {
      report[["content"]][[a_elem]] <- list()
      report[["content"]][[a_elem]][["eval"]] <- FALSE

      if (a_elem == "tumor_purity" | a_elem == "tumor_ploidy") {
        report[["content"]][[a_elem]][["eval"]] <- TRUE
      }

      if (a_elem == "kataegis") {
        report[["content"]][[a_elem]][["events"]] <- data.frame()
      }

      if (a_elem == "expression") {
        report[["content"]][[a_elem]] <-
          init_expression_content()
      }

      if (a_elem == "rainfall") {
        report[["content"]][[a_elem]] <-
          init_rainfall_content()
      }

      if (a_elem == "germline_classified") {
        report[["content"]][[a_elem]][['callset']] <- list()
        report[["content"]][[a_elem]][['panel_info']] <- list()
        report[["content"]][[a_elem]]$sample_id <- "NA"
      }

      if (a_elem == "snv_indel" | a_elem == "cna") {
        report[["content"]][[a_elem]][['callset']] <- list()

        if (a_elem == "snv_indel") {
          report[["content"]][[a_elem]][['vstats']] <-
            init_snv_indel_vstats()
        }
        if (a_elem == "cna") {
          report[["content"]][[a_elem]][['vstats']] <-
            init_cna_vstats()
        }
      }
      if (a_elem == "clinicaltrials") {
        report[["content"]][[a_elem]][["trials"]] <- data.frame()
        report[["content"]][[a_elem]][["missing_data"]] <- F
      }

      if (a_elem == "mutational_signatures") {
        report[["content"]][[a_elem]] <-
          init_m_signature_content()
      }
      if (a_elem == "tmb") {
        report[["content"]][[a_elem]] <-
          init_tmb_content(ref_data = report[["ref_data"]])
      }
      if (a_elem == "msi") {
        report[["content"]][[a_elem]][["missing_data"]] <- FALSE
        report[["content"]][[a_elem]][["prediction"]] <- list()
      }
      if (a_elem == "tumor_only") {
        report[["content"]][[a_elem]] <-
          init_tumor_only_content()

      }
    }
  }

  if(!is.null(report$settings$conf$assay_properties)){
    report[['content']][['assay_properties']] <-
      report[['settings']]$conf$assay_properties
  }

  if(!is.null(report$settings$conf$sample_properties)){
    report[['content']][['sample_properties']] <-
      report[['settings']]$conf$sample_properties
  }

  return(report)
}

#' Function that updates a PCGR/CPSR report object structure
#'
#' @param report PCGR/CPSR report
#' @param report_data Object with report contents
#' @param a_elem section of PCGR report
#' @export
update_report <- function(report, report_data,
                          a_elem = "snv_indel") {

  if (!is.null(report_data) &
      !is.null(report[["content"]][[a_elem]])) {
    for (report_elem in names(report_data)) {
      if (!is.null(report_data[[report_elem]]) &
          !is.null(report[["content"]][[a_elem]][[report_elem]])) {
         report[["content"]][[a_elem]][[report_elem]] <-
           report_data[[report_elem]]
      }
    }
  }
  return(report)
}


#' Function that initiates report element with TMB information
#'
#' @param ref_data PCGR reference data object
#
#' @return rep TMB report element
#'
#' @export
init_tmb_content <- function(ref_data = NULL) {

  invisible(assertthat::assert_that(!is.null(ref_data)))
  invisible(assertthat::assert_that(!is.null(ref_data$misc)))
  invisible(assertthat::assert_that(is.data.frame(
    ref_data$misc$tmb) &
      nrow(ref_data$misc$tmb) > 0))
  rep <- list()
  rep[["eval"]] <- FALSE
  rep[["sample_estimate"]] <- data.frame()
  rep[["tmb_reference"]] <- ref_data$misc$tmb

  return(rep)
}

#' Function that initiates report element with CNA information
#'
#' @export
init_cna_vstats <- function() {

  vstats <- list()
  for (t in c("n_tsg_loss",
              "n_oncogene_gain",
              "n_other_drugtarget_gain",
              "n_segments_loss",
              "n_segments_gain",
              "n_actionable_tier1",
              "n_actionable_tier2")) {
    vstats[[t]] <- 0
  }
  return(vstats)
}

#' Function that initiates report element with SNV/InDel statistics information
#'
#' @export
init_snv_indel_vstats <- function() {

  vstats <- list()
  for (t in c("n",
              "n_snv",
              "n_indel",
              "n_sub",
              "n_coding",
              "n_noncoding",
              "n_actionable_tier1",
              "n_actionable_tier2",
              "n_actionable_tier3",
              "n_actionable_tier3_tsg",
              "n_actionable_tier3_oncogene",
              "n_actionable_tier3_dualrole",
              "n_tier4",
              "n_tier5",
              "n_eitems_diagnostic_tier1",
              "n_eitems_predictive_tier1",
              "n_eitems_prognostic_tier1",
              "n_eitems_diagnostic_tier2",
              "n_eitems_predictive_tier2",
              "n_eitems_prognostic_tier2",
              "n_genes_tier1",
              "n_genes_tier2",
              "n_genes_tier3")) {
    vstats[[t]] <- 0
  }
  return(vstats)
}

#' Function that initiates report element with mutational signatures information
#'
#' @return rep Report structure initialized for signature data
#' @export
init_m_signature_content <- function() {

  rep <- list()
  rep[["eval"]] <- FALSE

  rep[["variant_set"]] <- list()
  rep[["variant_set"]][["all"]] <- data.frame()
  rep[["missing_data"]] <- FALSE
  rep[["result"]] <- list()
  rep[["result"]][["vr"]] <- NULL
  rep[["result"]][["mut_mat"]] <- NULL
  rep[["result"]][["chromosomes"]] <- NULL
  rep[["result"]][["no_site_prevalence"]] <- FALSE
  rep[["result"]][["tsv"]] <- data.frame()
  rep[["result"]][["contributions"]] <- list()
  rep[["result"]][["contributions"]][["per_signature"]] <- data.frame()
  rep[["result"]][["contributions"]][["per_group"]] <- data.frame()
  rep[["result"]][["scale_fill_values"]] <- c()
  rep[["result"]][["scale_fill_names"]] <- c()
  rep[["result"]][["reference_data"]] <- list()
  rep[["result"]][["reference_data"]][["aetiology"]] <- data.frame()
  rep[["result"]][["reference_data"]][["matrix"]] <- NULL
  rep[["result"]][["goodness_of_fit"]] <- NA

  return(rep)
}

#' Function that initiates report element with MSI classification
#'
#' @export
init_msi_content <- function() {
  rep <- list()

  rep[["eval"]] <- FALSE
  rep[["missing_data"]] <- FALSE
  rep[["prediction"]] <- list()

  return(rep)

}

#' Function that initiates report element with kataegis information
#'
#' @export
init_kataegis_content <- function() {
  rep <- list()

  rep[["eval"]] <- FALSE
  rep[["events"]] <- data.frame()

  return(rep)

}

#' Function that initiates report element with expression information
#'
#' @return rep Report structure initialized for expression data
#' @export
init_expression_content <- function() {

  rep <- list()
  rep[["eval"]] <- FALSE
  rep[['similarity_analysis']] <- list()
  for(source in c('tcga','depmap','treehouse')){
    rep[['similarity_analysis']][[source]] <- data.frame()
  }
  rep[['expression']] <- data.frame()
  rep[['csq_expression']] <- data.frame()
  rep[['outliers']] <- data.frame()
  for(cat in c('immune_contexture','drug_targets')){
    rep[[cat]] <- data.frame()
  }

  return(rep)

}

#' Function that initiates report element with rainfall information
#'
#' @return rep Report structure initialized for rainfall data
#' @export
init_rainfall_content <- function() {

  rep <- list()

  rep[["eval"]] <- FALSE
  rep[["rfdata"]] <- list()
  for (e in c("intercept", "chr_cum", "cex")) {
    rep[["rfdata"]][[e]] <-
      numeric()
  }
  for (e in c("ylim", "cex_text")) {
    rep[["rfdata"]][[e]] <- integer()
  }
  rep[["rfdata"]][["labels"]] <-
    character()
  rep[["rfdata"]][["colors"]] <-
    character()
  rep[["rfdata"]][["data"]] <-
    data.frame()

  return(rep)

}

#' Function that initiates report element with tumor-only information
#'
#' @return rep Report structure initialized for tumor-only stats
#' @export
init_tumor_only_content <- function() {

  rep <- list()
  rep[["eval"]] <- FALSE
  rep[["upset_data"]] <- data.frame()
  rep[["upset_plot_valid"]] <- FALSE
  rep[["vfilter"]] <- list()

  for (successive_filter in
       c("unfiltered_n",
         "gnomad_n_remain",
         "clinvar_n_remain",
         "pon_n_remain",
         "hom_n_remain",
         "het_n_remain",
         "dbsnp_n_remain",
         "nonexonic_n_remain",
         "gnomad_frac_remain",
         "clinvar_frac_remain",
         "dbsnp_frac_remain",
         "pon_frac_remain",
         "hom_frac_remain",
         "het_frac_remain",
         "nonexonic_frac_remain")) {
    rep[["vfilter"]][[successive_filter]] <- 0
  }
  return(rep)
}

#' Function that initiates report element with variant data
#'
#' @return rep Report structure initialized for variant data
#' @export
init_var_content <- function() {

  rep <- list()

  rep[["eval"]] <- FALSE
  rep[["clin_eitem"]] <- list()
  rep[["disp"]] <- list()
  rep[["variant_set"]] <- list()
  rep[["v_stat"]] <- list()
  rep[["zero"]] <- FALSE
  for (tumorclass in c("any_ttype", "other_ttype", "query_ttype")) {
    rep[["clin_eitem"]][[tumorclass]] <- list()
    for (e_type in c("prognostic", "diagnostic", "predictive")) {
      for (e_level in c("A_B", "C_D_E", "any")) {
        rep[["clin_eitem"]][[tumorclass]][[e_type]][[e_level]] <-
          data.frame()
      }
    }
  }
  return(rep)
}

#' Function that initiates report element with germline variant information (CPSR)
#'
#' @return rep Report structure initialized for germline data (CPSR)
#' @export
init_germline_content <- function() {
  rep <- list()

  rep[['callset']] <- list()
  rep[["max_dt_rows"]] <- 0
  rep[["eval"]] <- FALSE
  rep[['callset']][["variant"]] <- list()
  rep[['callset']][["variant_display"]] <- list()
  rep[['callset']][['biomarker_evidence']] <- list()
  rep[["zero"]] <- FALSE
  for (t in c("all",
              "cpg_non_sf",
              "gwas",
              "bm",
              "sf")) {
    rep[["callset"]][["variant"]][[t]] <-
      data.frame()
    rep[["callset"]][["variant_display"]][[t]] <-
      data.frame()
  }

  for (cl in c("v_stat",
               "v_stat_cpg",
               "v_stat_sf")) {
    rep[[cl]] <- list()
    for (t in c("n",
                "n_snv",
                "n_indel",
                "n_coding",
                "n_noncoding",
                "n_p",
                "n_lp",
                "n_vus",
                "n_lb",
                "n_b")) {
      rep[[cl]][[t]] <- 0
    }
  }

  rep[['v_stat_bm']] <- list()
  for (cl in c("n_var_eitems",
               "n_eitems_predictive",
               "n_eitems_prognostic",
               "n_eitems_diagnostic",
               "n_eitems_predisposing")) {
    rep[['v_stat_bm']][[cl]] <- 0
  }


  return(rep)

}

#' Function that loads YAML data with settings and file paths
#' to annotated molecular profiles
#'
#' @param yml_fname Yaml file with configuration settings and paths
#' to annotated molecular datasets
#' @param report_mode Type of report ('PCGR' or 'CPSR')
#'
#' @export
load_yaml <- function(yml_fname, report_mode = "CPSR") {

  if (!file.exists(yml_fname)) {
    log4r_fatal(
      paste0("YAML file '",yml_fname,"' does not exist - exiting"))
  }
  report_settings <- yaml::read_yaml(yml_fname)
  missing_yaml_info <- F
  for(t in c('sample_id',
             'genome_assembly',
             'workflow',
             'output_dir')) {
    if (is.null(report_settings[[t]])) {
      missing_yaml_info <- T
    }else{
      if (identical(typeof(report_settings[[t]]),"character") == F) {
        missing_yaml_info <- T
      }
    }
  }
  for(t in c('conf',
             'molecular_data',
             'reference_data',
             'software')) {
    if (is.null(report_settings[[t]])) {
      missing_yaml_info <- T
    }else{
      if (identical(typeof(report_settings[[t]]),"list") == F) {
        missing_yaml_info <- T
      }
    }
  }
  if (missing_yaml_info == F) {
    log4r_info(paste0(
      "Successfully parsed YAML configuration file - reporting mode: ", report_mode))
  }else{
    log4r_fatal(paste0(
      "YAML configuration lacks necessary properties - reporting mode: ", report_mode))
  }

  ## check the validity of the yaml file
  ## check that it matches the report_mode
  ## return it

  if (report_settings[['workflow']] != report_mode) {
    log4r_fatal(
      paste0("Cannot read YAML file from ",
             report_settings[['workflow']],
             "when report_mode is set to '",
             report_mode,"'"))
  }

  ref_data <- list()
  if (dir.exists(
    report_settings[['reference_data']][['path']])) {
    ref_data <- load_reference_data(
      pcgr_db_assembly_dir = report_settings[['reference_data']][['path']],
      genome_assembly = report_settings[['genome_assembly']]
    )
  }else{
    log4r_fatal(
      paste0("Reference data directory ",
             report_settings[['reference_data']][['path']],
             " does not exist - exiting"))

  }

  if (identical(
    typeof(report_settings[['conf']][['sample_properties']][['phenotype']]),
               "character")) {
    report_settings[['conf']][['sample_properties']][['phenotype']] <-
      ref_data[['phenotype']][['oncotree']]

    if (report_mode == "PCGR") {
      report_settings[['conf']][['sample_properties']][['phenotype']] <-
        report_settings[['conf']][['sample_properties']][['phenotype']] |>
        dplyr::filter(!is.na(.data$PRIMARY_SITE))

    }

  }else{
    if (identical(
      typeof(
        report_settings[['conf']][['sample_properties']][['phenotype']]),
      "list")) {
      report_settings[['conf']][['sample_properties']][['phenotype']] <-
        as.data.frame(
          rrapply::rrapply(
            report_settings[['conf']][['sample_properties']][['phenotype']],
            how = "bind"))
    }

    for(col in c('do_id','do_name','efo_id','efo_name',
                 'icd10_code','ot_name','ot_primary_site',
                 'primary_site','ot_code','ot_code_path')) {

      if (NROW(report_settings[['conf']][['sample_properties']][['phenotype']][
        report_settings[['conf']][['sample_properties']][['phenotype']][[col]] == "NaN",]) > 0) {
          report_settings[['conf']][['sample_properties']][['phenotype']][
            report_settings[['conf']][['sample_properties']][['phenotype']][[col]] == "NaN",col] <-
          as.character(NA)
      }

      if (NROW(report_settings[['conf']][['sample_properties']][['phenotype']][
        is.nan(report_settings[['conf']][['sample_properties']][['phenotype']][[col]]),]) > 0) {
        report_settings[['conf']][['sample_properties']][['phenotype']][
          is.nan(report_settings[['conf']][['sample_properties']][['phenotype']]),col] <-
          as.character(NA)
      }
    }

    for(col in c('do_cancer_slim','ot_level')) {
      if (NROW(report_settings[['conf']][['sample_properties']][['phenotype']][
        is.nan(report_settings[['conf']][['sample_properties']][['phenotype']][[col]]),]) > 0) {
        report_settings[['conf']][['sample_properties']][['phenotype']][
          is.nan(report_settings[['conf']][['sample_properties']][['phenotype']][[col]]),col] <-
          as.numeric(NA)
      }
    }
  }

  report_settings[['reference_data']][['source_metadata']] <-
    as.data.frame(
      rrapply::rrapply(
        report_settings$reference_data$source_metadata,
        how = "bind"))

  for(col in c('source_version',
               'source_license',
               'source_license_url',
               'source_url',
               'source_citation')) {

    if (NROW(report_settings[['reference_data']][['source_metadata']][
      report_settings[['reference_data']][['source_metadata']][[col]] == "NaN",]) > 0) {
      report_settings[['reference_data']][['source_metadata']][
        report_settings[['reference_data']][['source_metadata']][[col]] == "NaN",col] <-
        as.character(NA)
    }
  }

  ## temporary ACMG url fix
  for(i in 1:NROW(report_settings[['reference_data']][['source_metadata']])){
    if(report_settings[['reference_data']][['source_metadata']][i,"source_abbreviation"] == "acmg_sf"){
      report_settings[['reference_data']][['source_metadata']][i,"source_url"] <-
        "https://pubmed.ncbi.nlm.nih.gov/37347242/"
    }
  }

  if (report_mode == "CPSR") {
    report_settings[['conf']][['gene_panel']][['panel_genes']] <-
      as.data.frame(
        rrapply::rrapply(
          report_settings$conf$gene_panel$panel_genes,
          how = "bind"))

    if (NROW(report_settings[['conf']][['gene_panel']][['panel_genes']]) == 1) {
      for(e in c('panel_id','panel_url','panel_version')) {
        report_settings[['conf']][['gene_panel']][['panel_genes']][,e] <- NA
      }
      for(e in c('mod','moi')) {
        if (is.nan(report_settings$conf$gene_panel$panel_genes[,e])) {
          report_settings[['conf']][['gene_panel']][['panel_genes']][,e] <- NA
        }
      }
    }else{

      for(col in c('panel_id','panel_version')) {

        if (NROW(report_settings[['conf']][['gene_panel']][['panel_genes']][
          is.nan(report_settings[['conf']][['gene_panel']][['panel_genes']][[col]]),]) > 0) {
          report_settings[['conf']][['gene_panel']][['panel_genes']][
            is.nan(report_settings[['conf']][['gene_panel']][['panel_genes']][[col]]),col] <-
            as.numeric(NA)
        }
      }
      for(col in c('mod','moi')) {

        if (NROW(report_settings[['conf']][['gene_panel']][['panel_genes']][
          is.nan(report_settings[['conf']][['gene_panel']][['panel_genes']][[col]]),]) > 0) {
          report_settings[['conf']][['gene_panel']][['panel_genes']][
            is.nan(report_settings[['conf']][['gene_panel']][['panel_genes']][[col]]),col] <-
            as.character(NA)
        }
      }
    }
    colnames(
      report_settings[['conf']][['gene_panel']][['panel_genes']]) <-
      toupper(
        colnames(
          report_settings[['conf']][['gene_panel']][['panel_genes']])
      )
  }

  report_settings$conf$report_color <-
    pcgrr::color_palette[["report_color"]][["values"]][1]

  if(!is.null(ref_data$assembly$chrom_coordinates)){
    report_settings$chrom_coordinates <-
      ref_data$assembly$chrom_coordinates
  }

  if (report_mode == "PCGR" &
    !is.null(report_settings$conf$assay_properties)) {
   if (report_settings$conf$assay_properties$vcf_tumor_only == 1) {
     report_settings$conf$report_color <-
       pcgrr::color_palette[["report_color"]][["values"]][2]
   }
  }

  return(list('settings' = report_settings,
              'ref_data' = ref_data))

}
