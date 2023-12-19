
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
    yaml_data <- load_yaml(
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
        if(population != "global"){
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
        if(NROW(pop_desc_df) == 1){
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
                     "tmb",
                     "msi",
                     "cna",
                     "cna_plot",
                     "m_signature_mp",
                     "sequencing_mode",
                     "tumor_only",
                     "value_box",
                     "rainfall",
                     "kataegis",
                     "tumor_purity",
                     "tumor_ploidy",
                     "cpsr",
                     "report_display_config",
                     "clinicaltrials")) {
      report[["content"]][[a_elem]] <- list()
      report[["content"]][[a_elem]][["eval"]] <- FALSE

      if(a_elem == "tumor_purity" | a_elem == "tumor_ploidy"){
        report[["content"]][[a_elem]][["eval"]] <- TRUE
      }

      if (a_elem == "kataegis") {
        report[["content"]][[a_elem]][["events"]] <- data.frame()
      }

      if (a_elem == "rainfall") {
        report[["content"]][[a_elem]] <-
          init_rainfall_content()
      }

      if (a_elem == "cna_plot") {
        report[["content"]][[a_elem]][["png"]] <- NULL
      }

      if (a_elem == "report_display_config") {
        report[["content"]][[a_elem]] <-
          init_report_display_content()
      }

      if (a_elem == "value_box") {
        report[["content"]][[a_elem]] <-
          init_valuebox_content()
      }

      if (a_elem == "snv_indel" | a_elem == "cna") {
        report[["content"]][[a_elem]] <- init_var_content()

        if (a_elem == "snv_indel") {
          report[["content"]][[a_elem]] <-
            init_snv_indel_content(rep = report[["content"]][[a_elem]])
        }
        if (a_elem == "cna") {
          report[["content"]][[a_elem]] <-
            init_cna_content(rep = report[["content"]][[a_elem]])
        }
      }
      if (a_elem == "clinicaltrials") {
        report[["content"]][[a_elem]][["trials"]] <- data.frame()
        report[["content"]][[a_elem]][["missing_data"]] <- F
      }

      if (a_elem == "m_signature_mp") {
        report[["content"]][[a_elem]] <-
          init_m_signature_content()
      }
      if (a_elem == "tmb" & !is.null(report$ref_data$tmb)) {
        report[["content"]][[a_elem]] <-
          init_tmb_content(tcga_tmb = report$ref_data$tmb,
                           config = report$settings$conf)
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
  # if (!is.null(class)) {
  #   if (!is.null(report[["content"]][[class]])) {
  #     return(report[["content"]][[class]])
  #   }
  # }

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
#' @param tcga_tmb data frame with TMB distribution in TCGA samples
#' @param config PCGR configuration object
#
#' @return rep TMB report element
#'
#' @export
init_tmb_content <- function(tcga_tmb = NULL,
                             config = NULL){

  invisible(assertthat::assert_that(!is.null(tcga_tmb)))
  invisible(assertthat::assert_that(is.data.frame(tcga_tmb) &
                                      nrow(tcga_tmb) > 0))
  invisible(assertthat::assert_that(!is.null(config)))
  invisible(assertthat::assert_that(!is.null(config[['assay_properties']])))
  invisible(assertthat::assert_that("effective_target_size_mb" %in%
                                      names(config[['assay_properties']])))

  rep <- list()
  rep[["eval"]] <- FALSE

  rep[["algorithm"]] <-"all_coding"
  rep[["v_stat"]] <- list()
  rep[["v_stat"]][["n_tmb"]] <- 0
  rep[["v_stat"]][["tmb_estimate"]] <- 0
  rep[["v_stat"]][["effective_target_size_mb"]] <-
    config[["assay_properties"]][["effective_target_size_mb"]]
  #rep[["v_stat"]][["tmb_tertile"]] <-
    #"TMB - not determined"
  rep[["tcga_tmb"]] <- tcga_tmb

  return(rep)
}

#' Function that initiates report element with CNA information
#'
#' @param rep PCGR report structure
#
#' @return rep updated PCGR report structure - initialized for CNA content
#' @export
init_cna_content <- function(rep = NULL){

  invisible(assertthat::assert_that(!is.null(rep)))
  invisible(assertthat::assert_that(!is.null(rep[['disp']])))
  invisible(assertthat::assert_that(!is.null(rep[['variant_set']])))

  rep[["variant_set"]][["tsv"]] <-
    data.frame()
  rep[["variant_set"]][["tier1"]] <-
    data.frame()
  rep[["variant_set"]][["tier2"]] <-
    data.frame()
  for (t in c("n_cna_loss", "n_cna_gain")) {
    rep[["v_stat"]][[t]] <- 0
  }
  for (t in c("segment",
              "oncogene_gain",
              "tsgene_loss",
              "other_target",
              "biomarker",
              "tier1",
              "tier2")) {
    rep[["disp"]][[t]] <-
      data.frame()
  }
  return(rep)


}

#' Function that initiates report element with SNV/InDel information
#'
#' @param rep PCGR report structure
#
#' @return rep updated PCGR report structure - initialized for SNV/InDel content
#' @export
init_snv_indel_content <- function(rep = NULL){

  invisible(assertthat::assert_that(!is.null(rep)))
  invisible(assertthat::assert_that(!is.null(rep[['disp']])))
  invisible(assertthat::assert_that(!is.null(rep[['variant_set']])))

  for (t in c("tier1", "tier2", "tier3",
              "tier4", "noncoding")) {
    rep[["disp"]][[t]] <-
      data.frame()
    if (t == "tier3") {
      rep[["disp"]][[t]] <- list()
      for (c in c("proto_oncogene",
                  "tumor_suppressor")) {
        rep[["disp"]][[t]][[c]] <-
          data.frame()
      }
    }
  }
  for (t in c("tier1", "tier2", "tier3", "tier4", "noncoding",
              "tsv", "tsv_unfiltered", "maf", "coding", "all")) {
    rep[["variant_set"]][[t]] <-
      data.frame()
  }
  for (t in c("n", "n_snv", "n_indel", "n_coding", "n_noncoding",
              "n_tier1", "n_tier2", "n_tier3", "n_tier4")) {
    rep[["v_stat"]][[t]] <- 0
  }
  return(rep)
}

#' Function that initiates report element with mutational signatures information
#'
#' @return rep Report structure initialized for signature data
#' @export
init_m_signature_content <- function(){

  rep <- list()
  rep[["eval"]] <- FALSE

  rep[["variant_set"]] <- list()
  rep[["variant_set"]][["all"]] <- data.frame()
  rep[["missing_data"]] <- FALSE
  rep[["result"]] <- list()
  rep[["result"]][["vr"]] <- NULL
  rep[["result"]][["mut_mat"]] <- NULL
  rep[["result"]][["chromosomes"]] <- NULL
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

#init_msi_content <- function(){}
#init_kataegis_content <- function(){}

#' Function that initiates report element with rainfall information
#'
#' @return rep Report structure initialized for rainfall data
#' @export
init_rainfall_content <- function(){

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
#' @return rep Report structure initialized for tumor-only data
#' @export
init_tumor_only_content <- function(){

  rep <- list()
  rep[["eval"]] <- FALSE

  rep[["variant_set"]] <- list()
  rep[["variant_set"]][["tsv_unfiltered"]] <-
    data.frame()
  rep[["variant_set"]][["filtered"]] <-
    data.frame()
  rep[["upset_data"]] <- data.frame()
  rep[["upset_plot_valid"]] <- FALSE
  rep[["v_stat"]] <- list()

  for (successive_filter in c("unfiltered_n",
                              "gnomad_n_remain", "clinvar_n_remain",
                              "pon_n_remain", "hom_n_remain",
                              "het_n_remain", "dbsnp_n_remain",
                              "nonexonic_n_remain",
                              "gnomad_frac_remain", "clinvar_frac_remain",
                              "dbsnp_frac_remain", "pon_frac_remain",
                              "hom_frac_remain", "het_frac_remain",
                              "nonexonic_frac_remain")) {
    rep[["v_stat"]][[successive_filter]] <- 0
  }
  return(rep)
}

#' Function that initiates report element with value box information
#'
#' @return rep Report structure initialized for value box data
#' @export
init_valuebox_content <- function(){
  rep <- list()

  rep[["eval"]] <- FALSE

  rep[["tmb"]] <-
    "Not determined"
  rep[["msi"]] <-
    "Not determined"
  rep[["scna"]] <-
    "Not determined"
  rep[["tier1"]] <-
    "Not determined"
  rep[["tier2"]] <-
    "Not determined"
  rep[["signatures"]] <-
    "Not determined"
  rep[["tumor_ploidy"]] <-
    "Not provided"
  rep[["tumor_purity"]] <-
    "Not provided"
  rep[["kataegis"]] <-
    "Not determined"

  return(rep)

}

#' Function that initiates ranked report display information
#'
#' @return rep Report structure initialized for ranked display
#' @export
init_report_display_content <- function(){

  rep <- list()
  rep[["eval"]] <- FALSE

  rep[["opentargets_rank"]] <- list()
  rep[["opentargets_rank"]][["breaks"]] <-
    c(0.40, 0.55, 0.70, 0.85)
  rep[["opentargets_rank"]][["colors"]] <-
    c("#b8b8ba", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C")

  return(rep)

}


#' Function that initiates report element with variant data
#'
#' @return rep Report structure initialized for variant data
#' @export
init_var_content <- function(){

  rep <- list()

  rep[["eval"]] <- FALSE
  rep[["clin_eitem"]] <- list()
  rep[["disp"]] <- list()
  rep[["variant_set"]] <- list()
  rep[["v_stat"]] <- list()
  rep[["zero"]] <- FALSE
  for (tumorclass in c("any_ttype", "other_ttype", "specific_ttype")) {
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
init_germline_content <- function(){
  rep <- list()

  rep[["max_dt_rows"]] <- 0
  rep[["eval"]] <- FALSE
  rep[["disp"]] <- list()
  rep[["variant_set"]] <- list()
  rep[["zero"]] <- FALSE
  for (t in c("class1", "class2", "class3",
              "class4", "class5", "gwas",
              "secondary")) {
    rep[["disp"]][[t]] <-
      data.frame()
    rep[["variant_set"]][[t]] <-
      data.frame()
  }
  rep[["variant_set"]][["tsv"]] <-
    data.frame()
  rep[["clin_eitem"]] <- list()
  for (evidence_type in pcgrr::evidence_types) {
    rep[["clin_eitem"]][[evidence_type]] <- list()
    for(level in pcgrr::evidence_levels){
      rep[["clin_eitem"]][[evidence_type]][[level]] <-
        data.frame()
    }
    rep[['clin_eitem']][['all']] <- list()
    for(level in pcgrr::evidence_levels){
      rep[["clin_eitem"]][['all']][[level]] <-
        data.frame()
    }
  }

  for (cl in c("v_stat",
               "v_stat_cpg",
               "v_stat_secondary")) {
    rep[[cl]] <- list()
    for (t in c("n", "n_snv", "n_indel", "n_coding", "n_noncoding")) {
      rep[[cl]][[t]] <- 0
    }
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

  if(!file.exists(yml_fname)){
    log4r_fatal(
      paste0("YAML file '",yml_fname,"' does not exist - exiting"))
  }
  report_settings <- yaml::read_yaml(yml_fname)
  missing_yaml_info <- F
  for(t in c('sample_id',
             'genome_assembly',
             'workflow',
             'output_dir')){
    if(is.null(report_settings[[t]])){
      missing_yaml_info <- T
    }else{
      if(identical(typeof(report_settings[[t]]),"character") == F){
        missing_yaml_info <- T
      }
    }
  }
  for(t in c('conf',
             'molecular_data',
             'reference_data',
             'software')){
    if(is.null(report_settings[[t]])){
      missing_yaml_info <- T
    }else{
      if(identical(typeof(report_settings[[t]]),"list") == F){
        missing_yaml_info <- T
      }
    }
  }
  if(missing_yaml_info == F){
    log4r_info(paste0(
      "Successfully parsed YAML configuration file - reporting mode: ", report_mode))
  }else{
    log4r_fatal(paste0(
      "YAML configuration lacks necessary properties - reporting mode: ", report_mode))
  }

  ## check the validity of the yaml file
  ## check that it matches the report_mode
  ## return it

  if(report_settings[['workflow']] != report_mode){
    log4r_fatal(
      paste0("Cannot read YAML file from ",
             report_settings[['workflow']],
             "when report_mode is set to '",
             report_mode,"'"))
  }

  ref_data <- list()
  if(dir.exists(
    report_settings[['reference_data']][['path']]
  )){
    ref_data <- load_reference_data(
      pcgr_db_assembly_dir = report_settings[['reference_data']][['path']],
      genome_assembly = report_settings[['genome_assembly']]
    )
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
      "list")){
      report_settings[['conf']][['sample_properties']][['phenotype']] <-
        as.data.frame(
          rrapply::rrapply(
            report_settings[['conf']][['sample_properties']][['phenotype']],
            how = "bind"))
    }

    for(col in c('do_id','do_name','efo_id','efo_name',
                 'icd10_code','ot_name','ot_primary_site',
                 'primary_site','ot_code','ot_code_path')){

      if(NROW(report_settings[['conf']][['sample_properties']][['phenotype']][
        report_settings[['conf']][['sample_properties']][['phenotype']][[col]] == "NaN",]) > 0){
          report_settings[['conf']][['sample_properties']][['phenotype']][
            report_settings[['conf']][['sample_properties']][['phenotype']][[col]] == "NaN",col] <-
          as.character(NA)
      }

      if(NROW(report_settings[['conf']][['sample_properties']][['phenotype']][
        is.nan(report_settings[['conf']][['sample_properties']][['phenotype']][[col]]),]) > 0){
        report_settings[['conf']][['sample_properties']][['phenotype']][
          is.nan(report_settings[['conf']][['sample_properties']][['phenotype']]),col] <-
          as.character(NA)
      }
    }

    for(col in c('do_cancer_slim','ot_level')){
      if(NROW(report_settings[['conf']][['sample_properties']][['phenotype']][
        is.nan(report_settings[['conf']][['sample_properties']][['phenotype']][[col]]),]) > 0){
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
               'source_citation')){

    if(NROW(report_settings[['reference_data']][['source_metadata']][
      report_settings[['reference_data']][['source_metadata']][[col]] == "NaN",]) > 0){
      report_settings[['reference_data']][['source_metadata']][
        report_settings[['reference_data']][['source_metadata']][[col]] == "NaN",col] <-
        as.character(NA)
    }
  }

  if(report_mode == "CPSR"){
    report_settings[['conf']][['gene_panel']][['panel_genes']] <-
      as.data.frame(
        rrapply::rrapply(
          report_settings$conf$gene_panel$panel_genes,
          how = "bind"))

    if(NROW(report_settings[['conf']][['gene_panel']][['panel_genes']]) == 1){
      for(e in c('panel_id','panel_url','panel_version')){
        report_settings[['conf']][['gene_panel']][['panel_genes']][,e] <- NA
      }
      for(e in c('mod','moi')){
        if(is.nan(report_settings$conf$gene_panel$panel_genes[,e])){
          report_settings[['conf']][['gene_panel']][['panel_genes']][,e] <- NA
        }
      }
    }else{

      for(col in c('panel_id','panel_version')){

        if(NROW(report_settings[['conf']][['gene_panel']][['panel_genes']][
          is.nan(report_settings[['conf']][['gene_panel']][['panel_genes']][[col]]),]) > 0){
          report_settings[['conf']][['gene_panel']][['panel_genes']][
            is.nan(report_settings[['conf']][['gene_panel']][['panel_genes']][[col]]),col] <-
            as.numeric(NA)
        }
      }
      for(col in c('mod','moi')){

        if(NROW(report_settings[['conf']][['gene_panel']][['panel_genes']][
          is.nan(report_settings[['conf']][['gene_panel']][['panel_genes']][[col]]),]) > 0){
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

  report_settings$conf$visual_reporting[["color_palette"]] <-
    pcgrr::color_palette
  report_settings$conf$visual_reporting[["color_none"]] <-
    pcgrr::color_palette[["none"]][["values"]][1]
  report_settings$conf$visual_reporting[["color_value_box"]] <-
    pcgrr::color_palette[["report_color"]][["values"]][1]
  if(report_mode == "PCGR" &
     !is.null(report_settings$conf$assay_properties)) {
    if(report_settings$conf$assay_properties$vcf_tumor_only == 1) {
      report_settings$conf$visual_reporting[["color_value_box"]] <-
        pcgrr::color_palette[["report_color"]][["values"]][2]
    }
  }

  return(list('settings' = report_settings,
              'ref_data' = ref_data))

}
