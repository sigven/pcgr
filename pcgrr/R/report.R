
#' Function that initiates PCGR/CPSR report object
#'
#' @param config Object with configuration parameters
#' @param class report analysis section (NULL defaults to full report)
#' @param pcgr_data PCGR data bundle
#' @param type Type of report ('somatic' or 'germline')
#' @param virtual_panel_id identifier(s) for virtual panel id
#' @param custom_bed custom BED file with target loci for screening
#' @export
init_report <- function(config = NULL,
                        class = NULL,
                        pcgr_data = NULL,
                        type = "somatic",
                        virtual_panel_id = "-1",
                        custom_bed = NULL) {

  # invisible(assertthat::assert_that(
  #   !is.null(type), msg = "Argument 'type' must be of type character"))
  invisible(assertthat::assert_that(
    is.character(type), msg = "Argument 'type' must be of type character"))
  invisible(assertthat::assert_that(
    type == "germline" | type == "somatic", msg = "Argument 'type' must be of either 'germline' or 'somatic'"))

  if(!is.null(class)){
    invisible(assertthat::assert_that(
      is.character(class), msg = "Argument 'class' must be of type character"))
  }

  report <- list()
  for (elem in c("metadata", "content")) {
    report[[elem]] <- list()
  }

  if(!is.null(config) & !is.null(pcgr_data)){
    report_metadata <- pcgrr::set_report_metadata(
      config,
      pcgr_data,
      report_type = type,
      virtual_panel_id = virtual_panel_id,
      custom_bed = custom_bed
    )

    if (!is.null(report_metadata)) {
      report[["metadata"]] <- report_metadata
    }
  }

  if (type == "germline") {
    if (!is.null(pcgr_data)) {
      report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]] <-
        dplyr::filter(pcgr_data[["phenotype_ontology"]][["oncotree"]],
                      is.na(.data$primary_site))
    }

    report[["content"]][["snv_indel"]] <-
      init_germline_content()

    if (!is.null(report[["metadata"]][["config"]][["popgen"]])) {
      if (report[["metadata"]][["config"]][["popgen"]][["pop_gnomad"]] != "") {
        pop_tag_info <-
          pcgrr::get_population_tag(config[["popgen"]][["pop_gnomad"]],
                                    db = "GNOMAD", subset = "non_cancer")
        report[["metadata"]][["config"]][["popgen"]][["vcftag_gnomad"]] <-
          pop_tag_info$vcf_tag
        report[["metadata"]][["config"]][["popgen"]][["popdesc_gnomad"]] <-
          pop_tag_info$pop_description
      }
    }

  }else{

    if (!is.null(pcgr_data)) {
      if (config[["t_props"]][["tumor_type"]] != "Cancer, NOS") {
        tumor_group_entry <-
          dplyr::filter(
            pcgr_data[["phenotype_ontology"]][["cancer_groups"]],
            .data$primary_site == config[["t_props"]][["tumor_type"]])
        if (nrow(tumor_group_entry) == 1) {
          report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]] <-
            dplyr::filter(
              pcgr_data[["phenotype_ontology"]][["oncotree"]],
              .data$primary_site == config[["t_props"]][["tumor_type"]])
        }else{
          log4r_info(paste0("Cannot find tumor type ",
                                   config[["t_props"]][["tumor_type"]],
                                   " in list of primary sites"))
          report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]] <-
            pcgr_data[["phenotype_ontology"]][["oncotree"]]
        }
      }else{
        report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]] <-
          pcgr_data[["phenotype_ontology"]][["oncotree"]]
      }
    }
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
      if (a_elem == "tmb" & !is.null(pcgr_data)) {
        report[["content"]][[a_elem]] <-
          init_tmb_content(tcga_tmb = pcgr_data$tcga$tmb,
                           config = config)
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
  if (!is.null(class)) {
    if (!is.null(report[["content"]][[class]])) {
      return(report[["content"]][[class]])
    }
  }

  return(report)
}

#' Function that initates PCGR report object
#'
#' @param report PCGR final report
#' @param report_data Object with PCGR report data
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

#' Function that set PCGR report metadata
#'
#' @param config PCGR config object
#' @param pcgr_data PCGR data list object
#' @param virtual_panel_id identifier for virtual panel
#' @param report_type type of report ('germline' or 'somatic')
#' @param custom_bed  custom BED file with target loci for screening
#'
#' @return report_metadata Object with PCGR report metadata
#'
#' @export
set_report_metadata <- function(config,
                                pcgr_data,
                                virtual_panel_id = "-1",
                                report_type = NULL,
                                custom_bed = NULL) {

  if(is.null(pcgr_data)) {
    return(NULL)
  }

  report_metadata <- list()
  report_metadata[["pcgr_db_release"]] <- pcgr_data[["release_notes"]]
  report_metadata[["pcgr_version"]] <- config[['required_args']][['pcgr_version']]
  report_metadata[["cpsr_version"]] <- config[['required_args']][['cpsr_version']]
  report_metadata[["genome_assembly"]] <-
    pcgr_data[["assembly"]][["grch_name"]]
  report_metadata[["sample_name"]] <- config[['required_args']][['sample_name']]
  report_metadata[["report_type"]] <- report_type
  report_metadata[["config"]] <- config
  report_metadata[["color_palette"]] <- pcgrr::color_palette
  report_metadata[["color_none"]] <-
    pcgrr::color_palette[["none"]][["values"]][1]
  report_metadata[["color_value_box"]] <-
    pcgrr::color_palette[["report_color"]][["values"]][1]
  if(report_type == "somatic" &
     !is.null(report_metadata[["config"]][["assay_props"]])) {
    if(report_metadata[["config"]][["assay_props"]][["vcf_tumor_only"]] == T) {
      report_metadata[["color_value_box"]] <-
        pcgrr::color_palette[["report_color"]][["values"]][2]
    }
  }
  report_metadata[["phenotype_ontology"]] <- list()
  report_metadata[["phenotype_ontology"]][["oncotree"]] <-
    pcgr_data[["phenotype_ontology"]][["oncotree"]]
  report_metadata[["phenotype_ontology"]][["oncotree_query"]] <- NULL


  if (virtual_panel_id != "-1") {
    report_metadata[["gene_panel"]] <- list()
    report_metadata[["gene_panel"]][["genes"]] <- data.frame()

    panel_ids <- stringr::str_split(virtual_panel_id, ",")[[1]]
    for(pid in panel_ids){
      p <- as.integer(pid)
      panel_genes <- pcgr_data[["virtual_gene_panels"]] %>%
        dplyr::filter(.data$id == p) %>%
        dplyr::select(.data$symbol,
                      .data$ensembl_gene_id,
                      .data$confidence_level,
                      .data$panel_name,
                      .data$panel_id,
                      .data$id,
                      .data$panel_version,
                      .data$panel_url,
                      .data$entrezgene,
                      .data$genename)

      report_metadata[["gene_panel"]][["genes"]] <-
        report_metadata[["gene_panel"]][["genes"]] %>%
        dplyr::bind_rows(panel_genes) %>%
        dplyr::arrange(.data$symbol)
    }

    if (config[["diagnostic_grade_only"]] == T) {
      report_metadata[["gene_panel"]][["genes"]] <-
        report_metadata[["gene_panel"]][["genes"]] %>%
        dplyr::filter(.data$confidence_level == 3 | .data$confidence_level == 4)
    }


    ## Single gene panel from Genomics England
    if(!(stringr::str_detect(virtual_panel_id,","))){

      report_metadata[["gene_panel"]][["name"]] <-
        paste0(unique(report_metadata[["gene_panel"]][["genes"]]$panel_name),
               " - v",
               unique(report_metadata[["gene_panel"]][["genes"]]$panel_version))
      report_metadata[["gene_panel"]][["url"]] <-
        unique(report_metadata[["gene_panel"]][["genes"]]$panel_url)
      report_metadata[["gene_panel"]][["confidence"]] <-
        "Diagnostic-grade/Borderline/Low evidence (GREEN/AMBER/RED)"
      if (config[["diagnostic_grade_only"]] == T & virtual_panel_id != "0") {
        report_metadata[["gene_panel"]][["confidence"]] <-
          "Diagnostic-grade only (GREEN)"
      }
      if (virtual_panel_id == "0") {
        report_metadata[["gene_panel"]][["confidence"]] <-
          "Exploratory virtual gene panel (research)"
      }
    }else{

      panels_included <-
        report_metadata$gene_panel$genes %>%
        dplyr::select(.data$id, .data$panel_url, .data$panel_name) %>%
        dplyr::arrange(.data$id) %>%
        dplyr::distinct() %>%
        dplyr::mutate(link = paste0("<li><a href='",.data$panel_url,
                                    "' target='_blank'>",
                                    .data$panel_name,"</a></li>"))

      report_metadata[["gene_panel"]][["genes"]] <-
        report_metadata[["gene_panel"]][["genes"]] %>%
        dplyr::select(-c(.data$panel_id, .data$panel_name, .data$confidence_level,
                         .data$panel_version, .data$panel_url,
                         .data$id)) %>%
        dplyr::mutate(panel_version = NA, confidence_level = 5,
                      panel_id = NA, panel_url = NA) %>%
        dplyr::distinct()

      report_metadata[["gene_panel"]][["url"]] <- "https://panelapp.genomicsengland.co.uk"
      report_metadata[["gene_panel"]][["name"]] <-
        "Genomics England PanelApp - custom combination"
      report_metadata[["gene_panel"]][["confidence"]] <-
        paste0("<b>N = ", NROW(report_metadata[["gene_panel"]][["genes"]]),
               "</b> unique genes combined from the following virtual panels:<br> <ul>",
               paste(panels_included$link, collapse="\n"),"\n</ul><br>\n")

    }
  }else{

    ## virtual_panel_id is "-1", indicating custom gene list
    if (!is.null(custom_bed)) {
      target_genes <-
        pcgrr::custom_bed_genes(custom_bed, pcgr_data = pcgr_data)
      if (nrow(target_genes) > 0) {

        target_genes <- target_genes %>%
          dplyr::mutate(
            confidence_level = -1,
            panel_name =
              report_metadata[["config"]][["custom_panel"]][["name"]],
            panel_id = NA, panel_version = NA,
            panel_url = NA)

        report_metadata[["gene_panel"]] <- list()
        report_metadata[["gene_panel"]][["genes"]] <- target_genes
        report_metadata[["gene_panel"]][["name"]] <-
          report_metadata[["config"]][["custom_panel"]][["name"]]
        report_metadata[["gene_panel"]][["confidence"]] <-
          "User-defined panel (custom geneset from panel 0)"
      }else{
        log4r_warn(
          paste0("Custom BED file (", custom_bed,
                 ") does not contain regions that overlap ",
                 "protein-coding transcripts (regulatory/exonic)"))
        log4r_info("Quitting report generation")
      }
    }
  }
  return(report_metadata)
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
  invisible(assertthat::assert_that(!is.null(config[['assay_props']])))
  invisible(assertthat::assert_that("target_size_mb" %in%
                                      names(config[['assay_props']])))

  rep <- list()
  rep[["eval"]] <- FALSE

  rep[["algorithm"]] <-"all_coding"
  rep[["v_stat"]] <- list()
  rep[["v_stat"]][["n_tmb"]] <- 0
  rep[["v_stat"]][["tmb_estimate"]] <- 0
  rep[["v_stat"]][["target_size_mb"]] <-
    config[["assay_props"]][["target_size_mb"]]
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
    rep[["v_stat"]][[t]] <-
      data.frame()
  }
  for (t in c("segment", "oncogene_gain", "tsgene_loss",
              "other_target", "biomarker", "tier1", "tier2")) {
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
      for (c in c("proto_oncogene", "tumor_suppressor")) {
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

  for (successive_filter in c("unfiltered_n", "onekg_n_remain",
                              "gnomad_n_remain", "clinvar_n_remain",
                              "pon_n_remain", "hom_n_remain",
                              "het_n_remain", "dbsnp_n_remain",
                              "nonexonic_n_remain", "onekg_frac_remain",
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
  for (evidence_type in c("prognostic", "diagnostic", "
                            predictive", "predisposing")) {
    rep[["clin_eitem"]][[evidence_type]] <-
      data.frame()
  }

  for (cl in c("v_stat", "v_stat_cpg",
               "v_stat_secondary")) {
    rep[[cl]] <- list()
    for (t in c("n", "n_snv", "n_indel", "n_coding", "n_noncoding")) {
      rep[[cl]][[t]] <- 0
    }
  }

  return(rep)

}
