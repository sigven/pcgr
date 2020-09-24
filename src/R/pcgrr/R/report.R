
#' Function that initiates PCGR/CPSR report object
#'
#' @param config Object with configuration parameters
#' @param sample_name sample identifier
#' @param class report analysis section (NULL defaults to full report)
#' @param pcgr_data PCGR data bundle
#' @param pcgr_version PCGR software version
#' @param cpsr_version CPSR software version
#' @param type somatic or germline
#' @param virtual_panel_id identifier for virtual panel id
#' @param custom_bed custom BED file with target loci for screening
#' @param diagnostic_grade_only choose only clinical grade genes from genomics england panels

init_report <- function(config = NULL, sample_name = "SampleX",
                            class = NULL, pcgr_data = NULL, pcgr_version = "dev",
                            cpsr_version = "dev", type = "somatic",
                            virtual_panel_id = -1, custom_bed = NULL) {

  report <- list()
  for (elem in c("metadata", "content")) {
    report[[elem]] <- list()
  }

  report_metadata <- pcgrr::set_report_metadata(config, pcgr_data,
                                       cpsr_version = cpsr_version,
                                       pcgr_version = pcgr_version,
                                       report_type = type,
                                       virtual_panel_id = virtual_panel_id,
                                       custom_bed = custom_bed,
                                       sample_name = sample_name)

  if(!is.null(report_metadata)){
    report[['metadata']] <- report_metadata
  }

  if (type == "germline") {
    if (!is.null(pcgr_data)) {
      report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]] <-
        dplyr::filter(pcgr_data[["phenotype_ontology"]][["oncotree"]],
                      is.na(primary_site))
    }
    analysis_element <- "snv_indel"
    report[["content"]][[analysis_element]] <- list()
    report[["content"]][[analysis_element]][["max_dt_rows"]] <- 0
    report[["content"]][[analysis_element]][["eval"]] <- FALSE
    report[["content"]][[analysis_element]][["variant_display"]] <- list()
    report[["content"]][[analysis_element]][["variant_set"]] <- list()
    report[["content"]][[analysis_element]][["zero"]] <- FALSE
    for (t in c("class1", "class2", "class3", "class4", "class5", "gwas", "sf")) {
      report[["content"]][[analysis_element]][["variant_display"]][[t]] <- data.frame()
      report[["content"]][[analysis_element]][["variant_set"]][[t]] <- data.frame()
    }
    report[["content"]][[analysis_element]][["variant_set"]][["tsv"]] <- data.frame()
    report[["content"]][[analysis_element]][["clin_eitem"]] <- list()
    for (evidence_type in c("prognostic", "diagnostic", "predictive", "predisposing")) {
      report[["content"]][[analysis_element]][["clin_eitem"]][[evidence_type]] <- data.frame()
    }

    for (cl in c("variant_statistic", "variant_statistic_cpg", "variant_statistic_sf")) {
      report[["content"]][[analysis_element]][[cl]] <- list()
      for (t in c("n", "n_snv", "n_indel", "n_coding", "n_noncoding")) {
        report[["content"]][[analysis_element]][[cl]][[t]] <- 0
      }
    }

    if (!is.null(report[["metadata"]][["config"]][["popgen"]])) {
      if (report[["metadata"]][["config"]][["popgen"]][["pop_gnomad"]] != "") {
        pop_tag_info <- pcgrr::get_population_tag(config[["popgen"]][["pop_gnomad"]], db = "GNOMAD", subset = "non_cancer")
        report[["metadata"]][["config"]][["popgen"]][["vcftag_gnomad"]] <- pop_tag_info$vcf_tag
        report[["metadata"]][["config"]][["popgen"]][["popdesc_gnomad"]] <- pop_tag_info$pop_description
      }
    }

  }else{

    if (!is.null(pcgr_data)) {
      if (config[['tumor_properties']][['tumor_type']] != "Cancer, NOS") {
        tumor_group_entry <- dplyr::filter(pcgr_data[["phenotype_ontology"]][["cancer_groups"]],
                                           primary_site == config[['tumor_properties']][['tumor_type']])
        if (nrow(tumor_group_entry) == 1) {
          report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]] <-
            dplyr::filter(pcgr_data[["phenotype_ontology"]][["oncotree"]],
                          primary_site == config[['tumor_properties']][['tumor_type']])
        }else{
          rlogging::message(paste0("Cannot find tumor type ", config[['tumor_properties']][['tumor_type']], " in list of primary sites"))
          report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]] <-
            pcgr_data[["phenotype_ontology"]][["oncotree"]]
        }
      }else{
        report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]] <-
          pcgr_data[["phenotype_ontology"]][["oncotree"]]
      }
    }
    for (analysis_element in c("snv_indel", "tmb", "msi", "cna", "cna_plot",
                               "m_signature_mp", "sequencing_mode", "tumor_only", "value_box",
                               "rainfall", "kataegis", "tumor_purity", "tumor_ploidy",
                               "report_display_config","clinicaltrials")) {
      report[["content"]][[analysis_element]] <- list()
      report[["content"]][[analysis_element]][["eval"]] <- FALSE

      if (analysis_element == "kataegis") {
        report[["content"]][[analysis_element]][["events"]] <- data.frame()
      }

      if (analysis_element == "rainfall") {
        report[["content"]][[analysis_element]][["rfdata"]] <- list()
        for(e in c('intercept','chr_cum','cex')){
          report[['content']][[analysis_element]][['rfdata']][[e]] <- numeric()
        }
        for(e in c('ylim','cex_text')){
          report[['content']][[analysis_element]][['rfdata']][[e]] <- integer()
        }
        report[['content']][[analysis_element]][['rfdata']][['labels']] <- character()
        report[['content']][[analysis_element]][['rfdata']][['colors']] <- character()
        report[['content']][[analysis_element]][['rfdata']][['data']] <- data.frame()

      }

      if (analysis_element == "cna_plot") {
        report[["content"]][[analysis_element]][["png"]] <- NULL
      }

      if (analysis_element == "report_display_config") {
        report[["content"]][[analysis_element]][["opentargets_rank"]] <- list()
        report[["content"]][[analysis_element]][["opentargets_rank"]][["breaks"]] <- c(0.40, 0.55, 0.70, 0.85)
        report[["content"]][[analysis_element]][["opentargets_rank"]][["colors"]] <-
          c("#b8b8ba", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C")
      }
      if (analysis_element == "tumor_purity" | analysis_element == "tumor_ploidy") {
        report[["content"]][[analysis_element]][["estimate"]] <- "NA"
        if (!is.null(report[["metadata"]][["config"]][["tumor_properties"]])) {
          if (!is.null(report[["metadata"]][["config"]][["tumor_properties"]][[analysis_element]])) {
            report[["content"]][[analysis_element]][["estimate"]] <-
              report[["metadata"]][["config"]][["tumor_properties"]][[analysis_element]]
            report[["content"]][[analysis_element]][["eval"]] <- TRUE
          }
        }
      }

      if (analysis_element == "value_box") {
        report[["content"]][[analysis_element]][["tmb"]] <- "TMB:\nNot determined"
        report[["content"]][[analysis_element]][["msi"]] <- "MSI:\nNot determined"
        report[["content"]][[analysis_element]][["scna"]] <- "SCNA:\nNot determined"
        report[["content"]][[analysis_element]][["tier1"]] <- "Tier 1 variants:\nNot determined"
        report[["content"]][[analysis_element]][["tier2"]] <- "Tier 2 variants:\nNot determined"
        report[["content"]][[analysis_element]][["signatures"]] <- "Mutational signatures:\nNot determined"
        report[["content"]][[analysis_element]][["tumor_ploidy"]] <- "Tumor ploidy:\nNot provided/determined"
        report[["content"]][[analysis_element]][["tumor_purity"]] <- "Tumor purity:\nNot provided/determined"
        report[["content"]][[analysis_element]][["kataegis"]] <- "Kataegis events:\nNot determined"
      }

      if (analysis_element == "snv_indel" | analysis_element == "cna") {
        report[["content"]][[analysis_element]][["clin_eitem"]] <- list()
        report[["content"]][[analysis_element]][["variant_display"]] <- list()
        report[["content"]][[analysis_element]][["variant_set"]] <- list()
        report[["content"]][[analysis_element]][["variant_statistic"]] <- list()
        report[["content"]][[analysis_element]][["zero"]] <- FALSE
        for (tumorclass in c("any_ttype", "other_ttype", "specific_ttype")) {
          report[["content"]][[analysis_element]][["clin_eitem"]][[tumorclass]] <- list()
          for (evidence_type in c("prognostic", "diagnostic", "predictive")) {
            for (evidence_level in c("A_B", "C_D_E", "any")) {
              report[["content"]][[analysis_element]][["clin_eitem"]][[tumorclass]][[evidence_type]][[evidence_level]] <- data.frame()
            }
          }
        }
        if (analysis_element == "snv_indel") {
          for (t in c("tier1", "tier2", "tier3", "tier4", "noncoding")) {
            report[["content"]][[analysis_element]][["variant_display"]][[t]] <- data.frame()
            if (t == "tier3") {
              report[["content"]][[analysis_element]][["variant_display"]][[t]] <- list()
              for (c in c("proto_oncogene", "tumor_suppressor")) {
                report[["content"]][[analysis_element]][["variant_display"]][[t]][[c]] <- data.frame()
              }
            }
          }
          for (t in c("tier1", "tier2", "tier3", "tier4", "noncoding",
                      "tsv", "tsv_unfiltered", "maf", "coding", "all")) {
            report[["content"]][[analysis_element]][["variant_set"]][[t]] <- data.frame()
          }
          for (t in c("n", "n_snv", "n_indel", "n_coding", "n_noncoding",
                      "n_tier1", "n_tier2", "n_tier3", "n_tier4")) {
            report[["content"]][[analysis_element]][["variant_statistic"]][[t]] <- 0
          }
        }

        if (analysis_element == "cna") {
          report[["content"]][[analysis_element]][["variant_set"]][["tsv"]] <- data.frame()
          report[["content"]][[analysis_element]][["variant_set"]][['tier1']] <- data.frame()
          report[["content"]][[analysis_element]][["variant_set"]][['tier2']] <- data.frame()
          for (t in c("n_cna_loss", "n_cna_gain")) {
            report[["content"]][[analysis_element]][["variant_statistic"]][[t]] <- data.frame()
          }
          for (t in c("segment", "oncogene_gain", "tsgene_loss",
                      "other_target", "biomarker", "tier1", "tier2")) {
            report[["content"]][[analysis_element]][["variant_display"]][[t]] <- data.frame()
          }

        }
      }
      if (analysis_element == "clinicaltrials"){
        report[["content"]][[analysis_element]][['trials']] <- data.frame()
        report[["content"]][[analysis_element]][['missing_data']] <- F
      }

      # if (analysis_element == "m_signature") {
      #   report[["content"]][[analysis_element]][["variant_set"]] <- list()
      #   report[["content"]][[analysis_element]][["variant_set"]][["all"]] <- data.frame()
      #   report[["content"]][[analysis_element]][["missing_data"]] <- FALSE
      #   report[["content"]][[analysis_element]][["result"]] <- list()
      # }
      if (analysis_element == "m_signature_mp") {
        report[["content"]][[analysis_element]][["variant_set"]] <- list()
        report[["content"]][[analysis_element]][["variant_set"]][["all"]] <- data.frame()
        report[["content"]][[analysis_element]][["missing_data"]] <- FALSE
        report[["content"]][[analysis_element]][["result"]] <- list()
        report[["content"]][[analysis_element]][["result"]][["vr"]] <- NULL
        report[["content"]][[analysis_element]][["result"]][["mut_mat"]] <- NULL
        report[["content"]][[analysis_element]][["result"]][["chromosomes"]] <- NULL
        report[["content"]][[analysis_element]][["result"]][["tsv"]] <- data.frame()
        report[["content"]][[analysis_element]][["result"]][["contributions"]] <- list()
        report[["content"]][[analysis_element]][["result"]][["contributions"]][['per_signature']] <- data.frame()
        report[["content"]][[analysis_element]][["result"]][["contributions"]][['per_group']] <- data.frame()
        report[["content"]][[analysis_element]][["result"]][["scale_fill_values"]] <- c()
        report[["content"]][[analysis_element]][["result"]][["scale_fill_names"]] <- c()
        report[["content"]][[analysis_element]][["result"]][["reference_data"]] <- list()
        report[["content"]][[analysis_element]][["result"]][["reference_data"]][['aetiology']] <- data.frame()
        report[["content"]][[analysis_element]][["result"]][["reference_data"]][['matrix']] <- NULL
        report[["content"]][[analysis_element]][["result"]][["goodness_of_fit"]] <- NA
      }
      if (analysis_element == "tmb") {
        report[["content"]][[analysis_element]][["variant_statistic"]] <- list()
        report[["content"]][[analysis_element]][["variant_statistic"]][["n_tmb"]] <- 0
        report[["content"]][[analysis_element]][["variant_statistic"]][["tmb_estimate"]] <- 0
        report[["content"]][[analysis_element]][["variant_statistic"]][["target_size_mb"]] <- config[["assay_properties"]][["target_size_mb"]]
        report[["content"]][[analysis_element]][["variant_statistic"]][["tmb_tertile"]] <- "TMB - not determined"
        if (!is.null(pcgr_data)) {
          report[["content"]][[analysis_element]][["tcga_tmb"]] <- pcgr_data[["tcga"]][["tmb"]]
        }
      }
      if (analysis_element == "msi") {
        report[["content"]][[analysis_element]][["missing_data"]] <- FALSE
        report[["content"]][[analysis_element]][["prediction"]] <- list()
      }
      if (analysis_element == "tumor_only") {
        report[["content"]][[analysis_element]][["variant_set"]] <- list()
        report[["content"]][[analysis_element]][["variant_set"]][["unfiltered"]] <- data.frame()
        report[["content"]][[analysis_element]][["variant_set"]][["filtered"]] <- data.frame()
        report[["content"]][[analysis_element]][["upset_data"]] <- data.frame()
        report[["content"]][[analysis_element]][["upset_plot_valid"]] <- FALSE
        report[["content"]][[analysis_element]][["variant_statistic"]] <- list()

        for (successive_filter in c("unfiltered_n", "onekg_n_remain", "gnomad_n_remain", "clinvar_n_remain",
                                   "pon_n_remain", "hom_n_remain", "het_n_remain", "dbsnp_n_remain",
                                   "nonexonic_n_remain", "onekg_frac_remain", "gnomad_frac_remain", "clinvar_frac_remain",
                                   "dbsnp_frac_remain", "pon_frac_remain", "hom_frac_remain", "het_frac_remain",
                                   "nonexonic_frac_remain")) {
          report[["content"]][[analysis_element]][["variant_statistic"]][[successive_filter]] <- 0
        }
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
#' @param analysis_element section of PCGR report

update_report <- function(report, report_data, analysis_element = "snv_indel") {

  if (!is.null(report_data) & !is.null(report[["content"]][[analysis_element]])) {
    for (report_elem in names(report_data)) {
      if (!is.null(report_data[[report_elem]]) & !is.null(report[["content"]][[analysis_element]][[report_elem]])) {
         report[["content"]][[analysis_element]][[report_elem]] <- report_data[[report_elem]]
      }
    }
  }
  return(report)
}

set_report_metadata <- function(config, pcgr_data,
                                cpsr_version = NULL,
                                pcgr_version = NULL,
                                virtual_panel_id = -1,
                                sample_name = NULL,
                                report_type = NULL,
                                custom_bed = NULL){

  if(is.null(pcgr_data)){
    return(NULL)
  }

  report_metadata <- list()
  report_metadata[["pcgr_db_release"]] <- pcgr_data[["release_notes"]]
  report_metadata[["pcgr_version"]] <- pcgr_version
  report_metadata[["cpsr_version"]] <- cpsr_version
  report_metadata[["genome_assembly"]] <- pcgr_data[["assembly"]][["grch_name"]]
  report_metadata[["sample_name"]] <- sample_name
  report_metadata[["report_type"]] <- report_type
  report_metadata[["config"]] <- config
  report_metadata[['color_palette']] <- pcgrr::color_palette
  report_metadata[['color_none']] <- pcgrr::color_palette[['none']][['values']][1]
  report_metadata[['color_value_box']] <- pcgrr::color_palette[['report_color']][['values']][1]
  if(report_type == "somatic" & !is.null(report_metadata[['config']][['assay_properties']])){
    if(report_metadata[['config']][['assay_properties']][['vcf_tumor_only']] == T){
      report_metadata[['color_value_box']] <- pcgrr::color_palette[['report_color']][['values']][2]
    }
  }
  report_metadata[["phenotype_ontology"]] <- list()
  report_metadata[["phenotype_ontology"]][["oncotree"]] <- pcgr_data[["phenotype_ontology"]][["oncotree"]]
  report_metadata[["phenotype_ontology"]][["oncotree_query"]] <- NULL


  if (virtual_panel_id >= 0) {
    report_metadata[["gene_panel"]] <- list()
    report_metadata[["gene_panel"]][["genes"]] <- pcgr_data[["virtual_gene_panels"]] %>%
      dplyr::filter(id == virtual_panel_id) %>%
      dplyr::select(symbol, confidence_level, panel_name, panel_id, panel_version, panel_url, entrezgene, genename)
    if (config[['diagnostic_grade_only']] == T) {
      report_metadata[["gene_panel"]][["genes"]] <- report_metadata[["gene_panel"]][["genes"]] %>%
        dplyr::filter(confidence_level == 3 | confidence_level == 4)
    }
    report_metadata[["gene_panel"]][["name"]] <-
      paste0(unique(report_metadata[["gene_panel"]][["genes"]]$panel_name), " - v",
             unique(report_metadata[["gene_panel"]][["genes"]]$panel_version))
    report_metadata[["gene_panel"]][["url"]] <- unique(report_metadata[["gene_panel"]][["genes"]]$panel_url)
    report_metadata[["gene_panel"]][["confidence"]] <- "Diagnostic-grade/Borderline/Low evidence (GREEN/AMBER/RED)"
    if (diagnostic_grade_only == 1 & virtual_panel_id != 0) {
      report_metadata[["gene_panel"]][["confidence"]] <- "Diagnostic-grade only (GREEN)"
    }
    if (virtual_panel_id == 0) {
      report_metadata[["gene_panel"]][["confidence"]] <- "Exploratory geneset (research)"
    }
  }else{
    if (!is.null(custom_bed)) {
      target_genes <- pcgrr::custom_bed_genes(custom_bed, pcgr_data = pcgr_data)
      if (nrow(target_genes) > 0) {

        target_genes <- target_genes %>%
          dplyr::mutate(confidence_level = -1,
                        panel_name = report_metadata[["config"]][["custom_panel"]][["name"]],
                        panel_id = NA, panel_version = NA,
                        panel_url = report_metadata[["config"]][["custom_panel"]][["url"]])

        report_metadata[["gene_panel"]] <- list()
        report_metadata[["gene_panel"]][["genes"]] <- target_genes
        report_metadata[["gene_panel"]][["name"]] <- report_metadata[["config"]][["custom_panel"]][["name"]]
        report_metadata[["gene_panel"]][["confidence"]] <- "User-defined panel (custom geneset from panel 0)"
      }else{
        rlogging::warning(paste0("Custom BED file (", custom_bed, ") does not contain regions that overlap protein-coding transcripts (regulatory/exonic)"))
        rlogging::message("Quitting report generation")
      }
    }
  }
  return(report_metadata)
}
