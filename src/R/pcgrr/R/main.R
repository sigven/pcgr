#' Function that generates all contents of the cancer genome report (PCGR)
#'
#' @param project_directory name of project directory
#' @param pcgr_data List object with multiple PCGR data bundle annotations
#' @param config Object with PCGR configuration parameters
#' @param tier_model Variant tier model
#'
#' @export

generate_pcgr_report <-
  function(project_directory = NULL,
           pcgr_data = NULL,
           config = NULL,
           tier_model = "pcgr_acmg") {

    ## Assert argument values - config
    invisible(assertthat::assert_that(
      !is.null(config) & is.list(config),
      msg = "Argument 'config' must be a non-NULL list object"))

    ## Assert argument values - pcgr_data
    invisible(assertthat::assert_that(
      !is.null(pcgr_data) & is.list(pcgr_data),
      msg = "Argument 'pcgr_data' must be a non-NULL list object"))

    ## Assert argument values - tier_model
    invisible(assertthat::assert_that(
      is.character(tier_model),
      msg = "Argument 'tier_model' must be of type character"))

    ## Assert argument values - tier_model
    invisible(assertthat::assert_that(
      tier_model == "pcgr_acmg" | tier_model == "cpsr",
      msg = "Argument 'tier_model' must be either 'pcgr_acmg' or 'cpsr"))

    query_vcf2tsv <- config[['required_args']][['query_vcf2tsv']]
    sample_name <- config[['required_args']][['sample_name']]
    cna_segments_tsv <- config[['required_args']][['query_cna']]
    pcgr_version <- config[['required_args']][['pcgr_version']]
    cpsr_report_fname <- config[['required_args']][['cpsr_report']]

    log4r_info(paste0("Initializing PCGR report - sample ", sample_name))
    log4r_info("------")


    if (!is.null(cna_segments_tsv)) {
      invisible(assertthat::assert_that(
        file.exists(cna_segments_tsv),
        msg = paste0("Filename provided for argument 'cna_segments_tsv' (",
                     cna_segments_tsv, ") does not exist")))
      invisible(assertthat::assert_that(
        file.size(cna_segments_tsv) > 0,
        msg = paste0("File provided for argument 'cna_segments_tsv' (",
                     cna_segments_tsv, ") has a filesize of zero")))
    }
    if (!is.null(query_vcf2tsv)) {
      invisible(assertthat::assert_that(
        file.exists(query_vcf2tsv),
        msg = paste0("Filename provided for argument 'query_vcf2tsv' (",
                     query_vcf2tsv, ") does not exist")))
      invisible(assertthat::assert_that(
        file.size(query_vcf2tsv) > 0,
        msg = paste0("File provided for argument 'query_vcf2tsv' (",
                     query_vcf2tsv, ") has a filesize of zero")))

    }
    if (!is.null(cpsr_report_fname)) {
      invisible(assertthat::assert_that(
        file.exists(cpsr_report_fname),
        msg = paste0("Filename provided for argument 'cpsr_report' (",
                     cpsr_report_fname, ") does not exist")))
      invisible(assertthat::assert_that(
        file.size(cpsr_report_fname) > 0,
        msg = paste0("File provided for argument 'cpsr_report' (",
                     cpsr_report_fname, ") has a filesize of zero")))

    }

    pcg_report <- pcgrr::init_report(
      config = config,
      class = NULL,
      pcgr_data = pcgr_data
    )

    ## Set output and temporary file names
    sample_fname_pattern <-
      paste(sample_name, tier_model,
            pcgr_data[["assembly"]][["grch_name"]], sep = ".")

    fnames <- list()
    fnames[["maf"]] <-
      file.path(project_directory, paste0(sample_fname_pattern, ".maf"))
    fnames[["maf_tmp"]] <-
      file.path(project_directory, paste0(sample_fname_pattern, ".tmp.maf"))
    fnames[["vcf_mp"]] <-
      file.path(project_directory, paste0(sample_fname_pattern, ".mp_input.vcf"))


    ## Retrieve relevant clinical trials for the tumor type in question

    if (config[["clinicaltrials"]][["run"]] == T) {
      pcg_report_trials <-
        pcgrr::generate_report_data_trials(
          pcgr_data = pcgr_data,
          sample_name = sample_name,
          config = config)
      ## Update genome report with trial data
      pcg_report <-
        pcgrr::update_report(pcg_report, pcg_report_trials,
                             a_elem = "clinicaltrials")
    }



    ## Load sample calls (check integrity, filter variants for
    ## depth/allelic fraction, append links etc)
    sample_calls <-
      pcgrr::get_calls(
        query_vcf2tsv,
        pcgr_data,
        sample_name,
        config,
        oncotree =
          pcg_report[["metadata"]][["phenotype_ontology"]][["oncotree_query"]],
        maf_filenames = fnames)


    if (nrow(sample_calls) > 0) {

      assay_props <-
        pcg_report[["metadata"]][["config"]][["assay_props"]]

      ## Perform analyses in tumor-only mode
      if (assay_props[["vcf_tumor_only"]] == TRUE) {
        pcg_report_tumor_only <-
          pcgrr::generate_report_data_tumor_only(sample_calls,
                                                 sample_name, config)

        ## Generate data for SNVs/InDels
        ## -
        pcg_report_snv_indel_filtered <-
          pcgrr::generate_report_data_snv_indel(
            pcg_report_tumor_only[["variant_set"]][["filtered"]],
            pcgr_data,
            sample_name,
            config,
            callset = "germline-filtered callset",
            tier_model = tier_model)

        pcg_report_tumor_only[["upset_data"]] <-
          pcgrr::make_upset_plot_data(
            pcg_report_tumor_only$variant_set$tsv_unfiltered, config)
        num_upset_sources <- 0
        for (c in colnames(pcg_report_tumor_only[["upset_data"]])) {
          if (c != "VAR_ID") {
            if (sum(pcg_report_tumor_only[["upset_data"]][, c]) > 0) {
              num_upset_sources <- num_upset_sources + 1
            }
          }
        }
        if (num_upset_sources >= 2) {
          pcg_report_tumor_only[["upset_plot_valid"]] <- TRUE
        }

        ## Update genome report with SNV/InDels (display, tiers etc)
        pcg_report <-
          pcgrr::update_report(pcg_report, pcg_report_snv_indel_filtered,
                               a_elem = "snv_indel")
        pcg_report <-
          pcgrr::update_report(pcg_report, pcg_report_tumor_only,
                               a_elem = "tumor_only")

        ## Generate data for rainfall plot (SNVs)
        pcg_report_rainfall <-
          pcgrr::generate_report_data_rainfall(
            pcg_report$content$snv_indel$variant_set$tsv,
            build = pcg_report$metadata$genome_assembly)
        ## Update genome report
        pcg_report <-
          pcgrr::update_report(pcg_report, pcg_report_rainfall,
                               a_elem = "rainfall")

      }else{
        ## Generate report data for SNVs/InDels
        pcg_report_snv_indel <-
          pcgrr::generate_report_data_snv_indel(
            sample_calls,
            pcgr_data,
            sample_name,
            config, tier_model = tier_model)

        ## Update genome report
        pcg_report <- pcgrr::update_report(
          pcg_report, pcg_report_snv_indel,
          a_elem = "snv_indel")
      }
      rm(sample_calls)

      ## Estimate contribution of mutational signatures
      if (pcg_report[["metadata"]][["config"]][["msigs"]][["run"]] == T) {

        if(NROW(pcg_report$content$snv_indel$variant_set$tsv) > 0){
          pcgrr::write_processed_vcf(
            calls = pcg_report$content$snv_indel$variant_set$tsv,
            sample_name = sample_name,
            output_directory = project_directory,
            vcf_fname = fnames[["vcf_mp"]])

          pcg_report_signatures <-
            pcgrr::generate_report_data_signatures_mp(
              vcf_fname = paste0(fnames[["vcf_mp"]], ".gz"),
              pcgr_data,
              sample_name,
              config,
              type_specific =
                !config[["msigs"]][["all_reference_signatures"]])

          ## Update genome report with signature info
          pcg_report <- pcgrr::update_report(
            pcg_report, pcg_report_signatures,
            a_elem = "m_signature_mp")
        }

        ## Generate report data for rainfall plot
        pcg_report_rainfall <-
          pcgrr::generate_report_data_rainfall(
            variant_set = pcg_report$content$snv_indel$variant_set$tsv,
            build = pcg_report$metadata$genome_assembly)

        ## Update genome report
        pcg_report <-
          pcgrr::update_report(pcg_report,
                               pcg_report_rainfall,
                               a_elem = "rainfall")

        ## Generate report data for kataegis events (for WES/WGS runs)
        if (stringr::str_detect(
          assay_props[["type"]],
          "WGS|WES")) {
          pcg_report_kataegis <-
            pcgrr::generate_report_data_kataegis(
              variant_set = pcg_report$content$snv_indel$variant_set$tsv,
              sample_name = sample_name,
              build = pcg_report$metadata$genome_assembly)
          ## Update genome report
          pcg_report <- pcgrr::update_report(
            pcg_report, pcg_report_kataegis,
            a_elem = "kataegis")
        }
      }

      ## If assay is Tumor-Control and WES/WGS - perform MSI prediction
      if (pcg_report[["metadata"]][["config"]][["msi"]][["run"]] == T &
          stringr::str_detect(assay_props[["type"]], "WGS|WES") &
          assay_props[["vcf_tumor_only"]] == FALSE) {
        pcg_report_msi <-
          pcgrr::generate_report_data_msi(
            pcg_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]],
            pcgr_data, sample_name, config)

        ## Update genome report with MSI info
        pcg_report <-
          pcgrr::update_report(pcg_report,
                               pcg_report_msi,
                               a_elem = "msi")
      }

      ## Generate report contents for analysis of mutational burden (TMB)
      if (pcg_report[["metadata"]][["config"]][["tmb"]][["run"]] == T) {
        pcg_report_tmb <-
          pcgrr::generate_report_data_tmb(
            pcg_report[["content"]][["snv_indel"]][["variant_set"]][["tsv"]],
            pcgr_data, sample_name, config)

        ## Update genome report with TMB info
        pcg_report <- pcgrr::update_report(
          pcg_report, pcg_report_tmb,
          a_elem = "tmb")
      }
    }else{
      pcg_report[["content"]][["snv_indel"]][["zero"]] <- TRUE
      pcg_report[["metadata"]][["config"]][["other"]][["list_noncoding"]] <- FALSE
    }

    if (!is.null(cpsr_report_fname)) {
      pcg_report[["content"]][["cpsr"]][['eval']] <- TRUE

      pcg_report[['content']][['cpsr']][['report']] <-
        jsonlite::fromJSON(
          gzfile(cpsr_report_fname)
        )

      ## append report elements in pcg_report[['content']][['cpsr]][['cpsr_json']]
    }

    if (!is.null(cna_segments_tsv)) {
      pcg_report_cna <-
        pcgrr::generate_report_data_cna(
          cna_segments_tsv,
          pcgr_data,
          sample_name,
          config,
          transcript_overlap_pct = config[["cna"]][["cna_overlap_pct"]])
      pcg_report <-
        pcgrr::update_report(pcg_report,
                             pcg_report_cna,
                             a_elem = "cna")

    }

    pcg_report_value_box <- pcgrr::generate_report_data_value_box(
      pcg_report, pcgr_data, sample_name, config)
    pcg_report <- pcgrr::update_report(
      pcg_report, pcg_report_value_box,
      a_elem = "value_box")

    for (elem in c("tier1", "tier2", "tier3", "tier4")) {
      stat <- paste0("n_", elem)
      pcg_report[["content"]][["snv_indel"]][["v_stat"]][[stat]] <-
        nrow(pcg_report[["content"]][["snv_indel"]][["variant_set"]][[elem]])
      pcg_report[["content"]][["snv_indel"]][["variant_set"]][[elem]] <- NULL
    }
    pcg_report[["content"]][["snv_indel"]][["variant_set"]][["noncoding"]] <- NULL
    pcg_report[["content"]][["snv_indel"]][["variant_set"]][["coding"]] <- NULL
    pcg_report[["content"]][["snv_indel"]][["variant_set"]][["all"]] <- NULL
    if(!is.null(pcg_report[["content"]][["tumor_only"]])){
      pcg_report[["content"]][["snv_indel"]][["variant_set"]][["tsv_unfiltered"]] <-
        pcg_report[["content"]][["tumor_only"]][["variant_set"]][["tsv_unfiltered"]]
      pcg_report[["content"]][["tumor_only"]][["variant_set"]][["tsv_unfiltered"]] <- NULL
      pcg_report[["content"]][["tumor_only"]][["variant_set"]][["filtered"]] <- NULL
    }
    pcg_report[["content"]][["snv_indel"]][["variant_set"]][["all"]] <- NULL
    pcg_report[["content"]][["cna"]][["variant_set"]][["cna_print"]] <- NULL
    pcg_report[["metadata"]][["phenotype_ontology"]] <- list()
    gc()

    # if (!is.null(cna_plot) && cna_plot != "None") {
    #   pcg_report[["content"]][["cna_plot"]][["png"]] <- cna_plot
    #   pcg_report[["content"]][["cna_plot"]][["eval"]] <- TRUE
    # }
    return(pcg_report)
  }


#' Function that generates tiered variant sets for SNVs/InDels
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param config Object with PCGR configuration parameters
#' @param callset type of calls
#' @param biomarker_mapping_stringency quality level for biomarkers
#' @param tier_model tier model (pcgr_acmg)
#'
#' @return pcg_report_data data frame with all report elements
#'
#' @export
generate_report_data_snv_indel <- function(
  sample_calls,
  pcgr_data,
  sample_name,
  config,
  callset = "somatic calls",
  biomarker_mapping_stringency = 1,
  tier_model = "pcgr_acmg") {

  log4r_info("------")
  log4r_info(
    paste0("Generating data for tiered cancer genome report - ",
           callset, " tier model '", tier_model, "'"))

  pcg_report_snv_indel <- pcgrr::init_report(config = config,
                                             class = "snv_indel")
  pcg_report_snv_indel[["eval"]] <- TRUE
  pcg_report_snv_indel[["variant_set"]][["all"]] <- sample_calls

  ## Get basic variant statistics (type, coding status)
  call_stats <- pcgrr::variant_stats_report(sample_calls,
                                            name = "v_stat")
  for (stat in c("n", "n_snv", "n_indel", "n_coding", "n_noncoding")) {
    pcg_report_snv_indel[["v_stat"]][[stat]] <-
      call_stats[["v_stat"]][[stat]]
  }
  log4r_info(
    paste0("Number of protein-coding variants: ",
           pcg_report_snv_indel[["v_stat"]][["n_coding"]]))

  annotation_tags <- pcgr_data[["annotation_tags"]]

  ## Add custom annotation tags to lists of tags to display
  if (!is.null(config[["preserved_info_tags"]])) {
    if (config[["preserved_info_tags"]] != "None") {
      tags <- stringr::str_split(
        config[["preserved_info_tags"]],
        pattern = ",")[[1]]
      for (t in tags) {
        t <- stringr::str_trim(t)
        if (t %in% colnames(sample_calls)) {
          annotation_tags[["all"]] <-
            c(annotation_tags[["all"]], t)
        }
      }
    }
  }

  ## remove REGULATORY_ANNOTATION from display tags
  ## if regulatory annotation is not turned on
  if(!is.null(config)){
    if(config[['other']][['vep_regulatory']] == FALSE){
      for(e in c('all','tier4_display','tier5_display',
                 'tsv')){
        annotation_tags[[e]] <-
          annotation_tags[[e]][!annotation_tags[[e]] == "REGULATORY_ANNOTATION"]
      }
    }
  }

  if (pcg_report_snv_indel[["v_stat"]][["n"]] > 0) {

    tumor_type <- config[["t_props"]][["tumor_type"]]

    ## load all clinical evidence items ()
    eitems_any_tt <- pcgrr::load_eitems(
      eitems_raw = pcgr_data$biomarkers,
      alteration_type = "MUT",
      ontology =
        pcgr_data$phenotype_ontology$oncotree,
      origin = "Somatic",
      tumor_type_specificity = "any")

    ## Get all clinical evidence items that overlap
    ## query set (NOT tumor-type specific)
    biomarker_hits_snv_indels_any <-
      pcgrr::get_clin_assocs_snv_indel(
        sample_calls = pcg_report_snv_indel[["variant_set"]][["all"]],
        annotation_tags = annotation_tags,
        eitems = eitems_any_tt)

    ## Assign putative TIER 2 variant set
    pcg_report_snv_indel[["clin_eitem"]][["any_ttype"]] <-
      biomarker_hits_snv_indels_any$clin_eitem
    pcg_report_snv_indel[["variant_set"]][["tier2"]] <-
      biomarker_hits_snv_indels_any$variant_set

    ## Get all clinical evidence items that
    ## overlap query set (if tumor type is specified)
    if (tumor_type != "Cancer, NOS") {

      ## load tumor-type specific evidence items ()
      eitems_specific_tt <- pcgrr::load_eitems(
        eitems_raw = pcgr_data$biomarkers,
        alteration_type = "MUT",
        ontology =
          pcgr_data$phenotype_ontology$oncotree,
        origin = "Somatic",
        tumor_type_specificity = "specific",
        tumor_type = tumor_type)

      biomarker_hits_snv_indels_specific <-
        pcgrr::get_clin_assocs_snv_indel(
          sample_calls = pcg_report_snv_indel[["variant_set"]][["all"]],
          annotation_tags = annotation_tags,
          eitems = eitems_specific_tt)

      ## Assign putative TIER 1 variant set
      pcg_report_snv_indel[["clin_eitem"]][["specific_ttype"]] <-
        biomarker_hits_snv_indels_specific$clin_eitem
      pcg_report_snv_indel[["variant_set"]][["tier1"]] <-
        biomarker_hits_snv_indels_specific$variant_set
    }

    ## Remove potential overlap/redundancies and assign final
    ## TIER1/TIER2 classification
    pcg_report_snv_indel <- pcgrr::assign_tier1_tier2_acmg(pcg_report_snv_indel)
    tier12 <- rbind(pcg_report_snv_indel[["disp"]][["tier1"]],
                    pcg_report_snv_indel[["disp"]][["tier2"]])

    ## Determine TIER 3 variant set: coding mutations in
    ## oncogenes/tumor suppressors/cancer census genes
    pcg_report_snv_indel[["variant_set"]][["tier3"]] <-
      dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]],
                    dplyr::one_of(annotation_tags[["all"]])) %>%
      dplyr::filter(.data$CODING_STATUS == "coding") %>%
      dplyr::filter(.data$ONCOGENE == TRUE | .data$TUMOR_SUPPRESSOR == TRUE)
    if (nrow(tier12) > 0 &
        nrow(pcg_report_snv_indel[["variant_set"]][["tier3"]]) > 0) {
      pcg_report_snv_indel[["variant_set"]][["tier3"]] <-
        dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["tier3"]],
                         tier12, by = c("GENOMIC_CHANGE"))
    }
    tier123 <- tier12
    if (nrow(pcg_report_snv_indel[["variant_set"]][["tier3"]]) > 0) {
      pcg_report_snv_indel[["variant_set"]][["tier3"]] <-
        pcg_report_snv_indel[["variant_set"]][["tier3"]] %>%
        dplyr::arrange(dplyr::desc(.data$OPENTARGETS_RANK))
      tier123 <- tier12 %>%
        dplyr::bind_rows(
          dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier3"]],
                        .data$GENOMIC_CHANGE)) %>%
        dplyr::distinct()
      pcg_report_snv_indel[["disp"]][["tier3"]][["proto_oncogene"]] <-
        dplyr::select(
          pcg_report_snv_indel[["variant_set"]][["tier3"]],
          dplyr::one_of(annotation_tags[["tier3_display"]])) %>%
        dplyr::filter(.data$ONCOGENE == TRUE &
                        (is.na(.data$TUMOR_SUPPRESSOR) |
                           .data$TUMOR_SUPPRESSOR == FALSE))
      pcg_report_snv_indel[["disp"]][["tier3"]][["tumor_suppressor"]] <-
        dplyr::select(
          pcg_report_snv_indel[["variant_set"]][["tier3"]],
          dplyr::one_of(annotation_tags[["tier3_display"]])) %>%
        dplyr::filter(.data$TUMOR_SUPPRESSOR == TRUE)
    }

    ## Determine TIER 4: Other coding mutations
    pcg_report_snv_indel[["variant_set"]][["tier4"]] <-
      dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]],
                    dplyr::one_of(annotation_tags[["all"]])) %>%
      dplyr::filter(.data$CODING_STATUS == "coding")
    if (nrow(tier123) > 0 &
        nrow(pcg_report_snv_indel[["variant_set"]][["tier4"]]) > 0) {
      pcg_report_snv_indel[["variant_set"]][["tier4"]] <-
        dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["tier4"]],
                         tier123, by = c("GENOMIC_CHANGE"))
    }
    if (nrow(pcg_report_snv_indel[["variant_set"]][["tier4"]]) > 0) {
      pcg_report_snv_indel[["variant_set"]][["tier4"]] <-
        pcg_report_snv_indel[["variant_set"]][["tier4"]] %>%
        dplyr::arrange(dplyr::desc(.data$OPENTARGETS_RANK))
      pcg_report_snv_indel[["disp"]][["tier4"]] <-
        dplyr::select(
          pcg_report_snv_indel[["variant_set"]][["tier4"]],
          dplyr::one_of(annotation_tags[["tier4_display"]]))
    }

    ## Determine non-coding mutation set
    pcg_report_snv_indel[["variant_set"]][["noncoding"]] <-
      dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]],
                    dplyr::one_of(annotation_tags[["all"]])) %>%
      dplyr::filter(.data$CODING_STATUS == "noncoding")
    if (nrow(pcg_report_snv_indel[["variant_set"]][["noncoding"]]) > 0) {
      if (nrow(tier123) > 0) {
        pcg_report_snv_indel[["variant_set"]][["noncoding"]] <-
          dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["noncoding"]],
                           tier123,
                           by = c("GENOMIC_CHANGE"))
      }
      pcg_report_snv_indel[["variant_set"]][["noncoding"]] <-
        pcg_report_snv_indel[["variant_set"]][["noncoding"]] %>%
        dplyr::arrange(dplyr::desc(.data$OPENTARGETS_RANK))
      pcg_report_snv_indel[["disp"]][["noncoding"]] <-
        dplyr::select(
          pcg_report_snv_indel[["variant_set"]][["noncoding"]],
          dplyr::one_of(annotation_tags[["tier5_display"]]))
    }

    ## Make TSV content with variant set
    pcg_report_snv_indel[["v_stat"]][["n_noncoding"]] <-
      pcg_report_snv_indel[["variant_set"]][["noncoding"]] %>% nrow()
    pcg_report_snv_indel[["variant_set"]][["tsv"]] <-
      pcgrr::generate_tier_tsv(
        pcg_report_snv_indel[["variant_set"]],
        config,
        annotation_tags,
        sample_name = sample_name)

  }

  return(pcg_report_snv_indel)

}


#' Function that generates germline-filtered callset and PCGR
#' report statistics for a given tumor-only callsets
#'
#' @param unfiltered_sample_calls variant calls
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#'
#' @export
generate_report_data_tumor_only <-
  function(unfiltered_sample_calls,
           sample_name,
           pcgr_config) {

  sample_calls <- unfiltered_sample_calls
  gline_filter_stats <- list()
  for (m in c("remain_post_onekg",
              "remain_post_gnomad",
              "remain_post_clinvar",
              "remain_post_dbsnp",
              "remain_post_pon",
              "remain_post_nonexonic",
              "remain_post_hom",
              "remain_post_het")) {
    gline_filter_stats[m] <- 0
  }

  ## initiate report
  pcg_report_to <-
    pcgrr::init_report(config = pcgr_config,
                       class = "tumor_only")

  ## assign evidence tags for germline/somatic state of variants,
  ## partially based on user-defined criteria
  ## (population allele frequency thresholds)
  vcalls <-
    pcgrr::assign_somatic_germline_evidence(sample_calls, pcgr_config)

  ## assign somatic classification based on accumulation
  ## of evidence tags and user-defined options
  vcalls <-
    pcgrr::assign_somatic_classification(vcalls, pcgr_config)

  ## Assign statistics to successive filtering levels for
  ## different evidence criteria
  ## excluded germline calls found in 1000 Genomes Project
  gline_filter_stats[["remain_post_onekg"]] <-
    nrow(vcalls) -
    nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_1KG", ])
  log4r_info(paste0("Excluding coinciding germline variants in ",
                           "1000 Genomes Project populations"))
  log4r_info(paste0("Total sample calls remaining: ",
                           gline_filter_stats[["remain_post_onekg"]]))

  ## excluded germline calls found in gnomAD
  gline_filter_stats[["remain_post_gnomad"]] <-
    gline_filter_stats[["remain_post_onekg"]] -
    nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_GNOMAD", ])
  log4r_info(
    paste0("Excluding coinciding germline variants in any ",
           "population in the genome aggregation database (gnomAD)"))
  log4r_info(paste0("Total sample calls remaining: ",
                           gline_filter_stats[["remain_post_gnomad"]]))

  ## excluded germline calls found in ClinVar
  gline_filter_stats[["remain_post_clinvar"]] <-
    gline_filter_stats[["remain_post_gnomad"]] -
    nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_CLINVAR", ])
  log4r_info(paste0("Excluding coinciding germline variants in ClinVar"))
  log4r_info(paste0("Total sample calls remaining: ",
                           gline_filter_stats[["remain_post_clinvar"]]))


  ## excluded germline calls found in panel of normals (if provided)
  gline_filter_stats[["remain_post_pon"]] <-
    gline_filter_stats[["remain_post_clinvar"]]
  if (pcgr_config[["tumor_only"]][["exclude_pon"]] == TRUE) {
    gline_filter_stats[["remain_post_pon"]] <-
      gline_filter_stats[["remain_post_pon"]] -
      nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_PON", ])
    log4r_info(
      paste0("Excluding putative germline variants found in calls ",
      "from panel-of-normals (PON)"))
    log4r_info(
      paste0("Total sample calls remaining: ",
             gline_filter_stats[["remain_post_pon"]]))
  }

  ## excluded germline calls found with 100% allelic fraction
  ## (likely homozygous germline variants)
  gline_filter_stats[["remain_post_hom"]] <-
    gline_filter_stats[["remain_post_pon"]]
  if (pcgr_config[["tumor_only"]][["exclude_likely_hom_germline"]] == TRUE) {
    gline_filter_stats[["remain_post_hom"]] <-
      gline_filter_stats[["remain_post_hom"]] -
      nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_HOMOZYGOUS", ])
    log4r_info(
      paste0("Excluding likely homozygous germline variants found ",
             "as variants with 100% allelic fraction"))
    log4r_info(paste0("Total sample calls remaining: ",
                             gline_filter_stats[["remain_post_hom"]]))
  }

  ## excluded germline calls found as likely heterozygous germline variants
  gline_filter_stats[["remain_post_het"]] <-
    gline_filter_stats[["remain_post_hom"]]
  if (pcgr_config[["tumor_only"]][["exclude_likely_het_germline"]] == TRUE) {
    gline_filter_stats[["remain_post_het"]] <-
      gline_filter_stats[["remain_post_het"]] -
      nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_HETEROZYGOUS", ])
    log4r_info(paste0(
      "Excluding likely heterozygous germline variants found as variants ",
      "with 40-60% allelic fraction and recorded in gnomAD + dbSNP"))
    log4r_info(paste0("Total sample calls remaining: ",
                             gline_filter_stats[["remain_post_het"]]))
  }

  ## excluded calls with dbSNP germline status (if set in config)
  gline_filter_stats[["remain_post_dbsnp"]] <-
    gline_filter_stats[["remain_post_het"]]
  if (pcgr_config[["tumor_only"]][["exclude_dbsnp_nonsomatic"]] == TRUE) {

    log4r_info(
      paste0("Excluding non-somatically associated dbSNP variants ",
             "(dbSNP - not recorded as somatic in DoCM/ClinVar",
             "and not registered in COSMIC or found in TCGA"))

    gline_filter_stats[["remain_post_dbsnp"]] <-
      gline_filter_stats[["remain_post_dbsnp"]] -
      nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_DBSNP", ])
    log4r_info(paste0("Total sample calls remaining: ",
                             gline_filter_stats[["remain_post_dbsnp"]]))
  }

  unfiltered_sample_calls <- vcalls
  vcalls <- vcalls %>%
    dplyr::filter(.data$SOMATIC_CLASSIFICATION == "SOMATIC")

  gline_filter_stats[["remain_post_nonexonic"]] <-
    gline_filter_stats[["remain_post_dbsnp"]]
  if (pcgr_config[["tumor_only"]][["exclude_nonexonic"]] == TRUE) {
    log4r_info(paste0("Excluding non-exonic variants"))
    vcalls <- dplyr::filter(vcalls, .data$EXONIC_STATUS == "exonic")
    log4r_info(paste0("Total sample calls remaining: ",
                             nrow(vcalls)))
    gline_filter_stats[["remain_post_nonexonic"]] <- nrow(vcalls)
  }

  pcg_report_to[["eval"]] <- TRUE
  pcg_report_to[["variant_set"]][["tsv_unfiltered"]] <- unfiltered_sample_calls %>%
    dplyr::select(.data$GENOMIC_CHANGE,
                  .data$VAR_ID,
                  .data$DP_TUMOR,
                  .data$AF_TUMOR,
                  .data$SYMBOL,
                  .data$EXONIC_STATUS,
                  .data$CONSEQUENCE,
                  .data$STATUS_PON,
                  .data$STATUS_LIKELY_GERMLINE_HOMOZYGOUS,
                  .data$STATUS_LIKELY_GERMLINE_HETEROZYGOUS,
                  .data$STATUS_DBSNP_GERMLINE,
                  .data$STATUS_POPFREQ_1KG_ABOVE_TOLERATED,
                  .data$STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED,
                  .data$STATUS_CLINVAR_GERMLINE,
                  .data$SOMATIC_CLASSIFICATION)
  pcg_report_to[["variant_set"]][["filtered"]] <- vcalls
  pcg_report_to[["v_stat"]][["unfiltered_n"]] <-
    nrow(unfiltered_sample_calls)
  pcg_report_to[["v_stat"]][["onekg_n_remain"]] <-
    gline_filter_stats[["remain_post_onekg"]]
  pcg_report_to[["v_stat"]][["gnomad_n_remain"]] <-
    gline_filter_stats[["remain_post_gnomad"]]
  pcg_report_to[["v_stat"]][["clinvar_n_remain"]] <-
    gline_filter_stats[["remain_post_clinvar"]]
  pcg_report_to[["v_stat"]][["pon_n_remain"]] <-
    gline_filter_stats[["remain_post_pon"]]
  pcg_report_to[["v_stat"]][["hom_n_remain"]] <-
    gline_filter_stats[["remain_post_hom"]]
  pcg_report_to[["v_stat"]][["het_n_remain"]] <-
    gline_filter_stats[["remain_post_het"]]
  pcg_report_to[["v_stat"]][["dbsnp_n_remain"]] <-
    gline_filter_stats[["remain_post_dbsnp"]]
  pcg_report_to[["v_stat"]][["nonexonic_n_remain"]] <-
    gline_filter_stats[["remain_post_nonexonic"]]
  for (db_filter in c("onekg", "gnomad", "dbsnp", "pon",
                      "clinvar", "hom", "het", "nonexonic")) {
    if (pcg_report_to[["v_stat"]][[paste0(db_filter, "_n_remain")]] > 0 &
        pcg_report_to[["v_stat"]][["unfiltered_n"]] > 0) {
      pcg_report_to[["v_stat"]][[paste0(db_filter, "_frac_remain")]] <-
        round((as.numeric(pcg_report_to[["v_stat"]][[paste0(db_filter,
                                                            "_n_remain")]]) /
                 pcg_report_to[["v_stat"]][["unfiltered_n"]]) * 100, digits = 2)
    }
  }
  return(pcg_report_to)

}

#' Function that annotates CNV segment files
#'
#' @param cna_segments_tsv CNV file name with chromosomal log(2)-ratio segments
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param transcript_overlap_pct required aberration overlap fraction
#' (percent) for reported transcripts (default 100 percent)
#'
#' @export
generate_report_data_cna <-
  function(cna_segments_tsv,
           pcgr_data,
           sample_name,
           pcgr_config,
           transcript_overlap_pct = 100) {

    invisible(
      assertthat::assert_that(
        file.exists(cna_segments_tsv),
        msg = paste0("File 'cna_segments_tsv' (",
                     cna_segments_tsv, ") does not exist")))
    pcg_report_cna <- pcgrr::init_report(config = pcgr_config,
                                         class = "cna")
    log_r_homdel <- pcgr_config[["cna"]][["log_r_homdel"]]
    log_r_gain <- pcgr_config[["cna"]][["log_r_gain"]]
    tumor_type <- pcgr_config[["t_props"]][["tumor_type"]]
    MEGABASE <- 1000000

    log4r_info("------")
    log4r_info(paste0("Generating report data for copy number segment file ",
                      cna_segments_tsv))

    ## READ INPUT FILE, VALIDATE INPUT CHROMOSOMES AND SEGMENTS, ADD CYTOBAND INFO
    cna_df <- utils::read.table(file = cna_segments_tsv, header = T,
                         stringsAsFactors = F, sep = "\t",
                         comment.char = "", quote = "") %>%
      dplyr::rename(chromosome = .data$Chromosome, LogR = .data$Segment_Mean,
                    segment_start = .data$Start, segment_end = .data$End) %>%
      dplyr::distinct() %>%
      dplyr::select(.data$chromosome, .data$LogR, .data$segment_start, .data$segment_end) %>%
      dplyr::mutate(
        chromosome = stringr::str_replace(.data$chromosome, "^chr", "")) %>%
      pcgrr::get_valid_chromosomes(
        chromosome_column = "chromosome",
        bsg = pcgr_data[["assembly"]][["bsg"]]) %>%
      pcgrr::get_valid_chromosome_segments(
        genome_assembly = pcgr_data[["assembly"]][["grch_name"]],
        bsg = pcgr_data[["assembly"]][["bsg"]]) %>%
      dplyr::filter(!is.na(.data$LogR)) %>%
      dplyr::mutate(LogR = round(as.numeric(.data$LogR), digits = 3)) %>%
      dplyr::mutate(SEGMENT_ID = paste0(.data$chromosome, ":",
                                        .data$segment_start, "-",
                                        .data$segment_end)) %>%
      pcgrr::get_cna_cytoband(pcgr_data = pcgr_data) %>%
      dplyr::mutate(SAMPLE_ID = sample_name) %>%
      pcgrr::append_ucsc_segment_link(
        hgname = pcgr_data[["assembly"]][["hg_name"]],
        chrom = "chromosome",
        start = "segment_start",
        end = "segment_end") %>%
      dplyr::mutate(SEGMENT_LENGTH_MB =
                      round((as.numeric((.data$segment_end - .data$segment_start) /
                                          MEGABASE)),
                            digits = 5)) %>%
      dplyr::rename(SEGMENT = .data$SEGMENT_LINK, LOG_R = .data$LogR)

    ## MAKE SIMPLE SEGMENTS DATA FRAME FOR FILTERING IN REPORT
    cna_segments <- cna_df %>%
      dplyr::select(.data$SEGMENT, .data$SEGMENT_LENGTH_MB,
                    .data$CYTOBAND, .data$LOG_R, .data$EVENT_TYPE) %>%
      dplyr::distinct()

    #### FIND AND APPEND GENCODE TRANSCRIPTS THAT OVERLAP
    cna_transcript_df <-
      pcgrr::get_cna_overlapping_transcripts(cna_df, pcgr_data = pcgr_data)

    #### GENERATE DATAFRAME OF UNIQUE TRANSCRIPT-CNA SEGMENTS FOR OUTPUT TSV
    cna_transcript_df_print <- cna_transcript_df %>%
      dplyr::select(.data$chrom, .data$segment_start, .data$segment_end,
                    .data$SEGMENT_ID, .data$SEGMENT_LENGTH_MB,
                    .data$EVENT_TYPE, .data$CYTOBAND, .data$LOG_R, .data$SAMPLE_ID, .data$ensembl_gene_id,
                    .data$symbol, .data$ensembl_transcript_id, .data$transcript_start,
                    .data$transcript_end, .data$transcript_overlap_percent,
                    .data$name, .data$biotype,
                    .data$tumor_suppressor, .data$oncogene,
                    .data$intogen_drivers, .data$chembl_compound_id,
                    .data$gencode_tag, .data$gencode_release) %>%
      magrittr::set_colnames(tolower(names(.)))

    avg_transcript_overlap <- as.data.frame(
      cna_transcript_df %>%
        dplyr::filter(.data$biotype == "protein_coding") %>%
        dplyr::group_by(.data$SEGMENT_ID, .data$symbol) %>%
        dplyr::summarise(
          MEAN_TRANSCRIPT_CNA_OVERLAP = mean(.data$transcript_overlap_percent),
          TRANSCRIPTS = paste0(.data$ensembl_transcript_id, collapse = ", ")) %>%
        dplyr::rename(SYMBOL = .data$symbol) %>%
        dplyr::mutate(
          MEAN_TRANSCRIPT_CNA_OVERLAP =
            round(.data$MEAN_TRANSCRIPT_CNA_OVERLAP, digits = 2))
    )

    cna_transcript_df <- dplyr::select(cna_transcript_df, -.data$ensembl_transcript_id) %>%
      dplyr::filter(.data$biotype == "protein_coding") %>%
      dplyr::distinct() %>%
      dplyr::mutate(VAR_ID = as.character(rep(1:nrow(.)))) %>%
      magrittr::set_colnames(toupper(names(.))) %>%
      dplyr::select(.data$VAR_ID, .data$SEGMENT_ID, .data$SYMBOL, .data$ONCOGENE,
                    .data$ONCOGENE_EVIDENCE, .data$TUMOR_SUPPRESSOR,
                    .data$TUMOR_SUPPRESSOR_EVIDENCE, .data$CANCERGENE_SUPPORT,
                    .data$ENTREZGENE, .data$CHROM, .data$NAME, .data$EVENT_TYPE,
                    .data$SEGMENT_LENGTH_MB, .data$SEGMENT,
                    .data$TRANSCRIPT_OVERLAP_PERCENT, .data$LOG_R) %>%
      dplyr::mutate(ENTREZ_ID = as.character(.data$ENTREZGENE)) %>%
      dplyr::rename(GENENAME = .data$NAME,
                    TRANSCRIPT_OVERLAP = .data$TRANSCRIPT_OVERLAP_PERCENT,
                    CHROMOSOME = .data$CHROM) %>%
      dplyr::left_join(pcgr_data[["kegg"]][["pathway_links"]],
                       by = c("ENTREZ_ID" = "gene_id")) %>%
      dplyr::rename(KEGG_PATHWAY = .data$kegg_pathway_urls)

    ## Get gene annotation links
    entrezgene_annotation_links <-
      pcgrr::generate_annotation_link(
        cna_transcript_df,
        vardb = "GENE_NAME",
        group_by_var = "VAR_ID",
        link_key_var = "ENTREZ_ID",
        link_display_var = "GENENAME",
        url_prefix = "http://www.ncbi.nlm.nih.gov/gene/")

    cna_transcript_df <- cna_transcript_df %>%
      dplyr::left_join(dplyr::rename(entrezgene_annotation_links,
                                     GENE_NAME = .data$link),
                       by = c("VAR_ID")) %>%
      dplyr::select(.data$SEGMENT_ID, .data$CHROMOSOME, .data$SYMBOL, .data$GENE_NAME, .data$KEGG_PATHWAY,
                    .data$TUMOR_SUPPRESSOR, .data$TUMOR_SUPPRESSOR_EVIDENCE, .data$ONCOGENE,
                    .data$ONCOGENE_EVIDENCE, .data$CANCERGENE_SUPPORT, .data$SEGMENT_LENGTH_MB,
                    .data$SEGMENT, .data$EVENT_TYPE, .data$LOG_R) %>%
      dplyr::distinct() %>%
      dplyr::left_join(avg_transcript_overlap, by = c("SEGMENT_ID", "SYMBOL"))


    n_cna_loss <- dplyr::filter(cna_segments, .data$LOG_R <= log_r_homdel) %>%
      nrow()
    n_cna_gain <- dplyr::filter(cna_segments, .data$LOG_R >= log_r_gain) %>%
      nrow()
    cna_segments_filtered <- cna_segments %>%
      dplyr::filter(.data$LOG_R >= log_r_gain | .data$LOG_R <= log_r_homdel) %>%
      dplyr::arrange(dplyr::desc(.data$LOG_R))
    log4r_info(
      paste0("Detected ", nrow(cna_segments_filtered),
             " segments subject to amplification/deletion (",
             n_cna_loss, " deletions, ", n_cna_gain,
             " gains according to user-defined log(2) ratio thresholds)"))


    ## Get aberration sets related to tumor suppressor genes
    ## /oncogenes/drug targets
    onco_ts_sets <-
      pcgrr::get_oncogene_tsgene_target_sets(cna_transcript_df,
                                             log_r_homdel = log_r_homdel,
                                             log_r_gain = log_r_gain,
                                             tumor_type = tumor_type,
                                             pcgr_data = pcgr_data)

    ## load all clinical evidence items ()
    eitems_any_tt <- pcgrr::load_eitems(
      eitems_raw = pcgr_data$biomarkers,
      alteration_type = "CNA",
      ontology =
        pcgr_data$phenotype_ontology$oncotree,
      origin = "Somatic",
      tumor_type_specificity = "any")



    ## Get all clinical evidence items that are related to
    ## tumor suppressor genes/oncogenes/drug targets (NOT tumor-type specific)
    biomarker_hits_cna_any <-
      pcgrr::get_clin_assocs_cna(onco_ts_sets,
                                 annotation_tags = pcgr_data$annotation_tags,
                                 eitems = eitems_any_tt)

    pcg_report_cna[["clin_eitem"]][["any_ttype"]] <-
      biomarker_hits_cna_any[["clin_eitem"]]
    pcg_report_cna[["variant_set"]][["tier2"]] <-
      biomarker_hits_cna_any$variant_set

    ## Get all clinical evidence items that
    ## overlap query set (if tumor type is specified)
    if (tumor_type != "Cancer, NOS") {

      ## load tumor-type specific evidence items ()
      eitems_specific_tt <- pcgrr::load_eitems(
        eitems_raw = pcgr_data$biomarkers,
        alteration_type = "CNA",
        ontology =
          pcgr_data$phenotype_ontology$oncotree,
        origin = "Somatic",
        tumor_type_specificity = "specific",
        tumor_type = tumor_type)

      biomarker_hits_cna_specific <-
        pcgrr::get_clin_assocs_cna(onco_ts_sets,
                                   annotation_tags = pcgr_data$annotation_tags,
                                   eitems = eitems_specific_tt)

      ## Assign putative TIER 1 variant set
      pcg_report_cna[["clin_eitem"]][["specific_ttype"]] <-
        biomarker_hits_cna_specific$clin_eitem
      pcg_report_cna[["variant_set"]][["tier1"]] <-
        biomarker_hits_cna_specific$variant_set
    }

    pcg_report_cna[["eval"]] <- T
    pcg_report_cna[["variant_set"]][["tsv"]] <-
      cna_transcript_df_print
    pcg_report_cna[["v_stat"]][["n_cna_gain"]] <-
      n_cna_gain
    pcg_report_cna[["v_stat"]][["n_cna_loss"]] <-
      n_cna_loss
    pcg_report_cna[["disp"]][["segment"]] <-
      cna_segments_filtered
    pcg_report_cna[["disp"]][["oncogene_gain"]] <-
      onco_ts_sets[["oncogene_gain"]]
    pcg_report_cna[["disp"]][["tsgene_loss"]] <-
      onco_ts_sets[["tsgene_loss"]]
    pcg_report_cna[["disp"]][["other_target"]] <-
      onco_ts_sets[["other_target"]]


    pcg_report_cna <-
      pcgrr::assign_tier1_tier2_acmg_cna(pcg_report_cna)

    return(pcg_report_cna)
  }



#' Function that generates dense and tiered annotated variant datasets
#' @param variant_set List with tiered variants
#' @param config PCGR configuration settings
#' @param annotation_tags List with display columns
#' @param sample_name Sample identifier
#'
#' @return tsv_variants data frame with tier-annotated list of
#' variants for tab-separated output
#'
#' @export
generate_tier_tsv <- function(variant_set,
                              config,
                              annotation_tags,
                              sample_name = "test") {

  tags <- NULL
  if (!is.null(config[["preserved_info_tags"]])) {
    if (config[["preserved_info_tags"]] != "None") {
      tags <-
        stringr::str_split(
          config[["preserved_info_tags"]], pattern = ",")[[1]]
    }
  }
  log4r_info(paste0(
    "Generating tiered set of result variants for output",
    " in tab-separated values (TSV) file"))
  tsv_variants <- NULL
  for (tier in c("tier1", "tier2", "tier3", "tier4", "noncoding")) {
    if (nrow(variant_set[[tier]]) > 0) {
      tierset <- variant_set[[tier]]
      tierset$VCF_SAMPLE_ID <- sample_name
      tsv_columns <- annotation_tags[["tsv"]]
      if (!is.null(tags)) {
        for (t in tags) {
          t <- stringr::str_trim(t)
          if (t %in% colnames(tierset)) {
            tsv_columns <- c(tsv_columns, t)
          }
        }
      }

      if (tier == "tier1") {
        tierset$TIER_DESCRIPTION <- "Variants of strong clinical significance"
        tierset$TIER <- "TIER 1"
      }
      if (tier == "tier2") {
        tierset$TIER_DESCRIPTION <-
          "Variants of potential clinical significance"
        tierset$TIER <- "TIER 2"
      }
      if (tier == "tier3") {
        tierset$TIER_DESCRIPTION <- "Variants of uncertain significance"
        tierset$TIER <- "TIER 3"
      }
      if (tier == "tier4") {
        tierset$TIER_DESCRIPTION <- "Other coding mutation"
        tierset$TIER <- "TIER 4"
      }
      if (tier == "noncoding") {
        tierset$TIER_DESCRIPTION <- "Noncoding mutation"
        tierset$TIER <- "NONCODING"
      }
      tierset <- tierset %>%
        dplyr::select(dplyr::one_of(tsv_columns)) %>%
        dplyr::distinct()

      tsv_variants <- dplyr::bind_rows(tsv_variants, tierset)
    }
  }
  tsv_variants$GENE_NAME <-
    unlist(lapply(stringr::str_match_all(tsv_variants$GENE_NAME, ">.+<"),
                  paste, collapse = ","))
  tsv_variants$GENE_NAME <-
    stringr::str_replace_all(tsv_variants$GENE_NAME, ">|<", "")
  tsv_variants$CLINVAR <-
    unlist(lapply(stringr::str_match_all(tsv_variants$CLINVAR, ">.+<"),
                  paste, collapse = ","))
  tsv_variants$CLINVAR <-
    stringr::str_replace_all(tsv_variants$CLINVAR, ">|<", "")
  tsv_variants$PROTEIN_DOMAIN <-
    unlist(lapply(stringr::str_match_all(tsv_variants$PROTEIN_DOMAIN, ">.+<"),
                  paste, collapse = ","))
  tsv_variants$PROTEIN_DOMAIN <-
    stringr::str_replace_all(tsv_variants$PROTEIN_DOMAIN, ">|<", "")
  tsv_variants$TCGA_FREQUENCY <-
    stringr::str_replace_all(
      tsv_variants$TCGA_FREQUENCY,
      "<a href='https://portal.gdc.cancer.gov/projects/TCGA-[A-Z]{1,}' target=\"_blank\">|</a>", "")
  tsv_variants <- tsv_variants %>% dplyr::distinct()

  return(tsv_variants)
}


#' Function that writes contents of PCGR object to various output formats
#' (Rmarkdown/flexdashboard HTML reports, JSON, tab-separated etc)
#'
#' @param report List object with all report data (PCGR/CPSR)
#' @param config Configuration object with all key configurations
#' (directories, sample names, genome assembly etc)
#' @param tier_model type of tier model
#' @param output_format contents/file format of output
#' (html/json/tsv/cna_tsv etc)
#' @param flexdb logical indicating if HTML output should be dashboard

#' @export
write_report_output <- function(report,
                                config,
                                tier_model = "pcgr_acmg",
                                output_format = "html",
                                flexdb = F) {

  project_directory <- config[['required_args']][['output_dir']]
  sample_name <- config[['required_args']][['sample_name']]
  genome_assembly <- config[['required_args']][['genome_assembly']]

  sample_fname_pattern <-
    paste(sample_name, tier_model, genome_assembly, sep = ".")

  fnames <- list()
  fnames[["snv_tsv_unfiltered"]] <-
    file.path(project_directory,
              paste0(sample_fname_pattern,
                     ".snvs_indels.unfiltered.tsv"))
  fnames[["msigs_tsv"]] <-
    file.path(project_directory,
              paste0(sample_fname_pattern,
                     ".mutational_signatures.tsv"))
  fnames[["snv_tsv"]] <-
    file.path(project_directory,
              paste0(sample_fname_pattern,
                     ".snvs_indels.tiers.tsv"))
  fnames[["xlsx"]] <-
    file.path(project_directory,
              paste0(sample_fname_pattern,
                     ".snvs_indels.tiers.xlsx"))
  fnames[["cna_tsv"]] <-
    file.path(project_directory,
              paste0(sample_fname_pattern,
                     ".cna_segments.tsv"))
  fnames[["maf"]] <-
    file.path(project_directory,
              paste0(sample_fname_pattern, ".maf"))
  fnames[["maf_tmp"]] <-
    file.path(project_directory,
              paste0(sample_fname_pattern, ".tmp.maf"))
  fnames[["json"]] <-
    file.path(project_directory,
              paste0(sample_fname_pattern, ".json"))
  fnames[["html"]] <-
    file.path(project_directory,
              paste0(sample_fname_pattern, ".html"))
  if (flexdb == T) {
    fnames[["html"]] <-
      file.path(project_directory,
                paste0(sample_fname_pattern,
                       ".flexdb.html"))
  }



  ## Set to CPSR/germline settings as default
  sequencing_design <- "Germline"
  disclaimer <- system.file("templates",
                            "disclaimer_predisposition.md",
                            package = "pcgrr")
  markdown_input <- system.file("templates",
                                "cpsr_rmarkdown_report.Rmd",
                                package = "pcgrr")
  css_fname <- system.file("templates", "cpsr.css",
                           package = "pcgrr")
  report_theme <-
    report[["metadata"]][["config"]][["visual"]][["report_theme"]]

  ## Somatic/tumor report settings
  if(tier_model == "pcgr_acmg"){

    disclaimer <- system.file("templates",
                              "disclaimer.md",
                              package = "pcgrr")
    assay_props <-
      report[["metadata"]][["config"]][["assay_props"]]
    sequencing_assay <-
      assay_props[["type"]]

    ## Flexdashboard layout
    sequencing_design <- "Tumor-Control"
    markdown_input <- system.file("templates",
                                  "pcgr_flexdb_report.Rmd",
                                  package = "pcgrr")
    css_fname <- system.file("templates", "pcgr_flexdb_tumor_control.css",
                             package = "pcgrr")

    ## Rmarkdown layout
    if(flexdb == F){
      markdown_input <-
        system.file("templates",
                    "pcgr_rmarkdown_report.Rmd",
                    package = "pcgrr")
      css_fname <-
        system.file("templates", "pcgr_rmarkdown_tumor_control.css",
                    package = "pcgrr")
    }

    ## Tumor-only settings (CSS)
    if (assay_props[["vcf_tumor_only"]] == T) {
      sequencing_design <- "Tumor-Only"
      css_fname <-
        system.file("templates", "pcgr_flexdb_tumor_only.css",
                    package = "pcgrr")

      if (flexdb == F){
        css_fname <-
          system.file("templates", "pcgr_rmarkdown_tumor_only.css",
                      package = "pcgrr")
      }
    }
  }

  if (output_format == "html") {

    if (flexdb == T & tier_model == "pcgr_acmg") {
      log4r_info("------")
      log4r_info(
        "Writing HTML file (.html) with report contents - flexdashboard")
      navbar_items <- list()
      navbar_items[[1]] <-
        list("title" = paste0(
          "<b>", sample_name, "</b> | <i>",
          report[["metadata"]][["config"]][["t_props"]][["tumor_type"]],
          "</i> | ", sequencing_design, " | ", sequencing_assay),
          href = "", target = "_blank", align = "right")
      navbar_items[[2]] <-
        list("icon" = "fa-github",
             href = "https://github.com/sigven/pcgr", target = "_blank",
             align = "right")

      rmarkdown::render(
        markdown_input,
        output_format =
          flexdashboard::flex_dashboard(
            orientation = "rows",
            theme = "cosmo",
            css = css_fname,
            navbar = navbar_items),
        output_file = fnames[["html"]],
        output_dir = project_directory,
        clean = T,
        intermediates_dir = project_directory,
        quiet = T)
    }else{

      toc_float <-
        list(collapsed = TRUE,
             smooth_scroll = TRUE,
             print = TRUE)
      toc_depth <- 3

      ## Ignore collapsing menu for CPSR
      if(tier_model == 'cpsr'){
        toc_float <-
          list(collapsed = FALSE,
               smooth_scroll = FALSE,
               print = TRUE)
      }

      ## If nonfloating TOC is chosen (PCGR/CPSR), set toc_float to FALSE
      nonfloating_toc <-
        report[["metadata"]][["config"]][["visual"]][["nonfloating_toc"]]
      if(nonfloating_toc == T){
        toc_float <- F
      }

      log4r_info("------")
      log4r_info(paste0(
        "Writing HTML file (.html) with report contents - rmarkdown (theme = '",
        report_theme,"')"))
      rmarkdown::render(
        markdown_input,
        output_format =
          rmarkdown::html_document(
            theme = report_theme,
            fig_width = 5,
            fig_height = 4,
            toc = T,
            toc_depth = toc_depth,
            toc_float = toc_float,
            number_sections = F,
            css = css_fname,
            includes =
              rmarkdown::includes(
                after_body = disclaimer)),
        output_file = fnames[["html"]],
        output_dir = project_directory,
        clean = T,
        intermediates_dir = project_directory,
        quiet = T)
    }
  }
  if (output_format == "json") {
    if (!is.null(report[["cna_plot"]][["png"]])) {
      report[["cna_plot"]][["png"]] <- NULL
    }
    if (!is.null(report[["tmb"]][["tcga_tmb"]])) {
      report[["tmb"]][["tcga_tmb"]] <- NULL
    }
    log4r_info("------")
    log4r_info("Writing JSON file (.json) with report contents")
    pcgr_json <- jsonlite::toJSON(report, pretty = T, na = "string",
                                  null = "null", force = T)
    write(pcgr_json, fnames[["json"]])
    gzip_command <- paste0("gzip -f ", fnames[["json"]])
    system(gzip_command, intern = F)
  }

  if (output_format == "snv_tsv" | output_format == "snv_tsv_unfiltered") {
    output_format_slim <- stringr::str_replace(output_format, "snv_", "")
    if (NROW(
      report[["content"]][["snv_indel"]][["variant_set"]][[output_format_slim]]) > 0) {
      log4r_info("------")
      if (tier_model == "pcgr_acmg") {
        log4r_info(
          paste0("Writing SNV/InDel tab-separated output file with ",
                 "PCGR annotations - ('",
                 output_format_slim, "')"))
      }
      if (tier_model == "cpsr") {
        log4r_info(
          paste0("Writing SNV/InDel tab-separated output file ",
                 "with CPSR annotations - ('",
                 output_format_slim, "')"))
      }
      utils::write.table(
        report[["content"]][["snv_indel"]][["variant_set"]][[output_format_slim]],
        file = fnames[[output_format]], sep = "\t", col.names = T,
        row.names = F, quote = F)

      # if(tier_model == "pcgr_acmg"){
      #   log4r_info(
      #     paste0("Writing SNV/InDel Excel output file with ",
      #            "PCGR annotations"))
      #   workbook <- openxlsx::createWorkbook()
      #   openxlsx::addWorksheet(workbook,
      #                          sheetName = "SNV_INDELS")
      #
      #   ## set automatic column widths
      #   openxlsx::setColWidths(
      #     workbook,
      #     sheet = "SNV_INDELS",
      #     cols = 1:ncol(report[["content"]][["snv_indel"]][["variant_set"]][[output_format_slim]]),
      #     widths = "auto")
      #
      #   ## write with default Excel Table style
      #   openxlsx::writeDataTable(
      #     workbook,
      #     sheet = "SNV_INDELS",
      #     x = report[["content"]][["snv_indel"]][["variant_set"]][[output_format_slim]],
      #     startRow = 1,
      #     startCol = 1,
      #     colNames = TRUE,
      #     tableStyle = "TableStyleMedium15")
      #
      #   openxlsx::saveWorkbook(
      #     workbook,
      #     fnames[['excel']],
      #     overwrite = TRUE)
      # }
    }
  }

  if (output_format == "msigs_tsv") {
    if (
      NROW(report[["content"]][["m_signature_mp"]][["result"]][["tsv"]]) > 0) {
      log4r_info("------")
      log4r_info(paste0(
        "Writing tab-separated output file with details ",
        "of contributing mutational signatures - ('tsv')"))
      utils::write.table(report[["content"]][["m_signature_mp"]][["result"]][["tsv"]],
                  file = fnames[[output_format]], sep = "\t", col.names = T,
                  row.names = F, quote = F)
    }
  }

  if (output_format == "cna_tsv") {
    if (NROW(report[["content"]][["cna"]][["variant_set"]][["tsv"]]) > 0) {
      log4r_info("------")
      log4r_info(
        "Writing CNA tab-separated output file with PCGR annotations (.tsv.gz)")
      utils::write.table(report[["content"]][["cna"]][["variant_set"]][["tsv"]],
                  file = fnames[["cna_tsv"]], sep = "\t", col.names = T,
                  row.names = F, quote = F)
      gzip_command <- paste0("gzip -f ", fnames[["cna_tsv"]])
      system(gzip_command, intern = F)
    }
  }

}

