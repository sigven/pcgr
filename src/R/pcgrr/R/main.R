#' Function that generates all contents of the cancer genome report (PCGR)
#'
#' @param project_directory name of project directory
#' @param query_vcf2tsv name of gzipped TSV file (vcf2tsv) with annotated
#' query SNVs/InDels
#' @param pcgr_data List object with multiple PCGR data bundle annotations
#' @param pcgr_version PCGR software version
#' @param config Object with PCGR configuration parameters
#' @param sample_name sample identifier
#' @param cna_segments_tsv name of CNA segments file (tab-separated values)
#' @param cna_plot Path to PNG image with CNA plot
#' @param tier_model Variant tier model (pcgr_acmg/cpsr)
#'

generate_pcgr_report <- function(project_directory,
                                 query_vcf2tsv,
                                 pcgr_data = NULL,
                                 pcgr_version = "dev",
                                 config = NULL,
                                 sample_name = "SampleX",
                                 cna_segments_tsv = NULL,
                                 cna_plot = NULL,
                                 tier_model = "pcgr_acmg") {


  rlogging::message(paste0("Initializing PCGR report - sample ", sample_name))
  #rlogging::message(head(pcgr_data$tcga$tmb))

  ## Assert argument values
  invisible(assertthat::assert_that(
    !is.null(config) & is.list(config),
    msg = "Argument 'config' must be a non-NULL list object"))
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
  pcg_report <- pcgrr::init_report(
    config,
    sample_name,
    class = NULL,
    pcgr_data = pcgr_data,
    pcgr_version = pcgr_version)

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
      pcgrr::generate_report_data_trials(pcgr_data = pcgr_data,
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
          pcgr_data, sample_name,
          config, callset = "germline-filtered callset",
          tier_model = tier_model)

      pcg_report_tumor_only[["upset_data"]] <-
        pcgrr::make_upset_plot_data(
          pcg_report_tumor_only$variant_set$unfiltered, config)
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
          sample_calls, pcgr_data, sample_name,
          config, tier_model = tier_model)

      ## Update genome report
      pcg_report <- pcgrr::update_report(
        pcg_report, pcg_report_snv_indel,
        a_elem = "snv_indel")
    }
    rm(sample_calls)

    ## Estimate contribution of mutational signatures
    if (pcg_report[["metadata"]][["config"]][["msigs"]][["run"]] == T) {

      pcgrr::write_processed_vcf(
        calls = pcg_report$content$snv_indel$variant_set$tsv,
        sample_name = sample_name,
        output_directory = project_directory,
        vcf_fname = fnames[["vcf_mp"]])

      pcg_report_signatures <-
        pcgrr::generate_report_data_signatures_mp(
          vcf_fname = paste0(fnames[["vcf_mp"]], ".gz"), pcgr_data,
          sample_name, config,
          type_specific =
            !config[["msigs"]][["all_reference_signatures"]])

      ## Update genome report with signature info
      pcg_report <- pcgrr::update_report(
        pcg_report, pcg_report_signatures,
        a_elem = "m_signature_mp")

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
  pcg_report[["content"]][["cna"]][["variant_set"]][["cna_print"]] <- NULL
  pcg_report[["metadata"]][["phenotype_ontology"]] <- list()

  if (!is.null(cna_plot) && cna_plot != "None") {
    pcg_report[["content"]][["cna_plot"]][["png"]] <- cna_plot
    pcg_report[["content"]][["cna_plot"]][["eval"]] <- TRUE
  }
  return(pcg_report)
}


#' Function that generates tiered variant sets for SNVs/InDels
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param biomarker_mapping_stringency quality level for biomarkers
#' @param callset type of calls
#' @param tier_model tier model (pcgr_acmg)
#'
#' @return pcg_report_data data frame with all report elements
#'
generate_report_data_snv_indel <- function(sample_calls,
                                           pcgr_data,
                                           sample_name,
                                           config,
                                           callset = "somatic calls",
                                           biomarker_mapping_stringency = 1,
                                           tier_model = "pcgr_acmg") {

  rlogging::message("------")
  rlogging::message(
    paste0("Generating data for tiered cancer genome report - ",
           callset, " tier model '", tier_model, "'"))

  pcg_report_snv_indel <- pcgrr::init_report(config,
                                             sample_name,
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
  rlogging::message(
    paste0("Number of protein-coding variants: ",
           pcg_report_snv_indel[["v_stat"]][["n_coding"]]))

  ## Add custom annotation tags to lists of tags to display
  if (!is.null(config[["custom_tags"]])) {
    if (config[["custom_tags"]][["custom_tags"]] != "") {
      tags <- stringr::str_split(config[["custom_tags"]][["custom_tags"]],
                                 pattern = ",")[[1]]
      for (t in tags) {
        t <- stringr::str_trim(t)
        if (t %in% colnames(sample_calls)) {
          pcgr_data[["annotation_tags"]][["all"]] <-
            c(pcgr_data[["annotation_tags"]][["all"]], t)
        }
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
        annotation_tags = pcgr_data$annotation_tags,
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
          annotation_tags = pcgr_data$annotation_tags,
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
                    dplyr::one_of(pcgr_data[["annotation_tags"]][["all"]])) %>%
      dplyr::filter(CODING_STATUS == "coding") %>%
      dplyr::filter(ONCOGENE == TRUE | TUMOR_SUPPRESSOR == TRUE)
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
        dplyr::arrange(desc(OPENTARGETS_RANK))
      tier123 <- tier12 %>%
        dplyr::bind_rows(
          dplyr::select(pcg_report_snv_indel[["variant_set"]][["tier3"]],
                        GENOMIC_CHANGE)) %>%
        dplyr::distinct()
      pcg_report_snv_indel[["disp"]][["tier3"]][["proto_oncogene"]] <-
        dplyr::select(
          pcg_report_snv_indel[["variant_set"]][["tier3"]],
          dplyr::one_of(pcgr_data[["annotation_tags"]][["tier3_display"]])) %>%
        dplyr::filter(ONCOGENE == TRUE &
                        (is.na(TUMOR_SUPPRESSOR) |
                           TUMOR_SUPPRESSOR == FALSE))
      pcg_report_snv_indel[["disp"]][["tier3"]][["tumor_suppressor"]] <-
        dplyr::select(
          pcg_report_snv_indel[["variant_set"]][["tier3"]],
          dplyr::one_of(pcgr_data[["annotation_tags"]][["tier3_display"]])) %>%
        dplyr::filter(TUMOR_SUPPRESSOR == TRUE)
    }

    ## Determine TIER 4: Other coding mutations
    pcg_report_snv_indel[["variant_set"]][["tier4"]] <-
      dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]],
                    dplyr::one_of(pcgr_data[["annotation_tags"]][["all"]])) %>%
      dplyr::filter(CODING_STATUS == "coding")
    if (nrow(tier123) > 0 &
        nrow(pcg_report_snv_indel[["variant_set"]][["tier4"]]) > 0) {
      pcg_report_snv_indel[["variant_set"]][["tier4"]] <-
        dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["tier4"]],
                         tier123, by = c("GENOMIC_CHANGE"))
    }
    if (nrow(pcg_report_snv_indel[["variant_set"]][["tier4"]]) > 0) {
      pcg_report_snv_indel[["variant_set"]][["tier4"]] <-
        pcg_report_snv_indel[["variant_set"]][["tier4"]] %>%
        dplyr::arrange(desc(OPENTARGETS_RANK))
      pcg_report_snv_indel[["disp"]][["tier4"]] <-
        dplyr::select(
          pcg_report_snv_indel[["variant_set"]][["tier4"]],
          dplyr::one_of(pcgr_data[["annotation_tags"]][["tier4_display"]]))
    }

    ## Determine non-coding mutation set
    pcg_report_snv_indel[["variant_set"]][["noncoding"]] <-
      dplyr::select(pcg_report_snv_indel[["variant_set"]][["all"]],
                    dplyr::one_of(pcgr_data[["annotation_tags"]][["all"]])) %>%
      dplyr::filter(CODING_STATUS == "noncoding")
    if (nrow(pcg_report_snv_indel[["variant_set"]][["noncoding"]]) > 0) {
      if (nrow(tier123) > 0) {
        pcg_report_snv_indel[["variant_set"]][["noncoding"]] <-
          dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["noncoding"]],
                           tier123,
                           by = c("GENOMIC_CHANGE"))
      }
      pcg_report_snv_indel[["variant_set"]][["noncoding"]] <-
        pcg_report_snv_indel[["variant_set"]][["noncoding"]] %>%
        dplyr::arrange(desc(OPENTARGETS_RANK))
      pcg_report_snv_indel[["disp"]][["noncoding"]] <-
        dplyr::select(
          pcg_report_snv_indel[["variant_set"]][["noncoding"]],
          dplyr::one_of(pcgr_data[["annotation_tags"]][["tier5_display"]]))
    }

    ## Make TSV content with variant set
    pcg_report_snv_indel[["v_stat"]][["n_noncoding"]] <-
      pcg_report_snv_indel[["variant_set"]][["noncoding"]] %>% nrow()
    pcg_report_snv_indel[["variant_set"]][["tsv"]] <-
      pcgrr::generate_tier_tsv(
        pcg_report_snv_indel[["variant_set"]],
        pcgr_data = pcgr_data, config, sample_name = sample_name)

  }

  return(pcg_report_snv_indel)

}


#' Function that generates germline-filtered callset and PCGR
#' report statistics for a given tumor-only callsets
#'
#' @param unfiltered_sample_calls variant calls
#' @param pcgr_version PCGR software version
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param genome_seq BSgenome object
#' @param genome_assembly human genome assembly version
#'
generate_report_data_tumor_only <- function(unfiltered_sample_calls,
                                            sample_name, pcgr_config) {

  sample_calls <- unfiltered_sample_calls
  gline_filter_stats <- list()
  for (m in c("remain_post_onekg", "remain_post_gnomad", "remain_post_clinvar",
              "remain_post_dbsnp", "remain_post_pon",
              "remain_post_nonexonic", "remain_post_hom", "remain_post_het")) {
    gline_filter_stats[m] <- 0
  }

  ## initiate report
  pcg_report_to <- pcgrr::init_report(pcgr_config,
                                      sample_name, class = "tumor_only")

  ## assign evidence tags for germline/somatic state of variants,
  ## partially based on user-defined criteria
  ## (population allele frequency thresholds)
  vcalls <- pcgrr::assign_somatic_germline_evidence(sample_calls, pcgr_config)

  ## assign somatic classification based on accumulation
  ## of evidence tags and user-defined options
  vcalls <- pcgrr::assign_somatic_classification(vcalls, pcgr_config)

  ## Assign statistics to successive filtering levels for
  ## different evidence criteria
  ## excluded germline calls found in 1000 Genomes Project
  gline_filter_stats[["remain_post_onekg"]] <-
    nrow(vcalls) -
    nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_1KG", ])
  rlogging::message(paste0("Excluding coinciding germline variants in ",
                           "1000 Genomes Project populations"))
  rlogging::message(paste0("Total sample calls remaining: ",
                           gline_filter_stats[["remain_post_onekg"]]))

  ## excluded germline calls found in gnomAD
  gline_filter_stats[["remain_post_gnomad"]] <-
    gline_filter_stats[["remain_post_onekg"]] -
    nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_GNOMAD", ])
  rlogging::message(
    paste0("Excluding coinciding germline variants in any ",
           "population in the genome aggregation database (gnomAD)"))
  rlogging::message(paste0("Total sample calls remaining: ",
                           gline_filter_stats[["remain_post_gnomad"]]))

  ## excluded germline calls found in ClinVar
  gline_filter_stats[["remain_post_clinvar"]] <-
    gline_filter_stats[["remain_post_gnomad"]] -
    nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_CLINVAR", ])
  rlogging::message(paste0("Excluding coinciding germline variants in ClinVar"))
  rlogging::message(paste0("Total sample calls remaining: ",
                           gline_filter_stats[["remain_post_clinvar"]]))


  ## excluded germline calls found in panel of normals (if provided)
  gline_filter_stats[["remain_post_pon"]] <-
    gline_filter_stats[["remain_post_clinvar"]]
  if (pcgr_config[["tumor_only"]][["exclude_pon"]] == TRUE) {
    gline_filter_stats[["remain_post_pon"]] <-
      gline_filter_stats[["remain_post_pon"]] -
      nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_PON", ])
    rlogging::message(
      "Excluding putative germline variants found in calls ",
      "from panel-of-normals (PON)")
    rlogging::message(
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
    rlogging::message(
      paste0("Excluding likely homozygous germline variants found ",
             "as variants with 100% allelic fraction"))
    rlogging::message(paste0("Total sample calls remaining: ",
                             gline_filter_stats[["remain_post_hom"]]))
  }

  ## excluded germline calls found as likely heterozygous germline variants
  gline_filter_stats[["remain_post_het"]] <-
    gline_filter_stats[["remain_post_hom"]]
  if (pcgr_config[["tumor_only"]][["exclude_likely_het_germline"]] == TRUE) {
    gline_filter_stats[["remain_post_het"]] <-
      gline_filter_stats[["remain_post_het"]] -
      nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_HETEROZYGOUS", ])
    rlogging::message(paste0(
      "Excluding likely heterozygous germline variants found as variants ",
      "with 40-60% allelic fraction and recorded in gnomAD + dbSNP"))
    rlogging::message(paste0("Total sample calls remaining: ",
                             gline_filter_stats[["remain_post_het"]]))
  }

  ## excluded calls with dbSNP germline status (if set in config)
  gline_filter_stats[["remain_post_dbsnp"]] <-
    gline_filter_stats[["remain_post_het"]]
  if (pcgr_config[["tumor_only"]][["exclude_dbsnp_nonsomatic"]] == TRUE) {

    rlogging::message(
      paste0("Excluding non-somatically associated dbSNP variants ",
             "(dbSNP - not recorded as somatic in DoCM/ClinVar",
             "and not registered in COSMIC or found in TCGA"))

    gline_filter_stats[["remain_post_dbsnp"]] <-
      gline_filter_stats[["remain_post_dbsnp"]] -
      nrow(vcalls[vcalls$SOMATIC_CLASSIFICATION == "GERMLINE_DBSNP", ])
    rlogging::message(paste0("Total sample calls remaining: ",
                             gline_filter_stats[["remain_post_dbsnp"]]))
  }

  unfiltered_sample_calls <- vcalls
  vcalls <- vcalls %>%
    dplyr::filter(SOMATIC_CLASSIFICATION == "SOMATIC")

  gline_filter_stats[["remain_post_nonexonic"]] <-
    gline_filter_stats[["remain_post_dbsnp"]]
  if (pcgr_config[["tumor_only"]][["exclude_nonexonic"]] == TRUE) {
    rlogging::message(paste0("Excluding non-exonic variants"))
    vcalls <- dplyr::filter(vcalls, EXONIC_STATUS == "exonic")
    rlogging::message(paste0("Total sample calls remaining: ",
                             nrow(vcalls)))
    gline_filter_stats[["remain_post_nonexonic"]] <- nrow(vcalls)
  }


  pcg_report_to[["eval"]] <- TRUE
  pcg_report_to[["variant_set"]][["unfiltered"]] <- unfiltered_sample_calls
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
    pcg_report_cna <- pcgrr::init_report(pcgr_config,
                                         sample_name, class = "cna")
    log_r_homdel <- pcgr_config[["cna"]][["log_r_homdel"]]
    log_r_gain <- pcgr_config[["cna"]][["log_r_gain"]]
    tumor_type <- pcgr_config[["t_props"]][["tumor_type"]]
    MEGABASE <- 1000000

    rlogging::message("------")
    rlogging::message("Generating report data for copy number segment file ",
                      cna_segments_tsv)

    ## READ INPUT FILE, VALIDATE INPUT CHROMOSOMES AND SEGMENTS, ADD CYTOBAND INFO
    cna_df <- read.table(file = cna_segments_tsv, header = T,
                         stringsAsFactors = F, sep = "\t",
                         comment.char = "", quote = "") %>%
      dplyr::rename(chromosome = Chromosome, LogR = Segment_Mean,
                    segment_start = Start, segment_end = End) %>%
      dplyr::distinct() %>%
      dplyr::select(chromosome, LogR, segment_start, segment_end) %>%
      dplyr::mutate(
        chromosome = stringr::str_replace(chromosome, "^chr", "")) %>%
      pcgrr::get_valid_chromosomes(
        chromosome_column = "chromosome",
        bsg = pcgr_data[["assembly"]][["bsg"]]) %>%
      pcgrr::get_valid_chromosome_segments(
        genome_assembly = pcgr_data[["assembly"]][["grch_name"]],
        bsg = pcgr_data[["assembly"]][["bsg"]]) %>%
      dplyr::filter(!is.na(LogR)) %>%
      dplyr::mutate(LogR = round(as.numeric(LogR), digits = 3)) %>%
      dplyr::mutate(SEGMENT_ID = paste0(chromosome, ":",
                                        segment_start, "-",
                                        segment_end)) %>%
      pcgrr::get_cna_cytoband(pcgr_data = pcgr_data) %>%
      dplyr::mutate(SAMPLE_ID = sample_name) %>%
      pcgrr::append_ucsc_segment_link(
        hgname = pcgr_data[["assembly"]][["hg_name"]],
        chrom = "chromosome",
        start = "segment_start",
        end = "segment_end") %>%
      dplyr::mutate(SEGMENT_LENGTH_MB =
                      round((as.numeric((segment_end - segment_start) /
                                          MEGABASE)),
                            digits = 5)) %>%
      dplyr::rename(SEGMENT = SEGMENT_LINK, LOG_R = LogR)

    ## MAKE SIMPLE SEGMENTS DATA FRAME FOR FILTERING IN REPORT
    cna_segments <- cna_df %>%
      dplyr::select(SEGMENT, SEGMENT_LENGTH_MB,
                    CYTOBAND, LOG_R, EVENT_TYPE) %>%
      dplyr::distinct()

    #### FIND AND APPEND GENCODE TRANSCRIPTS THAT OVERLAP
    cna_transcript_df <-
      pcgrr::get_cna_overlapping_transcripts(cna_df, pcgr_data = pcgr_data)

    #### GENERATE DATAFRAME OF UNIQUE TRANSCRIPT-CNA SEGMENTS FOR OUTPUT TSV
    cna_transcript_df_print <- cna_transcript_df %>%
      dplyr::select(chrom, segment_start, segment_end,
                    SEGMENT_ID, SEGMENT_LENGTH_MB,
                    EVENT_TYPE, CYTOBAND, LOG_R, SAMPLE_ID, ensembl_gene_id,
                    symbol, ensembl_transcript_id, transcript_start,
                    transcript_end, transcript_overlap_percent,
                    name, biotype,
                    tumor_suppressor, oncogene,
                    intogen_drivers, chembl_compound_id,
                    gencode_tag, gencode_release) %>%
      magrittr::set_colnames(tolower(names(.)))

    avg_transcript_overlap <- as.data.frame(
      cna_transcript_df %>%
        dplyr::filter(biotype == "protein_coding") %>%
        dplyr::group_by(SEGMENT_ID, symbol) %>%
        dplyr::summarise(
          MEAN_TRANSCRIPT_CNA_OVERLAP = mean(transcript_overlap_percent),
          TRANSCRIPTS = paste0(ensembl_transcript_id, collapse = ", ")) %>%
        dplyr::rename(SYMBOL = symbol) %>%
        dplyr::mutate(
          MEAN_TRANSCRIPT_CNA_OVERLAP =
            round(MEAN_TRANSCRIPT_CNA_OVERLAP, digits = 2))
    )

    cna_transcript_df <- dplyr::select(cna_transcript_df, -ensembl_transcript_id) %>%
      dplyr::filter(biotype == "protein_coding") %>%
      dplyr::distinct() %>%
      dplyr::mutate(VAR_ID = as.character(rep(1:nrow(.)))) %>%
      magrittr::set_colnames(toupper(names(.))) %>%
      dplyr::select(VAR_ID, SEGMENT_ID, SYMBOL, ONCOGENE,
                    ONCOGENE_EVIDENCE, TUMOR_SUPPRESSOR,
                    TUMOR_SUPPRESSOR_EVIDENCE, CANCERGENE_SUPPORT,
                    ENTREZGENE, CHROM, NAME, EVENT_TYPE,
                    SEGMENT_LENGTH_MB, SEGMENT,
                    TRANSCRIPT_OVERLAP_PERCENT, LOG_R) %>%
      dplyr::mutate(ENTREZ_ID = as.character(ENTREZGENE)) %>%
      dplyr::rename(GENENAME = NAME,
                    TRANSCRIPT_OVERLAP = TRANSCRIPT_OVERLAP_PERCENT,
                    CHROMOSOME = CHROM) %>%
      dplyr::left_join(pcgr_data[["kegg"]][["pathway_links"]],
                       by = c("ENTREZ_ID" = "gene_id")) %>%
      dplyr::rename(KEGG_PATHWAY = kegg_pathway_urls)

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
                                     GENE_NAME = link),
                       by = c("VAR_ID")) %>%
      dplyr::select(SEGMENT_ID, CHROMOSOME, SYMBOL, GENE_NAME, KEGG_PATHWAY,
                    TUMOR_SUPPRESSOR, TUMOR_SUPPRESSOR_EVIDENCE, ONCOGENE,
                    ONCOGENE_EVIDENCE, CANCERGENE_SUPPORT, SEGMENT_LENGTH_MB,
                    SEGMENT, EVENT_TYPE, LOG_R) %>%
      dplyr::distinct() %>%
      dplyr::left_join(avg_transcript_overlap, by = c("SEGMENT_ID", "SYMBOL"))


    n_cna_loss <- dplyr::filter(cna_segments, LOG_R <= log_r_homdel) %>%
      nrow()
    n_cna_gain <- dplyr::filter(cna_segments, LOG_R >= log_r_gain) %>%
      nrow()
    cna_segments_filtered <- cna_segments %>%
      dplyr::filter(LOG_R >= log_r_gain | LOG_R <= log_r_homdel) %>%
      dplyr::arrange(desc(LOG_R))
    rlogging::message(
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
#' @param pcgr_data List of data frames with PCGR data annotations
#' @param sample_name Sample identifier
#'
#' @return tsv_variants data frame with tier-annotated list of
#' variants for tab-separated output
#'
generate_tier_tsv <- function(variant_set, pcgr_data, config,
                              sample_name = "test") {

  tags <- NULL
  if (!is.null(config[["custom_tags"]])) {
    if (config[["custom_tags"]][["custom_tags"]] != "") {
      tags <-
        stringr::str_split(
          config[["custom_tags"]][["custom_tags"]], pattern = ",")[[1]]
    }
  }
  rlogging::message(
    "Generating tiered set of result variants for output",
    " in tab-separated values (TSV) file")
  tsv_variants <- NULL
  for (tier in c("tier1", "tier2", "tier3", "tier4", "noncoding")) {
    if (nrow(variant_set[[tier]]) > 0) {
      tierset <- variant_set[[tier]]
      tierset$VCF_SAMPLE_ID <- sample_name
      tsv_columns <- pcgr_data[["annotation_tags"]][["tsv"]]
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
#' @param project_directory working directory
#' @param report List object with all report data (PCGR/CPSR)
#' @param sample_name sample name
#' @param genome_assembly genome assembly (grch37/grch38)
#' @param tier_model type of tier model
#' @param output_format contents/file format of output
#' (html/json/tsv/cna_tsv etc)
#' @param flexdb logical indicating if HTML output should be dashboard

write_report_output <- function(project_directory,
                                report,
                                sample_name,
                                genome_assembly,
                                tier_model = "pcgr_acmg",
                                output_format = "html",
                                flexdb = F) {

  sample_fname_pattern <-
    paste(sample_name, tier_model, genome_assembly, sep = ".")


  fnames <- list()
  fnames[["snv_tsv_unfiltered"]] <-
    file.path(project_directory, paste0(sample_fname_pattern,
                                        ".snvs_indels.tiers.unfiltered.tsv"))
  fnames[["msigs_tsv"]] <-
    file.path(project_directory, paste0(sample_fname_pattern,
                                        ".mutational_signatures.tsv"))
  fnames[["snv_tsv"]] <-
    file.path(project_directory, paste0(sample_fname_pattern,
                                        ".snvs_indels.tiers.tsv"))
  fnames[["cna_tsv"]] <-
    file.path(project_directory, paste0(sample_fname_pattern,
                                        ".cna_segments.tsv"))
  fnames[["maf"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".maf"))
  fnames[["maf_tmp"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".tmp.maf"))
  fnames[["json"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".json"))
  fnames[["html"]] <-
    file.path(project_directory, paste0(sample_fname_pattern, ".html"))
  if (flexdb == T) {
    fnames[["html"]] <-
      file.path(project_directory, paste0(sample_fname_pattern,
                                          ".flexdb.html"))
  }

  disclaimer <- "disclaimer.md"
  report_theme <-
    report[["metadata"]][["config"]][["visual"]][["report_theme"]]
  assay_props <-
    report[["metadata"]][["config"]][["assay_props"]]


  sequencing_assay <- "Unknown"
  sequencing_design <- "Tumor-Control"
  css_fname <- system.file("templates", "pcgr_flexdb_tumor_control.css",
                           package = "pcgrr")

  if (tier_model == "cpsr") {
    sequencing_design <- "Germline"
    disclaimer <- "disclaimer_predisposition.md"
  }else{
    sequencing_assay <-
      assay_props[["type"]]
    if (assay_props[["vcf_tumor_only"]] == T) {
      sequencing_design <- "Tumor-Only"
      css_fname <- system.file("templates", "pcgr_flexdb_tumor_only.css",
                               package = "pcgrr")
    }
  }

  if (output_format == "html") {

    markdown_input <- system.file("templates",
                                  "pcgr_rmarkdown_report.Rmd",
                                  package = "pcgrr")
    if (flexdb == T) {
      rlogging::message("------")
      rlogging::message(
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

      markdown_input <- system.file("templates",
                                    "pcgr_flexdb_report.Rmd",
                                    package = "pcgrr")
      rmarkdown::render(markdown_input,
                        output_format =
                          flexdashboard::flex_dashboard(orientation = "rows",
                                                        theme = "cosmo",
                                                        css = css_fname,
                                                        navbar = navbar_items),
                        output_file = fnames[["html"]],
                        output_dir = project_directory,
                        clean = T,
                        intermediates_dir = project_directory, quiet = T)
    }else{
      rlogging::message("------")
      rlogging::message(
        "Writing HTML file (.html) with report contents - rmarkdown")
      css_fname <- NULL
      if (tier_model == "cpsr") {
        markdown_input <- system.file("templates",
                                      "cpsr_rmarkdown_report.Rmd",
                                      package = "pcgrr")
        css_fname <- system.file("templates", "cpsr.css",
                                 package = "pcgrr")
      }
      rmarkdown::render(
        markdown_input,
        output_format =
          rmarkdown::html_document(theme = report_theme,
                                   toc = T, toc_depth = 3,
                                   toc_float = T,
                                   number_sections = F,
                                   css = css_fname,
                                   includes =
                                     rmarkdown::includes(
                                       after_body = disclaimer)),
        output_file = fnames[["html"]],
        output_dir = project_directory,
        clean = T,
        intermediates_dir = project_directory, quiet = T)
    }
  }
  if (output_format == "json") {
    if (!is.null(report[["cna_plot"]][["png"]])) {
      report[["cna_plot"]][["png"]] <- NULL
    }
    if (!is.null(report[["tmb"]][["tcga_tmb"]])) {
      report[["tmb"]][["tcga_tmb"]] <- NULL
    }
    rlogging::message("------")
    rlogging::message("Writing JSON file (.json) with report contents")
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
      rlogging::message("------")
      if (tier_model == "pcgr_acmg") {
        rlogging::message(
          paste0("Writing SNV/InDel tab-separated output file with ",
                 "PCGR annotations - ('",
                 output_format_slim, "')"))
      }
      if (tier_model == "cpsr") {
        rlogging::message(
          paste0("Writing SNV/InDel tab-separated output file ",
                 "with CPSR annotations - ('",
                 output_format_slim, "')"))
      }
      write.table(
        report[["content"]][["snv_indel"]][["variant_set"]][[output_format_slim]],
        file = fnames[[output_format]], sep = "\t", col.names = T,
        row.names = F, quote = F)
    }
  }

  if (output_format == "msigs_tsv") {
    if (
      NROW(report[["content"]][["m_signature_mp"]][["result"]][["tsv"]]) > 0) {
      rlogging::message("------")
      rlogging::message(
        "Writing tab-separated output file with details ",
        "of contributing mutational signatures - ('tsv')")
      write.table(report[["content"]][["m_signature_mp"]][["result"]][["tsv"]],
                  file = fnames[[output_format]], sep = "\t", col.names = T,
                  row.names = F, quote = F)
    }
  }

  if (output_format == "cna_tsv") {
    if (NROW(report[["content"]][["cna"]][["variant_set"]][["tsv"]]) > 0) {
      rlogging::message("------")
      rlogging::message(
        "Writing CNA tab-separated output file with PCGR annotations (.tsv.gz)")
      write.table(report[["content"]][["cna"]][["variant_set"]][["tsv"]],
                  file = fnames[["cna_tsv"]], sep = "\t", col.names = T,
                  row.names = F, quote = F)
      gzip_command <- paste0("gzip -f ", fnames[["cna_tsv"]])
      system(gzip_command, intern = F)
    }
  }

}

