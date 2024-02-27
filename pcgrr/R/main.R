#' Function that generates all contents of the cancer
#' genome report (PCGR)
#'
#' @param yaml_fname Name of PCGR settings file produced by
#' pre-reporting Python workflow (yaml)
#'
#' @export

generate_pcgr_report <-
  function(yaml_fname = NULL) {

    invisible(assertthat::assert_that(
      !is.null(yaml_fname),
      msg = "Object 'yaml_fname' cannot be NULL"
    ))
    pcgrr::check_file_exists(yaml_fname)

    ## Initialize report list object -
    ## 1) settings - run configurations
    ## 2) content (for display in HTML)
    ## 3) reference data
    rep <- pcgrr::init_report(
      yaml_fname = yaml_fname,
      report_mode = "PCGR")

    rep[['content']][['assay_properties']] <-
      rep[['settings']]$conf$assay_properties

    rep[['content']][['sample_properties']] <-
      rep[['settings']]$conf$sample_properties

    settings <- rep$settings
    ref_data <- rep$ref_data

    callset_snv <-
      pcgrr::load_somatic_snv_indel(
        fname = settings$molecular_data$fname_mut_tsv,
        ref_data = ref_data,
        settings = settings
      )

    callset_cna <- NULL
    if (settings$molecular_data$fname_cna_tsv != "None") {
      callset_cna <-
        pcgrr::load_somatic_cna(
          fname = settings$molecular_data$fname_cna_tsv,
          ref_data = ref_data,
          settings = settings
        )

      rep[['content']][['cna']][['callset']] <-
        callset_cna
      rep[['content']][['cna']][['vstats']] <-
        pcgrr::variant_stats_report(
          callset = callset_cna,
          vartype = "cna",
          name = "vstats")[['vstats']]
      rep[['content']][['cna']][['cnaqc']] <-
        pcgrr::make_cnaqc_object(
          callset_cna = callset_cna,
          callset_snv = callset_snv,
          settings = settings
        )
    }

    conf_somatic_snv <-
      settings$conf$somatic_snv
    conf_somatic_cna <-
      settings$conf$somatic_cna
    conf_other <-
      settings$conf$other
    assay_properties <-
      settings$conf$assay_properties
    sample_properties <-
      settings$conf$sample_properties

    ## Retrieve relevant clinical trials for the tumor type in question

    if (as.logical(settings$conf$clinicaltrials$run) == T) {
      # pcg_report_trials <-
      #   pcgrr::generate_report_data_trials(
      #     ref_data = ref_data,
      #     settings = settings)
      # ## Update genome report with trial data
      # pcg_report <-
      #   pcgrr::update_report(pcg_report, pcg_report_trials,
      #                        a_elem = "clinicaltrials")
    }

    if (NROW(callset_snv$variant) > 0) {

      rep[['content']][['snv_indel']][['callset']] <-
        callset_snv
      rep[['content']][['snv_indel']][['vstats']] <-
        pcgrr::variant_stats_report(
          callset = callset_snv,
          vartype = "snv_indel",
          name = "vstats")[['vstats']]

      ## Estimate contribution of mutational signatures
      if (conf_somatic_snv[["mutational_signatures"]][["run"]] == T) {

        ## Generate report data for mutational signatures assessment
        pcg_report_signatures <-
          pcgrr::generate_report_data_signatures(
            variant_set = callset_snv$variant,
            settings = settings,
            ref_data = ref_data)

        ## Update genome report with signature data
        rep <- pcgrr::update_report(
          rep,
          pcg_report_signatures,
          a_elem = "mutational_signatures")


        ## Generate report data for rainfall plot
        pcg_report_rainfall <-
          pcgrr::generate_report_data_rainfall(
            variant_set = callset_snv$variant,
            build = settings$genome_assembly)

        ## Update genome report with rainfall data
        rep <-
          pcgrr::update_report(
            rep,
            pcg_report_rainfall,
            a_elem = "rainfall")

        ## Generate report data for kataegis events
        ## (for WES/WGS runs)
        if (stringr::str_detect(
          assay_properties[["type"]],
          "WGS|WES")) {
          pcg_report_kataegis <-
            pcgrr::generate_report_data_kataegis(
              variant_set = callset_snv$variant,
              sample_name = settings$sample_id,
              build = settings$genome_assembly)

          ## Update genome report with kataegis data
          rep <- pcgrr::update_report(
            rep,
            pcg_report_kataegis,
            a_elem = "kataegis")
        }
      }

      ## If assay is Tumor-Control and WES/WGS -
      ## perform MSI prediction
      if (as.logical(conf_somatic_snv[['msi']][['run']]) == T &
          stringr::str_detect(
            assay_properties[["type"]], "WGS|WES") &
          as.logical(
            assay_properties[["vcf_tumor_only"]]) == FALSE) {

        ## Generate report data for MSI classification
        pcg_report_msi <-
          pcgrr::generate_report_data_msi(
            variant_set = callset_snv$variant,
            ref_data = ref_data,
            settings = settings)

        ## Update genome report with MSI classiciation
        rep <-
          pcgrr::update_report(
            rep,
            pcg_report_msi,
            a_elem = "msi")
      }

      ## Generate report contents for analysis of
      ## mutational burden (TMB)
      if (as.logical(conf_somatic_snv[['tmb']][['run']]) == T) {

        pcg_report_tmb <-
          pcgrr::generate_report_data_tmb(
            settings = settings)

        ## Update genome report with TMB data
        rep[["content"]][["tmb"]][["eval"]] <-
          pcg_report_tmb[["eval"]]
        rep[["content"]][["tmb"]][["sample_estimate"]] <-
          pcg_report_tmb[["sample_estimate"]]

      }
    }else{
      rep[["content"]][["snv_indel"]][["zero"]] <- TRUE
      rep[["metadata"]][["config"]][["other"]][["list_noncoding"]] <- FALSE
    }

    # if (!is.null(cpsr_report_fname)) {
    #   pcg_report[["content"]][["cpsr"]][['eval']] <- TRUE
    #
    #   pcg_report[['content']][['cpsr']][['report']] <-
    #     jsonlite::fromJSON(
    #       gzfile(cpsr_report_fname)
    #     )
    #
    #   ## append report elements in pcg_report[['content']][['cpsr]][['cpsr_json']]
    # }


    #}

    #pcg_report_value_box <- pcgrr::generate_report_data_value_box(
    #  pcg_report, pcgr_data, sample_name, config)
    #pcg_report <- pcgrr::update_report(
    #  pcg_report, pcg_report_value_box,
    #  a_elem = "value_box")

    return(rep)
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
  pcgrr::log4r_info(paste0(
    "Generating tiered set of result variants for output",
    " in tab-separated values (TSV) file"))
  tsv_variants <- NULL
  for (tier in c("tier1", "tier2", "tier3", "tier4", "tier5")) {
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

      tierset <- tierset |>
        dplyr::select(dplyr::any_of(tsv_columns)) |>
        dplyr::distinct()

      tsv_variants <- dplyr::bind_rows(tsv_variants, tierset)
    }
  }
  tsv_variants$OFFICIAL_GENENAME <-
    unlist(lapply(stringr::str_match_all(tsv_variants$OFFICIAL_GENENAME, ">.+<"),
                  paste, collapse = ","))
  tsv_variants$OFFICIAL_GENENAME <-
    stringr::str_replace_all(tsv_variants$OFFICIAL_GENENAME, ">|<", "")
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
  tsv_variants <- tsv_variants |> dplyr::distinct()

  return(tsv_variants)
}


#' Function that writes contents of PCGR object to various output formats
#' (Rmarkdown/flexdashboard HTML reports, JSON, tab-separated etc)
#'
#' @param report List object with all report data (PCGR/CPSR), settings etc.
#' @param tier_model type of tier model
#' @param output_format contents/file format of output
#' (html/json/tsv/cna_tsv etc)
#' @param flexdb logical indicating if HTML output should be dashboard

#' @export
write_report_output <- function(report,
                                tier_model = "pcgr_acmg",
                                output_format = "html",
                                flexdb = FALSE) {

  settings <- report[['settings']]
  project_directory <- settings[['output_dir']]
  sample_name <- settings[['sample_id']]
  genome_assembly <- settings[['genome_assembly']]

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
  # fnames[["cna_tsv"]] <-
  #   file.path(project_directory,
  #             paste0(sample_fname_pattern,
  #                    ".cna_segments.tsv"))
  # fnames[["json"]] <-
  #   file.path(project_directory,
  #             paste0(sample_fname_pattern, ".json"))
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
  cpsr_tmpl <- system.file("templates", package = "cpsr")
  disclaimer <- file.path(cpsr_tmpl, "disclaimer_predisposition.md")
  markdown_input <- file.path(cpsr_tmpl, "cpsr_rmarkdown_report.Rmd")
  css_fname <- file.path(cpsr_tmpl, "cpsr.css")
  report_theme <-
    settings[["conf"]][["visual_reporting"]][["visual_theme"]]

  ## Somatic/tumor report settings
  if (tier_model == "pcgr_acmg") {
    pcgrr_tmpl <- system.file("templates", package = "pcgrr")

    disclaimer <- file.path(pcgrr_tmpl, "disclaimer.md")
    assay_props <-
      settings[["conf"]][["assay_properties"]]
    sequencing_assay <-
      assay_props[["type"]]

    ## Flexdashboard layout
    sequencing_design <- "Tumor-Control"
    markdown_input <- file.path(pcgrr_tmpl, "pcgr_flexdb_report.Rmd")
    css_fname <- file.path(pcgrr_tmpl, "pcgr_flexdb_tumor_control.css")

    ## Rmarkdown layout
    if (flexdb == FALSE) {
      markdown_input <- file.path(pcgrr_tmpl, "pcgr_rmarkdown_report.Rmd")
      css_fname <- file.path(pcgrr_tmpl, "pcgr_rmarkdown_tumor_control.css")
    }

    ## Tumor-only settings (CSS)
    if (assay_props[["vcf_tumor_only"]] == T) {
      sequencing_design <- "Tumor-Only"
      css_fname <- file.path(pcgrr_tmpl, "pcgr_flexdb_tumor_only.css")

      if (flexdb == FALSE) {
        css_fname <- file.path(pcgrr_tmpl, "pcgr_rmarkdown_tumor_only.css")
      }
    }
  }

  if (output_format == "html") {

    if (flexdb == T & tier_model == "pcgr_acmg") {
      pcgrr::log4r_info("------")
      pcgrr::log4r_info(
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
            favicon = system.file(
              "templates","favicon-16x16.png",
              package = "pcgrr"),
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
      if (tier_model == 'cpsr') {
        toc_float <-
          list(collapsed = FALSE,
               smooth_scroll = FALSE,
               print = TRUE)
      }

      ## If nonfloating TOC is chosen (PCGR/CPSR), set toc_float to FALSE
      nonfloating_toc <-
        as.logical(settings[["conf"]][["visual_reporting"]][["nonfloating_toc"]])
      if (nonfloating_toc == T) {
        toc_float <- F
      }

      disclaimer <- system.file(
        "templates",
        "disclaimer.md",
        package = "pcgrr")

      header <- system.file(
         "templates",
         "_header.html",
         package = "pcgrr")
      if (tier_model == "cpsr") {
        header <- system.file(
          "templates",
          "_header.html",
          package = "cpsr")
      }

      pcgrr::log4r_info("------")
      pcgrr::log4r_info(paste0(
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
                in_header = header,
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
    pcgrr::log4r_info("------")
    pcgrr::log4r_info("Writing JSON file (.json) with key report contents")

    report_strip <- report

    if (tier_model != "cpsr") {
      if (!is.null(report_strip$content$rainfall)) {
        report_strip$content$rainfall <- NULL
      }
      if (!is.null(report_strip$content$tmb)) {
        report_strip$content$tmb$tcga_tmb <- NULL
      }
      if (!is.null(report_strip$content$clinicaltrials)) {
        report_strip$content$clinicaltrials <- NULL
      }
      if (!is.null(report_strip$content$msi)) {
        if (!is.null(report_strip$content$msi$prediction)) {
          report_strip$content$msi$prediction$tcga_dataset <- NULL
        }
      }

      if (!is.null(report_strip$content$snv_indel$disp)) {
        report_strip$content$snv_indel$disp <- NULL
      }

      if (!is.null(report_strip$content$snv_indel$variant_set)) {
        if (!is.null(report_strip$content$snv_indel$variant_set$maf)) {
          report_strip$content$snv_indel$variant_set$maf <- NULL
        }
      }

      key_tsv_cols <- c("GENOMIC_CHANGE",
                        "VARIANT_CLASS",
                        "SYMBOL",
                        "ENTREZGENE",
                        "ENSEMBL_TRANSCRIPT_ID",
                        "TUMOR_SUPPRESSOR",
                        "ONCOGENE",
                        "CONSEQUENCE",
                        "PROTEIN_CHANGE",
                        "PROTEIN_DOMAIN",
                        "CODING_STATUS",
                        "EXONIC_STATUS",
                        "HGVSp",
                        "MUTATION_HOTSPOT",
                        "DBSNPRSID",
                        "COSMIC_MUTATION_ID",
                        "CALL_CONFIDENCE",
                        "DP_TUMOR",
                        "AF_TUMOR",
                        "DP_CONTROL",
                        "AF_CONTROL",
                        "TIER")

      if (!is.null(report_strip$content$snv_indel$variant_set)) {

        for(o in c('tsv')) {

          if (!is.null(report_strip$content$snv_indel$variant_set[[o]])) {

            if (nrow(report_strip$content$snv_indel$variant_set[[o]]) == 0) {
              next
            }
            assertable::assert_colnames(
              report_strip$content$snv_indel$variant_set[[o]],
              colnames = key_tsv_cols,
              only_colnames = F,
              quiet = T
            )

            report_strip$content$snv_indel$variant_set[[o]] <-
              dplyr::select(
                report_strip$content$snv_indel$variant_set[[o]],
                dplyr::any_of(key_tsv_cols)
              )

          }
        }
      }

    } ## if tier_model != "cpsr"


    size <- format(utils::object.size(report_strip), units = "auto")
    #hsize <- R.utils::hsize.object_size(size)
    pcgrr::log4r_info(paste0("Size of PCGR report object for JSON output: ", size))


    ## NOTE: set max size of report object to 750 Mb - have not figured out
    ## what the exact size should be for jsonlite::toJSON to succeed/fail
    if (utils::object.size(report_strip) < 750000000) {

      pcgr_json <- jsonlite::toJSON(
        report_strip, pretty = T, na = "string",
        null = "null", force = T)
      write(pcgr_json, fnames[["json"]])
      gzip_command <- paste0("gzip -f ", fnames[["json"]])
      system(gzip_command, intern = F)
    }else{
      pcgrr::log4r_info("JSON output not possible - report contents too large (> 750Mb)")

    }
  }

  if (output_format == "snv_tsv" | output_format == "snv_tsv_unfiltered") {
    output_format_slim <- stringr::str_replace(output_format, "snv_", "")
    if (NROW(
      report[["content"]][["snv_indel"]][["variant_set"]][[output_format_slim]]) > 0) {
      pcgrr::log4r_info("------")
      if (tier_model == "pcgr_acmg") {
        pcgrr::log4r_info(
          paste0(
            "Writing SNV/InDel tab-separated output file with ",
            "PCGR annotations - ('",
            output_format_slim, "')"))
      }
      if (tier_model == "cpsr") {
        pcgrr::log4r_info(
          paste0(
            "Writing SNV/InDel tab-separated output file ",
            "with CPSR annotations - ('",
            output_format_slim, "')"))
      }
      utils::write.table(
        report[["content"]][["snv_indel"]][["variant_set"]][[output_format_slim]],
        file = fnames[[output_format]], sep = "\t", col.names = T,
        row.names = F, quote = F)

      # if (tier_model == "pcgr_acmg") {
      #   pcgrr::log4r_info(
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
      NROW(
        report[["content"]][["mutational_signatures"]][["result"]][["tsv"]]) > 0) {
      pcgrr::log4r_info("------")
      pcgrr::log4r_info(paste0(
        "Writing tab-separated output file with details ",
        "of contributing mutational signatures - ('tsv')"))
      utils::write.table(
        report[["content"]][["mutational_signatures"]][["result"]][["tsv"]],
        file = fnames[[output_format]], sep = "\t", col.names = T,
        row.names = F, quote = F)
    }
  }

  if (output_format == "cna_tsv") {
    if (NROW(report[["content"]][["cna"]][["variant_set"]][["tsv"]]) > 0) {
      pcgrr::log4r_info("------")
      pcgrr::log4r_info(
        "Writing CNA tab-separated output file with PCGR annotations (.tsv.gz)")
      utils::write.table(report[["content"]][["cna"]][["variant_set"]][["tsv"]],
                  file = fnames[["cna_tsv"]], sep = "\t", col.names = T,
                  row.names = F, quote = F)
      gzip_command <- paste0("gzip -f ", fnames[["cna_tsv"]])
      system(gzip_command, intern = F)
    }
  }

}

