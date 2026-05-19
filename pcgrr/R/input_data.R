#' Function that reads and validates fully annotated CNA data (segments and genes)
#' from PCGR pre-reporting pipeline
#'
#' @param fname_cna_segment Path to file with pre-processed CNA segments
#' @param fname_cna_gene Path to file with pre-processed CNA gene-level data
#' @param ref_data PCGR reference data object
#' @param settings PCGR run/configuration settings
#'
#' @export
load_somatic_cna <- function(
    fname_cna_segment = NULL,
    fname_cna_gene = NULL,
    ref_data = NULL,
    settings = NULL) {

  log4r_info("------")
  log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - somatic copy number aberrations"))

  hgname <- "hg38"
  if (settings$genome_assembly == "grch37") {
    hgname <- "hg19"
  }

  tumor_site <-
    settings[['conf']][['sample_properties']][['site']]

  ## read segments
  check_file_exists(fname_cna_segment)
  check_file_exists(fname_cna_gene)
  segments_raw <- suppressWarnings(
    as.data.frame(
      readr::read_tsv(
        file = fname_cna_segment,
        na = c(".","NA"),
        show_col_types = F,
        progress = F
      )
    )
  )

  compulsary_cols <-
    names(data_coltype_defs$cna_somatic_segment_raw$cols)

  raw_col_check <-
    rlang::has_name(segments_raw, compulsary_cols)
  if (FALSE %in% raw_col_check) {
    missing_cols <-
      compulsary_cols[!raw_col_check]
    log4r_fatal(
      paste0("Missing required columns in input file ",
             fname_cna_segment, " - ",
             paste(missing_cols, collapse=", ")))
  }

  segments <- suppressWarnings(
    as.data.frame(
      readr::read_tsv(
        file = fname_cna_segment,
        col_types =
          data_coltype_defs$cna_somatic_segment_raw$cols,
        na = c(".","NA"),
        progress = F
      )
    )) |>
    tidyr::separate(
      col = "SEGMENT_NAME",
      into = c("SEGMENT_ID",
               "CN_MAJOR",
               "CN_MINOR",
               "ARM",
               "CYTOBAND",
               "EVENT_TYPE",
               "VARIANT_CLASS",
               "LOH"),
      sep = "\\|",
      remove = T
    ) |>
    dplyr::mutate(
      CN_MAJOR = as.integer(.data$CN_MAJOR),
      CN_MINOR = as.integer(.data$CN_MINOR),
      CN_TOTAL = as.integer(.data$CN_MAJOR + .data$CN_MINOR),
      LOH = dplyr::if_else(
        nchar(.data$LOH) == 0,
        NA_character_,
        as.character(.data$LOH)
      )
    )


  callset_cna <- load_dna_variants(
    fname = fname_cna_gene,
    cols = data_coltype_defs$cna_somatic_gene_raw,
    ref_data = ref_data,
    settings = settings,
    vartype = 'cna',
    primary_site =
      tumor_site,
    retained_info_tags = "None",
    variant_origin = "Somatic")

  callset_cna[['segment']] <- segments
  if (NROW(callset_cna$variant) > 0) {
    callset_cna[['variant']] <- callset_cna[['variant']] |>
      dplyr::mutate(CN_TOTAL =
                      as.integer(.data$CN_MAJOR + .data$CN_MINOR)) |>
      dplyr::mutate(SAMPLE_ALTERATION = paste(
        .data$SYMBOL, .data$VARIANT_CLASS)) |>
      dplyr::mutate(CYTOBAND = paste0(
        "chr", .data$CHROM,":" ,.data$CYTOBAND
      )) |>
      append_cancer_association_ranks(
        ref_data = ref_data,
        primary_site = tumor_site) |>
      append_styled_cna_vclass(
        colname = "VARIANT_CLASS_STYLED") |>
      append_targeted_drug_annotations(
        ref_data = ref_data,
        primary_site = tumor_site) |>
      dplyr::arrange(
        .data$ACTIONABILITY_TIER,
        dplyr::desc(.data$TISSUE_ASSOC_RANK),
        dplyr::desc(.data$GLOBAL_ASSOC_RANK)) |>
      order_variants(pos_var = 'SEGMENT_START') |>
      exclude_non_chrom_variants()

    log4r_info(
      "Generating data frame with hyperlinked variant/gene annotations")

    callset_cna[['variant_display']] <- callset_cna[['variant']] |>
      dplyr::mutate(
        SEGMENT = glue::glue(
          "<a href='http://genome.ucsc.edu/cgi-bin/hgTracks?db={hgname}&position=",
          "chr{CHROM}:{SEGMENT_START}-{SEGMENT_END}' target='_blank'>",
          "chr{CHROM}:{SEGMENT_START}-{SEGMENT_END}</a>"
        )
      ) |>
      append_cancer_gene_evidence(
        ref_data = ref_data
      ) |>
      tidyr::separate_rows(
        "TRANSCRIPT_OVERLAP",
        sep=", "
      ) |>
      tidyr::separate(
        "TRANSCRIPT_OVERLAP",
        c("ENSEMBL_TRANSCRIPT_ID",
          "REFSEQ_TRANSCRIPT_ID",
          "TRANSCRIPT_START",
          "TRANSCRIPT_END",
          "TRANSCRIPT_OVERLAP_PERCENT"),
        sep = "\\|", remove = T
        ) |>
      append_annotation_links(
        vartype = "cna"
      ) |>
      tidyr::unite(
        "TRANSCRIPT_OVERLAP",
        c("ENSEMBL_TRANSCRIPT_ID",
          "REFSEQ_TRANSCRIPT_ID",
          "TRANSCRIPT_OVERLAP_PERCENT"),
        sep="|", remove = T
      ) |>
      order_variants(pos_var = 'SEGMENT_START') |>
      dplyr::select(
        -dplyr::ends_with(c("_RAW","_END","_START"))
      ) |>
      dplyr::group_by(
        dplyr::across(-c("TRANSCRIPT_OVERLAP"))) |>
      dplyr::summarise(
        TRANSCRIPT_OVERLAP = paste(
          .data$TRANSCRIPT_OVERLAP, collapse=", "),
        .groups = "drop"
      ) |>
      dplyr::distinct() |>
      dplyr::rename(
        TARGETED_INHIBITORS = "TARGETED_INHIBITORS2",
        TARGETED_INHIBITORS_ALL = "TARGETED_INHIBITORS_ALL2"
      ) |>
      dplyr::arrange(
        .data$ACTIONABILITY_TIER,
        dplyr::desc(.data$GLOBAL_ASSOC_RANK),
        dplyr::desc(.data$TISSUE_ASSOC_RANK))
  }

  return(callset_cna)

}

#' Function that reads and validates an annotated somatic SNV/InDel
#' file from PCGR pre-reporting pipeline
#'
#' @param fname Path to file with pre-processed somatic SNV/InDel variants
#' @param ref_data PCGR reference data object
#' @param settings PCGR run/configuration settings
#' @param simulate_vaf_dp Internal/test use only. If TRUE and VAF_TUMOR is
#'   entirely missing, replace it with random values drawn from
#'   Uniform(0.01, 0.99). Never set this in production runs.
#'
#' @export
load_somatic_snv_indel <- function(
    fname = NA,
    ref_data = NULL,
    settings = NULL,
    simulate_vaf_dp = TRUE) {

  log4r_info("------")
  log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - somatic SNV/InDels"))

  tumor_site <-
    settings[['conf']][['sample_properties']][['site']]

  callset <- load_dna_variants(
    fname = fname,
    cols = data_coltype_defs$snv_indel_somatic_raw,
    ref_data = ref_data,
    settings = settings,
    vartype = 'snv_indel',
    primary_site = tumor_site,
    retained_info_tags =
      settings[['conf']][['other']]$retained_vcf_info_tags,
    variant_origin = "Somatic")

  #callset[['variant_unfiltered']] <- data.frame()

  ## Internal/test only: simulate VAF_TUMOR when it is entirely absent
  if (isTRUE(simulate_vaf_dp) &&
      "VAF_TUMOR" %in% colnames(callset[['variant']]) &&
      all(is.na(callset[['variant']]$VAF_TUMOR)) &&
      NROW(callset[['variant']]) > 0) {
    set.seed(42L)
    callset[['variant']]$VAF_TUMOR <-
      round(stats::runif(NROW(callset[['variant']]), min = 0.01, max = 0.99), 3)
    log4r_info(
      "VAF_TUMOR simulated with uniform (0.01, 0.99) values [TEST MODE ONLY]")
  }

  ## simulate also DP_tumor when it is entirely absent,
  ## random number from 15 to 150
  if (isTRUE(simulate_vaf_dp) &&
      "DP_TUMOR" %in% colnames(callset[['variant']]) &&
      all(is.na(callset[['variant']]$DP_TUMOR)) &&
      NROW(callset[['variant']]) > 0) {
    set.seed(42L)
    callset[['variant']]$DP_TUMOR <-
      sample(15:150, NROW(callset[['variant']]), replace = TRUE)
    log4r_info(
      "DP_TUMOR simulated with random integer values from 15 to 150 [TEST MODE ONLY]")
  }

  callset[['variant']] <- callset[['variant']] |>
    append_protein_domains(
      ref_data = ref_data) |>
    append_cancer_association_ranks(
      ref_data = ref_data,
      primary_site = tumor_site) |>
    append_targeted_drug_annotations(
      ref_data = ref_data,
      primary_site = tumor_site) |>
    append_alteration_name() |>
    order_variants(pos_var = 'POS') |>
    exclude_non_chrom_variants() |>
    filter_read_support(config = settings$conf)


  ## Tumor-only input
  if (as.logical(settings$conf$assay_properties$vcf_tumor_only) == TRUE) {
    if (NROW(callset[['variant']]) > 0) {

      callset[['variant_unfiltered']] <-
        callset[['variant']] |>
        ## assign evidence tags for germline/somatic state of variants,
        ## partially based on user-defined options
        ## (population allele frequency thresholds)
          assign_somatic_germline_evidence(
            settings = settings) |>

        ## assign somatic variant classification/status based on accumulation
        ## of evidence tags and user-defined options
          assign_somatic_classification(
            settings = settings)

      ## Assign calls to filtered callset (SOMATIC_CLASSIFICATION = SOMATIC)
      if ("SOMATIC_CLASSIFICATION" %in%
         colnames(callset[['variant_unfiltered']])) {
        callset[['variant']] <-
          callset[['variant_unfiltered']] |>
          dplyr::filter(
            .data$SOMATIC_CLASSIFICATION == "SOMATIC")

        ## Issue warning if clinically actionable variants are filtered
        ## with current filtering settings
        actionable_filtered <-
          callset[['variant_unfiltered']] |>
          dplyr::filter(
            !is.na(.data$ACTIONABILITY_TIER) &
              .data$ACTIONABILITY_TIER <= 2 &
              .data$SOMATIC_CLASSIFICATION != "SOMATIC")

        if (NROW(actionable_filtered) > 0) {
          log4r_warn(
            paste0(
              "A total of n = ", NROW(actionable_filtered),
              " clinically actionable ",
              "variants were filtered as likely germline events"))
        }

        n_excluded_calls_germline <- as.character(
          formattable::comma(
            (NROW(callset$variant_unfiltered) - NROW(callset$variant)),
            digits = 0))

        log4r_info(paste0(
          "Tumor-only variant filtering based on multiple criteria - ",
          get_tumor_only_filtering_criteria(
            conf = settings$conf))
        )

        log4r_info(
          paste0(
            "Excluded n = ",
            n_excluded_calls_germline,
            " putative germline variants after applying the criteria above"))

        if (NROW(callset$variant) == 0) {
          log4r_warn(
            "NO (n = 0) somatic variants remain after filtering of putative germline events")
        }

        ## Option to exclude non-exonic variants
        if (as.logical(
          settings$conf$somatic_snv$tumor_only[["exclude_nonexonic"]]) == TRUE &
          NROW(callset[['variant']]) > 0 &
          "EXONIC_STATUS" %in% colnames(callset[['variant']])) {

          n_exonic_nonexonic <- NROW(callset[['variant']])
          log4r_info(paste0(
            "Tumor-only variant filtering based on exonic status only"
          ))

          callset[['variant']] <- callset[['variant']] |>
            dplyr::filter(
              .data$EXONIC_STATUS == "exonic")

          n_excluded_calls_nonexonic <- as.character(
            formattable::comma(
              (n_exonic_nonexonic - NROW(callset$variant)),
              digits = 0))

          log4r_info(
            paste0(
              "Excluded n = ",
              n_excluded_calls_nonexonic,
              " variants aftering filtering of non-exonic variants"))

          if (NROW(callset$variant) == 0) {
            log4r_warn(
              "NO (n = 0) somatic variants remain after filtering of non-exonic variants")
          }
        }

        ## filter also MAF file if provided
        filter_maf_file(
          callset = callset,
          settings = settings)

      }else{
        log4r_fatal(
          "Variant data.frame is lacking a 'SOMATIC_CLASSIFICATION' column")
      }
    }
  }

  if (NROW(callset[['variant']]) > 0) {
    callset[['variant']] <- callset[['variant']] |>
      dplyr::arrange(
        .data$ACTIONABILITY_TIER,
        dplyr::desc(.data$ONCOGENICITY_SCORE),
        dplyr::desc(.data$GLOBAL_ASSOC_RANK),
        dplyr::desc(.data$TISSUE_ASSOC_RANK))


    ## Make data frame with columns for display
    ## in HTML output

    log4r_info(
      "Generating data frame with hyperlinked variant/gene annotations")


    callset[['variant_display']] <- callset[['variant']] |>
      append_cancer_gene_evidence(
        ref_data = ref_data) |>
      append_oncogenicity_docs(
        ref_data = ref_data) |>
      append_dbmts_var_link() |>
      append_tcga_var_link() |>
      append_annotation_links() |>
      append_dbnsfp_var_link() |>
      append_tfbs_annotation() |>
      dplyr::select(
        -dplyr::contains("_RAW")
      ) |>

      dplyr::rename(
        TARGETED_INHIBITORS = "TARGETED_INHIBITORS2",
        TARGETED_INHIBITORS_ALL = "TARGETED_INHIBITORS_ALL2"
      ) |>
      dplyr::mutate(
        CONSEQUENCE = stringr::str_replace_all(
          .data$CONSEQUENCE,"&",", ")) |>
      dplyr::mutate(
        MUTATION_HOTSPOT_CANCERTYPE = stringr::str_replace_all(
          .data$MUTATION_HOTSPOT_CANCERTYPE,",",", ")) |>
      order_variants(pos_var = 'POS')

  }

  return(callset)


}

#' Function that reads CPSR-classified variants from a TSV file
#'
#' @param fname_cpsr_tsv Path to raw input file with CPSR-classified SNVs/InDels
#' @param fname_cpsr_yaml Path to YAML configuration file for CPSR analysis
#' @param cols column type definitions of raw input file
#' @param ignore_vus logical indicating if VUS should be ignored in report
#' @param ref_data PCGR reference data object
#' @param settings PCGR run/configuration settings
#'
#' @export
#'
load_cpsr_classified_variants <- function(
    fname_cpsr_tsv = NA,
    fname_cpsr_yaml = NA,
    cols = NULL,
    ignore_vus = FALSE,
    ref_data = NULL,
    settings = NULL) {

  log4r_info("------")
  log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - germline SNV/InDels (CPSR-classified)"))

  check_file_exists(fname_cpsr_tsv)
  check_file_exists(fname_cpsr_yaml)

  if (!file.exists(fname_cpsr_yaml)) {
    log4r_fatal(
      paste0("YAML file '", fname_cpsr_yaml, "' does not exist - exiting"))
  }
  cpsr_yaml <- yaml::read_yaml(fname_cpsr_yaml)
  if ("conf" %in% names(cpsr_yaml) == FALSE) {
    log4r_fatal(
      paste0(
        "YAML file '", fname_cpsr_yaml,
        "' does not contain a 'conf' section - exiting"))
  }
  if ("sample_id" %in% names(cpsr_yaml) == FALSE) {
    log4r_fatal(
      paste0(
        "YAML file '", fname_cpsr_yaml,
        "' does not contain a 'sample_id' variable - exiting"))
  }
  if ("gene_panel" %in% names(cpsr_yaml$conf) == FALSE) {
    log4r_fatal(
      paste0(
        "YAML file '", fname_cpsr_yaml,
        "' does not contain a 'conf->gene_panel' section - exiting"))
  }
  panel_info <- list()
  panel_info[['description']] <-
    cpsr_yaml$conf$gene_panel[['description']]
  panel_info[['description_trait']] <-
    cpsr_yaml$conf$gene_panel[['description_trait']]
  panel_info[['url']] <- cpsr_yaml$conf$gene_panel[['url']]
  panel_info[['panel_id']] <- cpsr_yaml$conf$gene_panel[['panel_id']]
  sample_id <- cpsr_yaml$sample_id

  callset <- load_dna_variants(
    fname = fname_cpsr_tsv,
    cols = cols,
    ref_data = ref_data,
    settings = settings,
    vartype = 'snv_indel',
    variant_origin = 'Germline')

  callset[['variant_display']] <- callset[['variant']] |>
    append_cancer_gene_evidence(
      ref_data = ref_data) |>
    dplyr::mutate(
      CLINVAR_TRAITS_ALL = paste(
        stringr::str_to_title(.data$CLINVAR_VARIANT_ORIGIN),
        .data$CLINVAR_PHENOTYPE,
        sep = " - ")) |>
    append_annotation_links() |>
    dplyr::select(
      -dplyr::contains("_RAW")
    ) |>
    dplyr::mutate(
      CONSEQUENCE = stringr::str_replace_all(
        .data$CONSEQUENCE,"&",", ")) |>
    dplyr::rename(
      #SOURCE = .data$CPSR_CLASSIFICATION_SOURCE,
      CLINICAL_SIGNIFICANCE = .data$CLASSIFICATION
    ) |>
    dplyr::select(
      -dplyr::any_of(
        c("CLINVAR_TRAITS_ALL",
         "CLINVAR_VARIANT_ORIGIN",
         "CLINVAR_PHENOTYPE",
         "PFAM_DOMAIN",
         "CANCERGENE_EVIDENCE",
         "PFAM_DOMAIN_NAME",
         "PROTEIN_CHANGE",
         "CLINVAR_MSID",
         "ENTREZGENE"))) |>
    dplyr::select(
      dplyr::any_of(
        c("SYMBOL","ALTERATION",
          "GENOTYPE","CONSEQUENCE",
        "CLINICAL_SIGNIFICANCE",
        "ASSERTION_AUTHORITY",
        "PROTEIN_DOMAIN",
        "HGVSc", "HGVSc_RefSeq",
        "HGVSp", "CDS_CHANGE",
        "EFFECT_PREDICTIONS","SPLICE_EFFECT",
        "CODING_STATUS",
        "LOSS_OF_FUNCTION", 
        "DP_CONTROL",
        "VARIANT_CLASS","GENENAME",
        "ENSEMBL_GENE_ID",
        "ENSEMBL_TRANSCRIPT_ID",
        "REFSEQ_TRANSCRIPT_ID",
        "DBSNP_RSID",
        "CLINVAR",
        "CLINVAR_CLASSIFICATION",
        "CLINVAR_CONFLICTED",
        "CLINVAR_GOLD_STARS",
        "ASSERTION_RATIONALE",
        "CPSR_PATHOGENICITY_SCORE",
        "CPSR_CLASSIFICATION",
        "ACMG_CODE")),
      dplyr::everything()
    )

  if (NROW(callset[['variant_display']]) > 0) {
    callset[['variant_display']] <- callset[['variant_display']] |>
      dplyr::filter((.data$CLINICAL_SIGNIFICANCE == "Pathogenic" |
                      .data$CLINICAL_SIGNIFICANCE == "Likely Pathogenic" |
                      .data$CLINICAL_SIGNIFICANCE == "VUS") &
                      .data$CODING_STATUS == "coding")|>
      dplyr::distinct()
    if (ignore_vus == TRUE) {
      callset[['variant_display']] <- callset[['variant_display']] |>
        dplyr::filter(.data$CLINICAL_SIGNIFICANCE != "VUS") |>
        dplyr::distinct()
    }
  }

  return(list('callset' = callset,
              'panel_info' = panel_info,
              'sample_id' = sample_id,
              eval = TRUE))


}

#' Function that reads and validates CNA or SNV/InDel TSV files
#' file from PCGR/CPSR pre-report (Python) pipeline
#'
#' @param fname Path to raw input file with DNA aberrations (PCGR/CPSR)
#' @param cols column type definitions of raw input file
#' @param ref_data PCGR reference data object
#' @param settings PCGR run/configuration settings
#' @param vartype type of DNA aberrations ('snv_indel','cna')
#' @param primary_site primary site of tumor
#' @param retained_info_tags VCF INFO tags to be retained in output (SNVs/InDels)
#' @param variant_origin Germline/Somatic
#'
#' @export
#'
load_dna_variants <- function(
    fname = NA,
    cols = NULL,
    ref_data = NULL,
    settings = NULL,
    vartype = 'snv_indel',
    primary_site = "Any",
    retained_info_tags = "None",
    variant_origin = "Somatic") {

  ## assert that ref_data is non-null, settings is non-null
  ## and that "molecular_data" is an element of the settings list
  if (is.null(ref_data)) {
    log4r_fatal(
      "Reference data object is NULL - cannot proceed with loading of DNA variants")
  }
  if (is.null(settings)) {
    log4r_fatal(
      "Settings object is NULL - cannot proceed with loading of DNA variants")
  }
  if ("molecular_data" %in% names(settings) == FALSE) {
    log4r_fatal(
      "Settings object is lacking a 'molecular_data' section - cannot proceed with loading of
      DNA variants")
  }

  check_file_exists(fname)
  calls_raw <- suppressWarnings(
    as.data.frame(
      readr::read_tsv(
        file = fname,
        na = c(".","NA"),
        show_col_types = F,
        progress = F
      )
    )
  )

  ## check that all columns are present among columns
  ## read from file
  compulsary_cols <-
    names(cols$cols)

  raw_col_check <-
    rlang::has_name(calls_raw, compulsary_cols)
  if (FALSE %in% raw_col_check) {
    missing_cols <-
      compulsary_cols[!raw_col_check]
    log4r_fatal(
      paste0("Missing required columns in input file ",
             fname, " - ",
             paste(missing_cols, collapse=", ")))
  }

  cols_including_retained <- cols
  retained_cols <- NULL
  calls <- calls_raw
  if (retained_info_tags != "None") {
    retained_cols <- stringr::str_split(
      retained_info_tags, pattern = ",")[[1]]
    for (c in retained_cols) {
      if (c %in% colnames(calls_raw)) {
        col_retain <- readr::cols_only(
          !!rlang::sym(c) := readr::col_character()
        )
        cols_including_retained$cols <-
          c(cols_including_retained$cols, col_retain$cols)
      }
    }

  }

  calls <- suppressWarnings(
    as.data.frame(
      readr::read_tsv(
        file = fname,
        col_types =
          cols_including_retained,
        na = c(".","NA"),
        progress = F
      )
    )
  )

  retained_cols_renamed <- c()
  if (!is.null(retained_cols)) {
    for (c in retained_cols) {
      if (c %in% colnames(calls)) {
        new_col <- paste0('VCF_INFO_', c)
        retained_cols_renamed <- c(
          retained_cols_renamed, new_col
        )
        calls[,new_col] <- calls[,c]
        calls[,c] <- NULL
      }
    }
  }

  results <- init_var_content()
  results[['variant']] <- calls

  ## Rename specific columns/annotations for more clarity
  colrename_map <- c("TSG" = "TUMOR_SUPPRESSOR",
                     "HGVSp_short" = "HGVSP",
                     "TSG_RANK" = "TUMOR_SUPPRESSOR_RANK",
                     "TSG_SUPPORT" = "TUMOR_SUPPRESSOR_SUPPORT",
                     "SPLICE_EFFECT_MUTSPLICEDB" = "SPLICE_EFFECT")

  for (old_name in names(colrename_map)) {
    new_name <- colrename_map[[old_name]]
    if (old_name %in% colnames(results[['variant']])) {
      results[['variant']] <-
        results[['variant']] |>
        dplyr::rename(
          !!rlang::sym(new_name) := !!rlang::sym(old_name)
        )
    }
  }

  ## Combine splice effect predictions
  ## from MutSpliceDB and MaxEntScan
  if ("MAXENTSCAN" %in% colnames(results[['variant']]) &
      "SPLICE_EFFECT" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(
        SPLICE_EFFECT = paste(
          .data$SPLICE_EFFECT,
          .data$MAXENTSCAN, sep = ", "
      )) |>
      dplyr::mutate(
        SPLICE_EFFECT = stringr::str_replace_all(
          .data$SPLICE_EFFECT,
          "(^, )|^(NA, )|^(NA, NA)$", ""
        )
      ) |>
      dplyr::mutate(
        SPLICE_EFFECT = dplyr::if_else(
          .data$SPLICE_EFFECT == "" |
            .data$SPLICE_EFFECT == "NA",
          NA_character_,
          as.character(.data$SPLICE_EFFECT)
        )
      )
  }

  if (vartype == "snv_indel") {
    results[['variant']] <- clean_gnomad_annotations(
      var_df = results[['variant']]
    )
  }

  ## Re-format annotations
  if ("VEP_ALL_CSQ" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(
        VEP_ALL_CSQ = stringr::str_replace_all(
          .data$VEP_ALL_CSQ, ",",", "
        )
      )
  }

  for (col in c('VAF_TUMOR','VAF_CONTROL','TPM','GERP_SCORE',
               'TPM_GENE','TPM_MIN','consTPM')) {
    if (col %in% colnames(results[['variant']])) {
      significant_digits = 3
      if (col == 'TPM' | col == 'TPM_GENE' |
         col == 'TPM_MIN' | col == 'consTPM') {
        significant_digits = 2
      }
      results[['variant']] <-
        results[['variant']] |>
        dplyr::mutate(
          !!col := round(as.numeric(.data[[col]]), significant_digits)
        )
    }
  }

  ## Still some bugs/artefacts in TPM assignment to genes
  ## (some genes lack 'TPM_GENE', but have 'TPM' values for transcripts?)
  ## that needs to be fixed - temporary clean up
  if ("TPM_GENE" %in% colnames(results[['variant']]) &
       "TPM" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(
        TPM_GENE = dplyr::case_when(
          is.na(.data$TPM_GENE) & !is.na(.data$TPM) ~ .data$TPM,
          is.na(.data$TPM_GENE) & is.na(.data$TPM) ~ 0,
          TRUE ~ .data$TPM_GENE
        )
      )
  }

  if ("TPM" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(
        TPM = dplyr::case_when(
          is.na(.data$TPM) ~ 0,
          TRUE ~ .data$TPM
        )
      )
  }

  if ("consTPM" %in% colnames(results[['variant']]) &
      "TPM" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(
        consTPM = dplyr::case_when(
          is.na(.data$consTPM) ~ .data$TPM,
          TRUE ~ .data$consTPM
        )
      )
  }

  if (vartype == 'cna') {

    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(REFSEQ_TRANSCRIPT_ID = dplyr::if_else(
        is.na(.data$REFSEQ_TRANSCRIPT_ID),
        "",
        as.character(.data$REFSEQ_TRANSCRIPT_ID)
      )) |>
      dplyr::mutate(
        TRANSCRIPT_OVERLAP = paste(
          .data$ENSEMBL_TRANSCRIPT_ID,
          .data$REFSEQ_TRANSCRIPT_ID,
          .data$TRANSCRIPT_START,
          .data$TRANSCRIPT_END,
          .data$TRANSCRIPT_OVERLAP_PERCENT, sep="|"
        )) |>
      dplyr::select(
        -dplyr::any_of(
          c("ENSEMBL_TRANSCRIPT_ID",
           "REFSEQ_TRANSCRIPT_ID",
           "TPM",
           "TRANSCRIPT_START",
           "TRANSCRIPT_END")
          )
        ) |>
      dplyr::group_by(
        dplyr::across(-c("TRANSCRIPT_OVERLAP",
                         "TRANSCRIPT_OVERLAP_PERCENT"))) |>
      dplyr::summarise(
        TRANSCRIPT_OVERLAP = paste(
          .data$TRANSCRIPT_OVERLAP, collapse=", "),
        MAX_TRANSCRIPT_OVERLAP_PERCENT =
          max(.data$TRANSCRIPT_OVERLAP_PERCENT,
              na.rm = T),
        .groups = "drop"
      ) |>
      dplyr::distinct()

  }

  for(c in c("ONCOGENICITY","CLINVAR_CLASSIFICATION")) {
    if (c %in% colnames(results[['variant']])) {
      results[['variant']] <-
        results[['variant']] |>
        dplyr::mutate(
          !!c := stringr::str_replace_all(
            .data[[c]], "_", " "
          )
        )
    }
  }

  ## Germline - CPSR, use CIViC for now
  if(variant_origin == "Germline"){
    results[['bm_evidence']][['eitems']] <-
      map_biomarker_data(
        varcalls = results[['variant']],
        ref_data = ref_data,
        variant_origin = variant_origin,
        vartype = vartype
      )
  }


  if(variant_origin == "Somatic"){

    ## A) no OncoKB token provided, or
    ## B) OncoKB token provided and OncoKB annotation
    ##    is to be run non-exclusively
    ## --> map also CGI/CIViC markers
    if("oncokb" %in% names(settings$conf) &
       (settings$conf$oncokb$run == 0 |
        (settings$conf$oncokb$run == 1 &
         settings$conf$oncokb$exclusive == 0))){
      results[['bm_evidence']][['eitems']] <-
        map_biomarker_data(
          varcalls = results[['variant']],
          ref_data = ref_data,
          variant_origin = variant_origin,
          vartype = vartype
        )

    }

    if(settings$conf$oncokb$run == 1 &
       vartype == "snv_indel"){
      if(!is.null(settings$molecular_data$fname_oncokb_output_maf_hgvsp) &
         file.exists(settings$molecular_data$fname_oncokb_output_maf_hgvsp)){

        log4r_info(
          "Querying OncoKB web API with OncoKB-annotated HGVSp/HGVSg MAF files")

        oncokb_results <-
          process_oncokb_maf(
            maf_file_hgvsp =
              settings$molecular_data$fname_oncokb_output_maf_hgvsp,
            maf_file_hgvsg =
              settings$molecular_data$fname_oncokb_output_maf_hgvsg,
            oncokb_token =
              settings$conf$oncokb$api_token,
            oncotree_code =
              settings$conf$oncokb$oncotree_code,
            var_calls = results[['variant']])

        if(NROW(oncokb_results$eitems) > 0){
          results$bm_evidence$eitems <- dplyr::bind_rows(
            results$bm_evidence$eitems,
            oncokb_results$eitems)
        }
      }
    }

    if(settings$conf$oncokb$run == 1 &
       vartype == "cna"){
      if(!is.null(settings$molecular_data$fname_oncokb_output_cna) &
         file.exists(settings$molecular_data$fname_oncokb_output_cna)){

        log4r_info(
          "Querying OncoKB web API with OncoKB-annotated CNA file")

        oncokb_results <-
          process_oncokb_cna(
            cna_file =
              settings$molecular_data$fname_oncokb_output_cna,
            oncokb_token =
              settings$conf$oncokb$api_token,
            oncotree_code =
              settings$conf$oncokb$oncotree_code,
            var_calls = results[['variant']])

        if(NROW(oncokb_results$eitems) > 0){
          results$bm_evidence$eitems <- dplyr::bind_rows(
            results$bm_evidence$eitems,
            oncokb_results$eitems)
        }
      }
    }
  }

  ## Assign each variant a tier according to
  ## AMP/ASCO/CAP guidelines for clinical actionability
  ## of somatic variants
  if (variant_origin == "Somatic") {

    if (vartype != "fusion") {
      log4r_info(
        paste0("Variant tier classification",
               " - somatic actionability guidelines (AMP/ASCO/CAP)"))

      amp_asco_cap_classified_variant <- list()
      for (clnsig in names(bm_categories)) {
        amp_asco_cap_classified_variant[[clnsig]] <-
          assign_amp_asco_cap_tiers(
            vartype = vartype,
            var_df = results$variant,
            primary_site = primary_site,
            clinical_significance = clnsig,
            biomarker_mapping_confidence = "medium",
            biomarker_items =
              results$bm_evidence$eitems
          )
      }

      for (clnsig in names(bm_categories)) {
        for (elem in c("classification",
                       "eitems")) {
          results$bm_evidence[[clnsig]][[elem]] <-
            amp_asco_cap_classified_variant[[clnsig]]$bm_evidence[[elem]]
        }
        if (clnsig == "therapeutic_sensitivity") {
          if (NROW(amp_asco_cap_classified_variant[[clnsig]][['variant']]) > 0) {
            results[['variant']] <-
              amp_asco_cap_classified_variant[[clnsig]][['variant']]
            results[['bm_evidence']][['classification']] <-
              amp_asco_cap_classified_variant[[clnsig]][['classification']]
          }
          ## issue warning if no variants are classified
          if (NROW(amp_asco_cap_classified_variant[[clnsig]][['variant']]) == 0) {
            log4r_warn(
              paste0(
                "No variants were classified with AMP/ASCO/CAP tiers of clinical",
                "significance for therapeutic sensitivity"))
          }
        }
      }
    }

  }

  results[['retained_info_tags']] <- paste(
    retained_cols_renamed, collapse=","
  )

  return(results)

}

#' Load expression similarity results
#'
#' @param settings PCGR run/configuration settings
#'
#' @export
load_expression_similarity <- function(settings = NULL) {

  ## Load expression similarity results for input sample
  ## against reference collections
  expression_sim <- list()

  if (settings$conf$expression$similarity_analysis == 0) {
    return(expression_sim)
  }

  if (file.exists(
    settings$molecular_data$fname_expression_similarity_tsv)) {
    log4r_info(
      paste0("Loading expression similarity results for sample ",
             settings$sample_id))

    expression_similarity <- suppressWarnings(readr::read_tsv(
      settings$molecular_data$fname_expression_similarity_tsv,
      show_col_types = F,
      na = ".", guess_max = 100000))

    for (source in unique(expression_similarity$EXT_DB)) {
      expression_sim[[source]] <- expression_similarity |>
        dplyr::filter(.data$EXT_DB == source) |>
        dplyr::select(-c("EXT_DB"))
    }
  }


  return(expression_sim)

}

#' Load expression outlier results
#'
#' @param settings PCGR run/configuration settings
#' @param ref_data PCGR reference data object
#' @param percentile_cutoff_high numeric, percentile cutoff for high expression
#' @param percentile_cutoff_low numeric, percentile cutoff for low expression
#' @param z_score_cutoff numeric, z-score cutoff for expression outliers
#'
#' @export
load_expression_outliers <- function(
    settings = NULL,
    ref_data = NULL,
    percentile_cutoff_high = 95,
    percentile_cutoff_low = 5,
    z_score_cutoff = 1.5) {

  ## Load expression outlier results for input sample
  ## against reference collections

  tumor_site <- settings$conf$sample_properties$site
  expression_outliers <- data.frame()
  if (file.exists(
    settings$molecular_data$fname_expression_outliers_tsv)) {
    log4r_info(
      paste0("Loading expression outlier results for sample ",
             settings$sample_id))

    ## Read raw expression outlier data from Python step of PCGR
    outlier_data <- suppressWarnings(readr::read_tsv(
      settings$molecular_data$fname_expression_outliers_tsv,
      show_col_types = F,
      na = ".", guess_max = 100000))

    if (!is.null(ref_data) &
       !is.null(ref_data$gene) &
       !is.null(ref_data$gene$gene_xref) &
       NROW(outlier_data) > 0 &
       "ENSEMBL_GENE_ID" %in% colnames(outlier_data)) {
      outlier_data <-
        outlier_data |>
        dplyr::filter(!is.na(.data$ENSEMBL_GENE_ID))

      ## Append gene annotations from reference data
      ## (SYMBOL, ENTREZGENE, GENENAME, TSG, GENE_BIOTYPE, ONCOGENE)
      if (NROW(outlier_data) > 0) {
        outlier_data <- outlier_data |>
          dplyr::inner_join(
            dplyr::select(
              ref_data$gene$gene_xref,
              c("SYMBOL","ENSEMBL_GENE_ID",
                "ENTREZGENE","GENENAME",
                "TSG","GENE_BIOTYPE",
                "ONCOGENE")),
            by = c("ENSEMBL_GENE_ID")
          )


        if (NROW(outlier_data) == 0) {
          return(expression_outliers)
        }

        ## Rename annotations for more clarity
        if ("TSG" %in% colnames(outlier_data)) {
          outlier_data <-
            outlier_data |>
            dplyr::rename(
              TUMOR_SUPPRESSOR = "TSG"
            )
        }

        ## Add links for display in datatable of HTML report
        expression_outliers <- outlier_data |>
          append_cancer_gene_evidence(
            ref_data = ref_data) |>
          append_cancer_association_ranks(
            ref_data = ref_data,
            primary_site = tumor_site) |>
          dplyr::mutate(VAR_ID = dplyr::row_number()) |> ## add unique ID
          append_drug_var_link(
            primary_site = tumor_site,
            ref_data = ref_data) |>
          append_annotation_links(
            vartype = "exp",
            skip = c("DBSNP_RSID",
                     "CLINVAR",
                     "PROTEIN_DOMAIN",
                     "COSMIC_ID",
                     "REFSEQ_TRANSCRIPT_ID",
                     "ENSEMBL_TRANSCRIPT_ID")) |>
          dplyr::select(
            -dplyr::contains("_RAW")) |>
          ## define criteria for expression outliers
          dplyr::mutate(EXPR_OUTLIER = dplyr::case_when(
            PERCENTILE >= percentile_cutoff_high &
              .data$IQR > 0 &
              TPM_LOG2_GENE > (.data$Q3 + (1.5 * .data$IQR)) &
              abs(.data$Z_SCORE) > z_score_cutoff ~ "Increased expression",
            PERCENTILE <= percentile_cutoff_low &
              .data$IQR > 0 &
              TPM_LOG2_GENE < (.data$Q1 - (1.5 * .data$IQR)) &
              abs(.data$Z_SCORE) > z_score_cutoff ~ "Reduced expression",
            TRUE ~ as.character(NA)
          )) |>
          dplyr::filter(
              .data$GENE_BIOTYPE == "protein_coding") |>
          dplyr::mutate(Z_SCORE = round(.data$Z_SCORE, 1))

        if (NROW(expression_outliers) == 0) {
          return(expression_outliers)
        }

        expression_outliers <- expression_outliers |>
          dplyr::filter(!is.na(.data$EXPR_OUTLIER))

        if (NROW(expression_outliers) == 0) {
          return(expression_outliers)
        }

        expression_outliers <- expression_outliers |>
          dplyr::filter(
            .data$GLOBAL_ASSOC_RANK > 0 |
              .data$TISSUE_ASSOC_RANK > 0 |
              .data$TUMOR_SUPPRESSOR == TRUE |
              .data$ONCOGENE == TRUE |
              !is.na(.data$CANCERGENE_EVIDENCE) |
              !is.na(.data$TARGETED_INHIBITORS_ALL2)
          )

        if (NROW(expression_outliers) == 0) {
          return(expression_outliers)
        }

        expression_outliers <- expression_outliers |>
          dplyr::arrange(
            dplyr::desc(abs(.data$Z_SCORE)),
            dplyr::desc(.data$GLOBAL_ASSOC_RANK),
            dplyr::desc(.data$TISSUE_ASSOC_RANK)) |>
          dplyr::mutate(EXPR_LEVEL = round(
            .data$TPM_LOG2_GENE, digits = 3)) |>
          dplyr::rename(REF_COHORT_IQR = "IQR") |>
          dplyr::mutate(REF_COHORT_QUARTILES = paste(
            paste0("Q1: ",.data$Q1),
            paste0("Q2: ",.data$Q2),
            paste0("Q3: ",.data$Q3),
            sep = " | ")) |>
          dplyr::select(
            c("SYMBOL",
              "GENENAME",
              "EXPR_LEVEL",
              "EXPR_OUTLIER",
              "Z_SCORE",
              "REF_COHORT",
              "REF_COHORT_SIZE",
              "REF_COHORT_IQR",
              "REF_COHORT_QUARTILES",
              "PERCENTILE",
              "TUMOR_SUPPRESSOR",
              "ONCOGENE",
              "TARGETED_INHIBITORS_ALL",
              "ENSEMBL_GENE_ID",
              "GENE_BIOTYPE",
              "GLOBAL_ASSOC_RANK",
              "TISSUE_ASSOC_RANK",
              "CANCERGENE_EVIDENCE")) |>
          dplyr::mutate(REF_COHORT = toupper(
            stringr::str_replace_all(
              .data$REF_COHORT,"_","-")
          ))

      }

    }
  }

  return(expression_outliers)

}

#' Load expression consequence settings
#'
#' @param settings PCGR run/configuration settings
#'
#' @export
load_expression_csq <- function(settings = NULL) {

  ## Load expression consequence results for input sample
  ## against reference collections
  expression_csq <- data.frame()

  if (!is.null(settings$conf$expression$consequence_db)) {

    if (!is.null(settings$molecular_data[['fname_csq_expression_tsv']]) &
       file.exists(settings$molecular_data[['fname_csq_expression_tsv']])) {
      expression_csq <- readr::read_tsv(
        settings$molecular_data[['fname_csq_expression_tsv']],
        show_col_types = F, na = ".")
    }
  }

  return(expression_csq)

}

#' Load RNA fusion results
#'
#' @param settings PCGR run/configuration settings
#' @param ref_data PCGR reference data object
#'
#' @export
load_rna_fusions <- function(
    settings = NULL,
    ref_data = NULL) {


  log4r_info("------")
  log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - somatic RNA fusions"))

  primary_site <-
    settings$conf$sample_properties$site

  hgname <- "hg38"
  if (settings$genome_assembly == "grch37") {
    hgname <- "hg19"
  }

  results <- list()
  results[['variant']] <- data.frame()
  results[['variant_recurrence']] <- data.frame()
  results[['variant_display']] <- data.frame()
  results[['bm_evidence']] <- list()
  results[['bm_evidence']][['eitems']] <-
    data.frame()
  results[['bm_evidence']][['classification']] <-
    data.frame()


  rna_fusions <- data.frame()
  #rna_fusions_mapped <- data.frame()

  rna_fusion_raw_fname <-
    settings$molecular_data$fname_rna_fusion_tsv

  check_file_exists(rna_fusion_raw_fname)
  rna_fusion_calls_raw <- suppressWarnings(
    as.data.frame(
      readr::read_tsv(
        file = rna_fusion_raw_fname,
        guess_max = 0,
        na = c(".","NA"),
        show_col_types = F,
        progress = F
      )
    )
  )

  ## check that all columns are present among columns
  ## read from file
  compulsary_cols <-
    names(data_coltype_defs$rna_fusion_raw$cols)

  raw_col_check <-
    rlang::has_name(rna_fusion_calls_raw, compulsary_cols)
  if (FALSE %in% raw_col_check) {
    missing_cols <-
      compulsary_cols[!raw_col_check]
    log4r_fatal(
      paste0("Missing required columns in input file ",
             rna_fusion_raw_fname, " - ",
             paste(missing_cols, collapse=", ")))
  }

  cols_including_optional <- data_coltype_defs$rna_fusion_raw
  if ("FUSION_SCORE" %in% colnames(rna_fusion_calls_raw)) {
    score_col <- readr::cols_only(
      FUSION_SCORE = readr::col_number())
    cols_including_optional$cols <-
      c(cols_including_optional$cols, score_col$cols)
  }


  callset_fusions <- load_dna_variants(
    fname = rna_fusion_raw_fname,
    cols = cols_including_optional,
    ref_data = ref_data,
    settings = settings,
    vartype = 'fusion',
    primary_site =
      primary_site,
    retained_info_tags = "None",
    variant_origin = "Somatic")

  results$variant <- callset_fusions$variant
  results$bm_evidence <-
    callset_fusions$bm_evidence

  rna_fusion_calls <- callset_fusions$variant

  if (NROW(rna_fusion_calls) == 0) {
    log4r_info(
      paste0("No RNA fusion calls found in input file ",
             rna_fusion_raw_fname))
    return(rna_fusions)
  }

  log4r_info(
    paste0("Loading RNA fusion data for sample ",
           settings$sample_id))

  rna_fusions <- rna_fusion_calls |>
    dplyr::mutate(
      FUSION_GENE2 = stringr::str_replace_all(
        .data$FUSION_GENE, "\\-\\-", "::"
      )
    ) |>
    tidyr::separate(
      "BREAKPOINT_5P",
      into = c("BP_CHROM_5P","BP_POSITION_5P"),
      sep = ":",
      remove = FALSE
    ) |>
    tidyr::separate(
      "BREAKPOINT_3P",
      into = c("BP_CHROM_3P","BP_POSITION_3P"),
      sep = ":",
      remove = FALSE
    ) |>
    dplyr::mutate(
      BP_POSITION_5P =
        as.integer(.data$BP_POSITION_5P),
      BP_POSITION_3P =
        as.integer(.data$BP_POSITION_3P)
    ) |>
    tidyr::separate(
      "FUSION_GENE",
      into = c("FUSION_GENE_5P",
               "FUSION_GENE_3P"),
      sep = "--",
      remove = FALSE
    )

  score_splitread_df <- data.frame()
  if (NROW(rna_fusions) > 0) {
    score_splitread_df <- rna_fusions |>
      dplyr::select(
        dplyr::any_of(
          c("SAMPLE_ID",
            "VAR_ID",
            "FUSION_SCORE",
            "SPLIT_READS")
        )
      )
  }

  # Get transcript overlap and rename suffix (5P or 3P)
  get_bp_transcript_overlap <- function(df, prefix, ref_data) {
    bp_cols <- paste0(c("BP_CHROM_", "BP_POSITION_"), prefix)
    bp_junctions <- df |>
      dplyr::select(dplyr::all_of(bp_cols)) |>
      dplyr::rename(BP_CHROM = dplyr::all_of(bp_cols[1]),
                    BP_POSITION = dplyr::all_of(bp_cols[2])) |>
      dplyr::filter(!is.na(.data$BP_CHROM) & !is.na(.data$BP_POSITION))

    bp_junction_transcript_overlap(bp_junctions, ref_data) |>
      dplyr::rename_with(~ paste0(.x, "_", prefix))
  }

  # refgene table
  # Get reference gene table for 5P or 3P
  ref_gene_table <- function(ref_data, suffix) {
    ref_data$gene$gene_xref |>
      dplyr::select(c("SYMBOL", "ENTREZGENE",
                    "ENSEMBL_GENE_ID",
                    "ONCOGENE", "TSG", "GENENAME")) |>
      dplyr::rename_with(~ paste0(.x, "_", suffix)) |>
      dplyr::distinct()
  }

  transcript_overlap_5p <- get_bp_transcript_overlap(
    df = rna_fusions, prefix = "5P", ref_data = ref_data)
  transcript_overlap_3p <- get_bp_transcript_overlap(
    df = rna_fusions, prefix = "3P", ref_data = ref_data)

  ref_data_gene <- list()
  ref_data_gene[['3P']] <- ref_gene_table(ref_data, "3P")
  ref_data_gene[['5P']] <- ref_gene_table(ref_data, "5P")

  rna_fusions <- rna_fusions |>
    dplyr::left_join(
      transcript_overlap_5p,
      by = c("BP_CHROM_5P" = "BP_CHROM_5P",
             "BP_POSITION_5P" = "BP_POSITION_5P"),
      relationship = "many-to-many"
    ) |>
    dplyr::left_join(
      transcript_overlap_3p,
      by = c("BP_CHROM_3P" = "BP_CHROM_3P",
             "BP_POSITION_3P" = "BP_POSITION_3P"),
      relationship = "many-to-many"
    ) |>
    dplyr::distinct()

  ## Join with reference gene table to get
  ## gene annotations (ENTREZGENE, ONCOGENE, TSG,
  ## GENENAME) for 5' and 3' fusion partner
  rna_fusions <- rna_fusions |>
    dplyr::left_join(
      ref_data_gene[['5P']], by = "ENSEMBL_GENE_ID_5P"
    ) |>
    dplyr::left_join(
      ref_data_gene[['3P']], by = "ENSEMBL_GENE_ID_3P"
    ) |>
    ## reveal proper 5' and 3' fusion partners if not
    ## provided directly in fusion call (e.g. "--ALK", "FGFR2--")
    dplyr::mutate(
      FUSION_GENE_5P = dplyr::if_else(
        is.na(.data$FUSION_GENE_5P) |
          .data$FUSION_GENE_5P == "",
        .data$SYMBOL_5P,
        .data$FUSION_GENE_5P
      ),
      FUSION_GENE_3P = dplyr::if_else(
        is.na(.data$FUSION_GENE_3P) |
          .data$FUSION_GENE_3P == "",
        .data$SYMBOL_3P,
        .data$FUSION_GENE_3P
      )
    ) |>
    dplyr::mutate(
      BP_TRANSCRIPT_MATCH_SYMBOL_3P = dplyr::if_else(
        .data$FUSION_GENE_3P == .data$SYMBOL_3P,
        TRUE, FALSE, FALSE),
      BP_TRANSCRIPT_MATCH_SYMBOL_5P = dplyr::if_else(
        .data$FUSION_GENE_5P == .data$SYMBOL_5P,
        TRUE, FALSE, FALSE)
    ) |>
    dplyr::distinct()


  ## Split fusion calls into those with
  ## breakpoint-transcript matches
  ## and those without breakpoint-transcript matches
  ## for 5' and 3' fusion partner gene
  transcript_matches <- list()
  transcript_non_matches <- list()
  transcript_matches[["5P"]] <-
    rna_fusions |>
    dplyr::filter(
      .data$BP_TRANSCRIPT_MATCH_SYMBOL_5P == TRUE)

  if (NROW(transcript_matches[["5P"]]) > 0) {
    transcript_matches[["5P"]] <-
      transcript_matches[["5P"]] |>
      dplyr::select(
        c("SAMPLE_ID",
          "VAR_ID",
          "VARIANT_CLASS",
          "FUSION_GENE",
          "FUSION_GENE2",
          "BREAKPOINT_5P",
          "BREAKPOINT_3P",
          "FUSION_GENE_5P",
          "ENSEMBL_GENE_ID_5P",
          "ENSEMBL_TRANSCRIPT_ID_5P",
          "ENTREZGENE_5P",
          "GENENAME_5P",
          "ONCOGENE_5P",
          "TSG_5P")) |>
      dplyr::group_by(
        dplyr::across(
          -c("ENSEMBL_TRANSCRIPT_ID_5P")
        )
      ) |>
      dplyr::summarise(
        ENSEMBL_TRANSCRIPT_ID_5P = paste(
          unique(.data$ENSEMBL_TRANSCRIPT_ID_5P),
          collapse = "; "),
        .groups = "drop"
      ) |>
      dplyr::distinct()
  }

  transcript_non_matches[['5P']] <-
    rna_fusions |>
    dplyr::filter(
      .data$BP_TRANSCRIPT_MATCH_SYMBOL_5P == FALSE)


  if (NROW(transcript_non_matches[['5P']]) > 0 &
     NROW(transcript_matches[['5P']]) > 0) {
    transcript_non_matches[['5P']] <-
      transcript_non_matches[['5P']] |>
      dplyr::anti_join(
        transcript_matches[['5P']],
        by = c("SAMPLE_ID",
               "VAR_ID",
               "VARIANT_CLASS",
               "FUSION_GENE",
               "FUSION_GENE_5P")
      )

  }

  if (NROW(transcript_non_matches[['5P']]) > 0) {

    gene_list_5p <- unique(transcript_non_matches[['5P']]$FUSION_GENE_5P)
    n_genes_5p <- length(gene_list_5p)
    gene_display_5p <- if (n_genes_5p <= 3) {
      paste0("'", paste(gene_list_5p, collapse = "', '"), "'")
    } else {
      paste0("'", paste(gene_list_5p[1:3], collapse = "', '"),
             "' (and ", n_genes_5p - 3, " more)")
    }

    log4r_warn(
      paste0("Breakpoints of 5' fusion partner genes do not intersect ",
             "database transcript annotations"))
    log4r_warn(
      paste0(gene_display_5p,
             " - check gene symbols or breakpoint coordinates."))

    transcript_non_matches[['5P']] <-
      transcript_non_matches[['5P']] |>
      dplyr::select(
        c("SAMPLE_ID",
          "VAR_ID",
          "VARIANT_CLASS",
          "FUSION_GENE",
          "FUSION_GENE2",
          "BREAKPOINT_5P",
          "BREAKPOINT_3P",
          "FUSION_GENE_5P"
        )) |>
      dplyr::distinct() |>
      dplyr::left_join(
        ref_data$gene$gene_xref |>
          dplyr::select(
            c("SYMBOL",
              "ENTREZGENE",
              "ENSEMBL_GENE_ID",
              "ONCOGENE",
              "TSG",
              "GENENAME")),
        by = c("FUSION_GENE_5P" = "SYMBOL"),
        relationship = "many-to-many"
      ) |>
      dplyr::rename(
        ENTREZGENE_5P = "ENTREZGENE",
        ENSEMBL_GENE_ID_5P = "ENSEMBL_GENE_ID",
        GENENAME_5P = "GENENAME",
        ONCOGENE_5P = "ONCOGENE",
        TSG_5P = "TSG"
      ) |>
      dplyr::filter(
        !is.na(.data$ENTREZGENE_5P))

    # if (NROW(transcript_non_matches[['5P']]) == 0) {
    #   log4r_warn(
    #     paste0(
    #       "No 5' fusion partner genes could be mapped to official gene symbols ",
    #       "in reference database - check gene symbols or breakpoint coordinates."))
    #   transcript_non_matches[['5P']] <-
    #     data.frame()
    # }
  }else{
    transcript_non_matches[['5P']] <-
      data.frame()
  }

  transcript_matches[["3P"]] <-
    rna_fusions |>
    dplyr::filter(
      .data$BP_TRANSCRIPT_MATCH_SYMBOL_3P == TRUE)

  if (NROW(transcript_matches[["3P"]]) > 0) {
    transcript_matches[["3P"]] <-
      transcript_matches[["3P"]] |>
      dplyr::select(
        c("SAMPLE_ID",
          "VAR_ID",
          "VARIANT_CLASS",
          "FUSION_GENE",
          "FUSION_GENE2",
          "BREAKPOINT_5P",
          "BREAKPOINT_3P",
          "FUSION_GENE_3P",
          "ENSEMBL_GENE_ID_3P",
          "ENSEMBL_TRANSCRIPT_ID_3P",
          "ENTREZGENE_3P",
          "GENENAME_3P",
          "ONCOGENE_3P",
          "TSG_3P")) |>
      dplyr::group_by(
        dplyr::across(
          -c("ENSEMBL_TRANSCRIPT_ID_3P")
        )
      ) |>
      dplyr::summarise(
        ENSEMBL_TRANSCRIPT_ID_3P = paste(
          unique(.data$ENSEMBL_TRANSCRIPT_ID_3P),
          collapse = "; "),
        .groups = "drop"
      ) |>
      dplyr::distinct()
  }

  transcript_non_matches[['3P']] <-
    rna_fusions |>
    dplyr::filter(
      .data$BP_TRANSCRIPT_MATCH_SYMBOL_3P == FALSE)

  if (NROW(transcript_non_matches[['3P']]) > 0 &
     NROW(transcript_matches[['3P']]) > 0) {
    transcript_non_matches[['3P']] <-
      transcript_non_matches[['3P']] |>
      dplyr::anti_join(
        transcript_matches[['3P']],
        by = c("SAMPLE_ID",
               "VAR_ID",
               "VARIANT_CLASS",
               "FUSION_GENE",
               "FUSION_GENE_3P")
      )
  }

  if (NROW(transcript_non_matches[['3P']]) > 0) {

    gene_list_3p <- unique(transcript_non_matches[['3P']]$FUSION_GENE_3P)
    n_genes_3p <- length(gene_list_3p)
    gene_display_3p <- if (n_genes_3p <= 3) {
      paste0("'", paste(gene_list_3p, collapse = "', '"), "'")
    } else {
      paste0("'", paste(gene_list_3p[1:3], collapse = "', '"),
             "' (and ", n_genes_3p - 3, " more)")
    }

    log4r_warn(
      paste0("Breakpoints of 3' fusion partner genes do not intersect ",
             "database transcript annotations"))
    log4r_warn(
      paste0(gene_display_3p,
             " - check gene symbols or breakpoint coordinates."))

    transcript_non_matches[['3P']] <-
      transcript_non_matches[['3P']] |>
      dplyr::select(
        c("SAMPLE_ID",
          "VAR_ID",
          "VARIANT_CLASS",
          "FUSION_GENE",
          "FUSION_GENE2",
          "BREAKPOINT_5P",
          "BREAKPOINT_3P",
          "FUSION_GENE_3P"
        )) |>
      dplyr::distinct() |>
      dplyr::left_join(
        ref_data$gene$gene_xref |>
          dplyr::select(
            c("SYMBOL",
              "ENTREZGENE",
              "ENSEMBL_GENE_ID",
              "ONCOGENE",
              "TSG",
              "GENENAME")),
        by = c("FUSION_GENE_3P" = "SYMBOL"),
        relationship = "many-to-many"
      ) |>
      dplyr::rename(
        ENTREZGENE_3P = "ENTREZGENE",
        ENSEMBL_GENE_ID_3P = "ENSEMBL_GENE_ID",
        GENENAME_3P = "GENENAME",
        ONCOGENE_3P = "ONCOGENE",
        TSG_3P = "TSG"
      ) |>
      dplyr::filter(
        !is.na(.data$ENTREZGENE_3P))

    # if (NROW(transcript_non_matches[['3P']]) == 0) {
    #   log4r_warn(
    #     paste0("No 3' fusion partner genes could be mapped to official gene symbols ",
    #            "in reference database - check gene symbols or breakpoint coordinates."))
    #   transcript_non_matches[['3P']] <-
    #     data.frame()
    # }
  }else{
    transcript_non_matches[['3P']] <-
      data.frame()
  }

  ## Combine fusion calls with and without
  ## breakpoint-transcript matches
  ## for 5' and 3' fusion partner gene

  if (NROW(transcript_matches[['5P']]) > 0 |
     NROW(transcript_non_matches[['5P']]) > 0 |
     NROW(transcript_matches[['3P']]) > 0 |
     NROW(transcript_non_matches[['3P']]) > 0) {

    transcript_non_matches_5P <- data.frame()
    if (NROW(transcript_non_matches[['5P']]) > 0) {
      transcript_non_matches_5P <-
      dplyr::select(
        transcript_non_matches[['5P']],
        c("SAMPLE_ID",
          "VAR_ID",
          "VARIANT_CLASS",
          "FUSION_GENE",
          "FUSION_GENE2",
          "BREAKPOINT_3P"),
        dplyr::ends_with("5P"))
    }

    transcript_non_matches_3P <- data.frame()
    if (NROW(transcript_non_matches[['3P']]) > 0) {
      transcript_non_matches_3P <-
        dplyr::select(
          transcript_non_matches[['3P']],
          c("SAMPLE_ID",
            "VAR_ID",
            "VARIANT_CLASS",
            "FUSION_GENE",
            "FUSION_GENE2",
            "BREAKPOINT_5P"),
          dplyr::ends_with("3P"))
    }

    results[['variant']] <-
      dplyr::bind_rows(
        transcript_matches[['5P']],
        transcript_non_matches_5P
      ) |>
      dplyr::left_join(
        dplyr::bind_rows(
          transcript_matches[['3P']],
          transcript_non_matches_3P
        ),
        by = c("SAMPLE_ID",
               "VAR_ID",
               "VARIANT_CLASS",
               "FUSION_GENE",
               "FUSION_GENE2",
               "BREAKPOINT_5P",
               "BREAKPOINT_3P")
      ) |>
      dplyr::select(
        dplyr::any_of(
          c("SAMPLE_ID",
            "VAR_ID",
            "VARIANT_CLASS",
            "FUSION_GENE",
            "FUSION_GENE2",
            "BREAKPOINT_5P",
            "BREAKPOINT_3P",
            "FUSION_GENE_5P",
            "FUSION_GENE_3P"
          )
        ),
        dplyr::everything()) |>
      dplyr::mutate(num_sort = as.numeric(
        stringr::str_extract(.data$VAR_ID, "\\d+"))) |>
      dplyr::arrange(.data$num_sort) |>
      dplyr::select(-c("num_sort"))


    if (NROW(results[['variant']]) == 0) {
      log4r_info(
        paste0("No RNA fusion calls with gene/transcript overlap found"))
      return(results)
    }

    ## Add split read support and score if available
    if (NROW(score_splitread_df) > 0 &
       "SPLIT_READS" %in% colnames(score_splitread_df)) {

      results[['variant']] <- results[['variant']] |>
        dplyr::left_join(
          score_splitread_df,
          by = c("SAMPLE_ID","VAR_ID"),
          relationship = "many-to-many"
        ) |>
        dplyr::mutate(
          SPLIT_READS = as.integer(.data$SPLIT_READS))

      if ("SCORE" %in% colnames(results[['variant']])) {
        results[['variant']] <- results[['variant']] |>
          dplyr::mutate(SCORE = as.numeric(.data$SCORE))
      }

      ## Filter fusions by minimum split read support
      min_split_reads <-
        settings$conf$rna_fusion$min_split_reads
      if (!is.null(min_split_reads) &
         "SPLIT_READS" %in% colnames(results[['variant']])) {
        n_before <- NROW(results[['variant']])
        results[['variant']] <- results[['variant']] |>
          dplyr::filter(
            .data$SPLIT_READS >= min_split_reads)
        n_filtered <- n_before - NROW(results[['variant']])
        if (n_filtered > 0) {
          log4r_info(
            paste0("Filtered out ", n_filtered,
                   " fusion event(s) with fewer than ",
                   min_split_reads, " split reads"))
        }
      }
    }

    results[['variant_recurrence']] <-
      rna_fusion_recurrence_mitdb(
        query_fusions = results[['variant']],
        ref_data = ref_data
      ) |>
      dplyr::mutate(
        MITDB_NUM_EVIDENCE = dplyr::if_else(
          is.na(.data$MITDB_NUM_EVIDENCE),
          0,
          as.integer(.data$MITDB_NUM_EVIDENCE)
      )) |>
      dplyr::arrange(
        .data$SAMPLE_ID,
        .data$VAR_ID,
        .data$FUSION_GENE2,
        dplyr::desc(.data$MITDB_NUM_EVIDENCE)
      )

    recurrency_data <- list()

    recurrency_data[['variant']] <- as.data.frame(
      results[['variant_recurrence']] |>
        dplyr::group_by(
          .data$SAMPLE_ID,
          .data$VAR_ID,
          .data$FUSION_GENE2
        ) |>
        dplyr::summarise(
          MITDB_NUM_EVIDENCE =
            sum(.data$MITDB_NUM_EVIDENCE),
          MITDB_EVIDENCE = paste(
            .data$MITDB_EVIDENCE2, collapse=" | "
          ), .groups = "drop"
        )
    )

    recurrency_data[['variant_display']] <- as.data.frame(
      results[['variant_recurrence']] |>
        dplyr::group_by(
          .data$SAMPLE_ID,
          .data$VAR_ID,
          .data$FUSION_GENE2
        ) |>
        dplyr::summarise(
          MITDB_NUM_EVIDENCE = sum(.data$MITDB_NUM_EVIDENCE),
          MITDB_EVIDENCE = paste(
            .data$MITDB_EVIDENCE, collapse=", "
          ), .groups = "drop"
        ) |>
        dplyr::mutate(
          MITDB_NUM_EVIDENCE = dplyr::if_else(
            is.na(.data$MITDB_NUM_EVIDENCE),
            0,
            as.integer(.data$MITDB_NUM_EVIDENCE)
          ))
    )

    results[['variant']] <- results[['variant']] |>
      dplyr::mutate(
        SAMPLE_ALTERATION = paste0(
          .data$FUSION_GENE, " fusion"
        )
      ) |>
      dplyr::left_join(
        recurrency_data[['variant']],
        by = c("SAMPLE_ID","VAR_ID","FUSION_GENE2"),
        relationship = "many-to-one"
      ) |>
      dplyr::mutate(
        MITDB_NUM_EVIDENCE = dplyr::if_else(
          is.na(.data$MITDB_NUM_EVIDENCE),
          0,
          as.integer(.data$MITDB_NUM_EVIDENCE)
        )) |>
      dplyr::mutate(ENTREZGENE = dplyr::case_when(
        !is.na(.data$ENTREZGENE_5P)  &
          !is.na(.data$ENTREZGENE_3P) ~ paste(
            .data$ENTREZGENE_5P, .data$ENTREZGENE_3P, sep = "::"),
        !is.na(.data$ENTREZGENE_5P) &
          is.na(.data$ENTREZGENE_3P) ~ paste(
            .data$ENTREZGENE_5P,".", sep = "::"),
        is.na(.data$ENTREZGENE_5P) &
          !is.na(.data$ENTREZGENE_3P) ~ paste(
            ".", .data$ENTREZGENE_3P, sep = "::"),
        TRUE ~ NA_character_
      )) |>
      get_druggable_fusion_partner(
        partner = "5P",
        ref_data = ref_data,
        variant_display = FALSE,
        primary_site = primary_site) |>
      get_druggable_fusion_partner(
        partner = "3P",
        ref_data = ref_data,
        variant_display = FALSE,
        primary_site = primary_site) |>
      dplyr::select(
        dplyr::any_of(
          c("SAMPLE_ID",
            "VAR_ID",
            "VARIANT_CLASS",
            "SAMPLE_ALTERATION",
            "ENTREZGENE",
            "FUSION_GENE",
            "FUSION_GENE2",
            "BREAKPOINT_5P",
            "BREAKPOINT_3P",
            "SPLIT_READS",
            "SCORE",
            "FUSION_GENE_5P",
            "FUSION_GENE_3P"
          )
        ),
        dplyr::everything())



    if(settings$conf$oncokb$run == 1){
      if(!is.null(settings$molecular_data$fname_oncokb_output_fusions) &
         file.exists(settings$molecular_data$fname_oncokb_output_fusions)){


        log4r_info(
          "Querying OncoKB web API with OncoKB-annotated fusion file")

        oncokb_results <-
          process_oncokb_fusion(
            fusion_file =
              settings$molecular_data$fname_oncokb_output_fusions,
            oncokb_token =
              settings$conf$oncokb$api_token,
            oncotree_code =
              settings$conf$oncokb$oncotree_code,
            var_calls = results[['variant']])

        if(NROW(oncokb_results$eitems) > 0){
          results$bm_evidence$eitems <- dplyr::bind_rows(
            results$bm_evidence$eitems,
            oncokb_results$eitems)
        }
      }
    }


    ## assign fusion events to tiers of significance (AMP/ASCO/CAP)
    etype_for_tiering <- c("predictive")

    log4r_info(
      paste0("Variant tier classification",
             " - somatic actionability guidelines (AMP/ASCO/CAP)"))

    amp_asco_cap_classified_variant <- list()
    for (clnsig in names(bm_categories)) {
      etype_for_tiering <- bm_categories[[clnsig]]$etype
      amp_asco_cap_classified_variant[[clnsig]] <-
        assign_amp_asco_cap_tiers(
          vartype = "fusion",
          var_df = results$variant,
          primary_site = primary_site,
          clinical_significance = clnsig,
          biomarker_mapping_confidence = "medium",
          biomarker_items =
            results$bm_evidence$eitems
        )
    }

    for (clnsig in names(bm_categories)) {
      for (elem in c("classification",
                    "eitems")) {
        results$bm_evidence[[clnsig]][[elem]] <-
          amp_asco_cap_classified_variant[[clnsig]]$bm_evidence[[elem]]
      }
      if (clnsig == "therapeutic_sensitivity") {
        if (NROW(amp_asco_cap_classified_variant[[clnsig]][['variant']]) > 0) {
          results[['variant']] <-
            amp_asco_cap_classified_variant[[clnsig]][['variant']]
          results[['bm_evidence']][['classification']] <-
            amp_asco_cap_classified_variant[[clnsig]][['classification']]
        }
        ## issue warning if no variants are classified
        if (NROW(amp_asco_cap_classified_variant[[clnsig]][['variant']]) == 0) {
          log4r_warn(
            paste0(
              "No variants were classified with AMP/ASCO/CAP tiers of clinical",
              "significance for therapeutic sensitivity"))
        }
      }
    }

    OTP_GENE_URL <- variant_db_url[
      which(variant_db_url[,"name"] == "GENENAME"),"url_prefix"]

    ENSEMBL_GENE_URL <- variant_db_url[
      which(variant_db_url[,"name"] == "ENSEMBL_GENE_ID"),"url_prefix"]

    results[['variant_display']] <- results[['variant']] |>
      dplyr::select(
        -dplyr::any_of(
          c("MITDB_NUM_EVIDENCE",
            "MITDB_EVIDENCE",
            "TARGETED_INHIBITORS_ALL2_5P",
            "TARGETED_INHIBITORS_ALL2_3P",
            "TARGETED_INHIBITORS2_5P",
            "TARGETED_INHIBITORS2_3P")
        )
      ) |>
      dplyr::left_join(
        recurrency_data[['variant_display']],
        by = c("SAMPLE_ID","VAR_ID","FUSION_GENE2"),
        relationship = "many-to-one"
      ) |>
      get_druggable_fusion_partner(
        partner = "5P",
        ref_data = ref_data,
        variant_display = TRUE,
        primary_site = primary_site) |>
      get_druggable_fusion_partner(
        partner = "3P",
        ref_data = ref_data,
        variant_display = TRUE,
        primary_site = primary_site) |>
      dplyr::mutate(
        FUSION_GENE_5P = dplyr::case_when(
          !is.na(.data$ENTREZGENE_5P) ~ glue::glue(
            "<a href='{OTP_GENE_URL}{ENTREZGENE_5P}'",
            " target='_blank'>{FUSION_GENE_5P}</a>"
          ),
          is.na(.data$ENTREZGENE_5P) &
            !is.na(.data$ENSEMBL_GENE_ID_5P) ~ glue::glue(
              "<a href='{ENSEMBL_GENE_URL}{ENSEMBL_GENE_ID_5P}'",
              " target='_blank'>{FUSION_GENE_5P}</a>"
            ),
          TRUE ~ as.character(FUSION_GENE_5P)
        )) |>
      dplyr::mutate(
        FUSION_GENE_3P = dplyr::case_when(
          !is.na(.data$ENTREZGENE_3P) ~ glue::glue(
            "<a href='{OTP_GENE_URL}{ENTREZGENE_3P}'",
            " target='_blank'>{FUSION_GENE_3P}</a>"
          ),
          is.na(.data$ENTREZGENE_3P) &
            !is.na(.data$ENSEMBL_GENE_ID_3P) ~ glue::glue(
              "<a href='{ENSEMBL_GENE_URL}{ENSEMBL_GENE_ID_3P}'",
              " target='_blank'>{FUSION_GENE_3P}</a>"
            ),
          TRUE ~ as.character(FUSION_GENE_3P)
        )) |>
      dplyr::mutate(
        BREAKPOINT_5P = glue::glue(
          "<a href='https://genome.ucsc.edu/cgi-bin/hgTracks?db={hgname}",
          "&position=chr{BREAKPOINT_5P}' target='_blank'>chr{BREAKPOINT_5P}</a>"
        ),
        BREAKPOINT_3P = glue::glue(
          "<a href='https://genome.ucsc.edu/cgi-bin/hgTracks?db={hgname}",
          "&position=chr{BREAKPOINT_3P}' target='_blank'>chr{BREAKPOINT_3P}</a>"
        )) |>
      dplyr::select(
        dplyr::any_of(
          c("VAR_ID",
            "VARIANT_CLASS",
            "ENTREZGENE",
            "FUSION_GENE",
            "FUSION_GENE2",
            "SPLIT_READS",
            "SCORE",
            "FUSION_GENE_5P",
            "FUSION_GENE_3P",
            "BREAKPOINT_5P",
            "BREAKPOINT_3P",
            "GENENAME_5P",
            "ONCOGENE_5P",
            "ENSEMBL_TRANSCRIPT_ID_5P",
            "TARGETED_INHIBITORS_5P",
            "TARGETED_INHIBITORS_ALL_5P",
            "GENENAME_3P",
            "ONCOGENE_3P",
            "ENSEMBL_TRANSCRIPT_ID_3P",
            "TARGETED_INHIBITORS_3P",
            "TARGETED_INHIBITORS_ALL_3P",
            "SAMPLE_ALTERATION",
            "MITDB_NUM_EVIDENCE",
            "MITDB_EVIDENCE",
            "ACTIONABILITY_TIER",
            "ACTIONABILITY"
          )
        )
      ) |>
      dplyr::distinct()


  }else{
    log4r_info(
      paste0("No RNA fusion calls with gene/transcript overlap found"))
    return(results)
  }

  return(results)

}
