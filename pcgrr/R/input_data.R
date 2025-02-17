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

  pcgrr::log4r_info("------")
  pcgrr::log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - somatic copy number aberrations"))

  hgname <- "hg38"
  if(settings$genome_assembly == "grch37"){
    hgname <- "hg19"
  }

  tumor_site <-
    settings[['conf']][['sample_properties']][['site']]

  ## read segments
  pcgrr::check_file_exists(fname_cna_segment)
  pcgrr::check_file_exists(fname_cna_gene)
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
    names(pcgrr::data_coltype_defs$cna_somatic_segment_raw$cols)

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
          pcgrr::data_coltype_defs$cna_somatic_segment_raw$cols,
        na = c(".","NA"),
        progress = F
      )
    )) |>
    tidyr::separate(
      col = "SEGMENT_NAME",
      into = c("SEGMENT_ID", "N_MAJOR","N_MINOR","ARM","CYTOBAND","EVENT_TYPE"),
      sep = "\\|",
      remove = T
    ) |>
    dplyr::mutate(
      CN_TOTAL = as.integer(
        as.integer(.data$N_MAJOR) + as.integer(.data$N_MINOR))
    ) |>
    dplyr::rename(CN_MINOR = "N_MINOR",
                  CN_MAJOR = "N_MAJOR") |>
    dplyr::mutate(CN_MINOR = as.integer(.data$CN_MINOR),
                  CN_MAJOR = as.integer(.data$CN_MAJOR))


  callset_cna <- pcgrr::load_dna_variants(
    fname = fname_cna_gene,
    cols = pcgrr::data_coltype_defs$cna_somatic_gene_raw,
    ref_data = ref_data,
    vartype = 'cna',
    primary_site =
      tumor_site,
    retained_info_tags = "None",
    variant_origin = "Somatic")

  callset_cna[['segment']] <- segments

  if (NROW(callset_cna$variant) > 0) {
    callset_cna[['variant']] <- callset_cna[['variant']] |>
      dplyr::mutate(CN_TOTAL =
                      as.integer(.data$N_MAJOR + .data$N_MINOR)) |>
      dplyr::rename(CN_MINOR = "N_MINOR",
                    CN_MAJOR = "N_MAJOR") |>
      dplyr::mutate(MOLECULAR_ALTERATION = paste(
        .data$SYMBOL, .data$VARIANT_CLASS)) |>
      dplyr::mutate(CYTOBAND = paste0(
        "chr", .data$CHROM,":" ,.data$CYTOBAND
      )) |>
      pcgrr::append_cancer_association_ranks(
        ref_data = ref_data,
        primary_site = tumor_site) |>
      pcgrr::append_targeted_drug_annotations(
        ref_data = ref_data,
        primary_site = tumor_site) |>
      dplyr::arrange(
        .data$ACTIONABILITY_TIER,
        dplyr::desc(.data$TISSUE_ASSOC_RANK),
        dplyr::desc(.data$GLOBAL_ASSOC_RANK)) |>
      pcgrr::order_variants(pos_var = 'SEGMENT_START') |>
      pcgrr::exclude_non_chrom_variants()

    pcgrr::log4r_info(
      "Generating data frame with hyperlinked variant/gene annotations")

    callset_cna[['variant_display']] <- callset_cna[['variant']] |>
      dplyr::mutate(
        SEGMENT = glue::glue(
          "<a href='http://genome.ucsc.edu/cgi-bin/hgTracks?db={hgname}&position=",
          "chr{CHROM}:{SEGMENT_START}-{SEGMENT_END}' target='_blank'>",
          "chr{CHROM}:{SEGMENT_START}-{SEGMENT_END}</a>"
        )
      ) |>
      pcgrr::append_cancer_gene_evidence(
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
      pcgrr::append_annotation_links(
        vartype = "cna"
      ) |>
      tidyr::unite(
        "TRANSCRIPT_OVERLAP",
        c("ENSEMBL_TRANSCRIPT_ID",
          "REFSEQ_TRANSCRIPT_ID",
          "TRANSCRIPT_OVERLAP_PERCENT"),
        sep="|", remove = T
      ) |>
      pcgrr::order_variants(pos_var = 'SEGMENT_START') |>
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
#'
#' @export
load_somatic_snv_indel <- function(
    fname = NA,
    ref_data = NULL,
    settings = NULL) {

  pcgrr::log4r_info("------")
  pcgrr::log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - somatic SNV/InDels"))

  tumor_site <-
    settings[['conf']][['sample_properties']][['site']]

  callset <- pcgrr::load_dna_variants(
    fname = fname,
    cols = pcgrr::data_coltype_defs$snv_indel_somatic_raw,
    ref_data = ref_data,
    vartype = 'snv_indel',
    primary_site = tumor_site,
    retained_info_tags =
      settings[['conf']][['other']]$retained_vcf_info_tags,
    variant_origin = "Somatic")


  callset[['variant_unfiltered']] <- data.frame()
  #pcgrr::log4r_info("Adding hyperlinks to variant/gene annotations")

  callset[['variant']] <- callset[['variant']] |>
    pcgrr::append_cancer_association_ranks(
      ref_data = ref_data,
      primary_site = tumor_site) |>
    pcgrr::append_targeted_drug_annotations(
         ref_data = ref_data,
         primary_site = tumor_site) |>
    dplyr::mutate(
      MOLECULAR_ALTERATION = dplyr::case_when(
        is.na(.data$SYMBOL) &
          !is.na(.data$CONSEQUENCE) ~
          as.character(.data$CONSEQUENCE),
        !is.na(.data$CONSEQUENCE) &
        stringr::str_detect(
          .data$CONSEQUENCE, "^(splice_acceptor|splice_donor)") &
          !is.na(.data$SYMBOL) &
          !is.na(.data$HGVSc) ~
        paste0(.data$SYMBOL," ",
              stringr::str_replace_all(
                .data$CONSEQUENCE,"&",", "),
              " - ",
              .data$HGVSc),
        .data$EXONIC_STATUS == "exonic" &
          !stringr::str_detect(
            .data$CONSEQUENCE, "^(splice_acceptor|splice_donor)") &
          !is.na(.data$CONSEQUENCE) &
          !is.na(.data$HGVSc) &
          !is.na(.data$SYMBOL) &
          !is.na(.data$HGVSP) ~
          paste0(.data$SYMBOL," ",
                stringr::str_replace_all(
                  .data$CONSEQUENCE, "&",", "),
                  " - ",
                .data$HGVSc,
                " - ",
                .data$HGVSP),
        TRUE ~ as.character(paste0(
          .data$SYMBOL," ",
          stringr::str_replace_all(
            .data$CONSEQUENCE,"&",", "))
        )
      )) |>
    pcgrr::order_variants(pos_var = 'POS') |>
    pcgrr::exclude_non_chrom_variants()


  ## Tumor-only input
  if (as.logical(settings$conf$assay_properties$vcf_tumor_only) == TRUE) {
    callset[['variant_unfiltered']] <-
      callset[['variant']] |>
      ## assign evidence tags for germline/somatic state of variants,
      ## partially based on user-defined options
      ## (population allele frequency thresholds)
        pcgrr::assign_somatic_germline_evidence(
          settings = settings) |>

      ## assign somatic variant classification/status based on accumulation
      ## of evidence tags and user-defined options
        pcgrr::assign_somatic_classification(
          settings = settings)

    ## Assign calls to filtered callset (SOMATIC_CLASSIFICATION = SOMATIC)
    if("SOMATIC_CLASSIFICATION" %in%
       colnames(callset[['variant_unfiltered']])){
      callset[['variant']] <-
        callset[['variant_unfiltered']] |>
        dplyr::filter(
          .data$SOMATIC_CLASSIFICATION == "SOMATIC")

      ## Issue warning if clinically actionable variants are filtered
      ## with current filtering settings
      n_actionable_filtered <-
        callset[['variant_unfiltered']] |>
        dplyr::filter(
          !is.na(.data$ACTIONABILITY_TIER) &
            .data$ACTIONABILITY_TIER <= 2 &
            .data$SOMATIC_CLASSIFICATION != "SOMATIC") |>
        NROW()

      if(n_actionable_filtered > 0){
        pcgrr::log4r_warn(
          paste0(
            "A total of n = ", n_actionable_filtered,
            " clinically actionable ",
            "variants were filtered as likely germline events"))
      }

      n_excluded_calls <- as.character(
        formattable::comma(
          (NROW(callset$variant_unfiltered) - NROW(callset$variant)),
          digits = 0))

      pcgrr::log4r_info(
        paste0(
          "Excluded n = ",
          n_excluded_calls,
          " variants aftering filtering of putative germline events"))

      if(NROW(callset$variant) == 0){
        pcgrr::log4r_warn(
          "NO (n = 0) somatic variants remain after filtering of putative germline events")
      }

      if(as.logical(
        settings$conf$somatic_snv$tumor_only[["exclude_nonexonic"]]) == TRUE &
        NROW(callset[['variant']]) > 0 &
        "EXONIC_STATUS" %in% colnames(callset[['variant']])){

        callset[['variant']] <- callset[['variant']] |>
          dplyr::filter(
            .data$EXONIC_STATUS == "exonic")

        if(NROW(callset$variant) == 0){
          pcgrr::log4r_warn(
            "NO (n = 0) somatic variants remain after filtering of putative germline events")
        }
      }

      ## filter also MAF file if provided
      pcgrr::filter_maf_file(
        callset = callset,
        settings = settings)

    }else{
      pcgrr::log4r_fatal(
        "Variant data.frame is lacking a 'SOMATIC_CLASSIFICATION' column")
    }
  }

  if(NROW(callset[['variant']]) > 0){
    callset[['variant']] <- callset[['variant']] |>
      dplyr::arrange(
        .data$ACTIONABILITY_TIER,
        dplyr::desc(.data$ONCOGENICITY_SCORE),
        dplyr::desc(.data$GLOBAL_ASSOC_RANK),
        dplyr::desc(.data$TISSUE_ASSOC_RANK))


    ## Make data frame with columns for display
    ## in HTML output

    pcgrr::log4r_info(
      "Generating data frame with hyperlinked variant/gene annotations")


    callset[['variant_display']] <- callset[['variant']] |>
      pcgrr::append_cancer_gene_evidence(
        ref_data = ref_data) |>
      pcgrr::append_dbmts_var_link() |>
      pcgrr::append_tcga_var_link() |>
      pcgrr::append_annotation_links() |>
      pcgrr::append_dbnsfp_var_link() |>
      pcgrr::append_tfbs_annotation() |>
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
      pcgrr::order_variants(pos_var = 'POS')

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
#'
#' @export
#'
load_cpsr_classified_variants <- function(
    fname_cpsr_tsv = NA,
    fname_cpsr_yaml = NA,
    cols = NULL,
    ignore_vus = FALSE,
    ref_data = NULL){

  pcgrr::log4r_info("------")
  pcgrr::log4r_info(paste0(
    "Reading annotated molecular dataset (DNA) - germline SNV/InDels (CPSR-classified)"))

  pcgrr::check_file_exists(fname_cpsr_tsv)
  pcgrr::check_file_exists(fname_cpsr_yaml)

  if (!file.exists(fname_cpsr_yaml)) {
    log4r_fatal(
      paste0("YAML file '", fname_cpsr_yaml, "' does not exist - exiting"))
  }
  cpsr_yaml <- yaml::read_yaml(fname_cpsr_yaml)
  if("conf" %in% names(cpsr_yaml) == FALSE){
    pcgrr::log4r_fatal(
      paste0(
        "YAML file '", fname_cpsr_yaml,
        "' does not contain a 'conf' section - exiting"))
  }
  if("sample_id" %in% names(cpsr_yaml) == FALSE){
    pcgrr::log4r_fatal(
      paste0(
        "YAML file '", fname_cpsr_yaml,
        "' does not contain a 'sample_id' variable - exiting"))
  }
  if("gene_panel" %in% names(cpsr_yaml$conf) == FALSE){
    pcgrr::log4r_fatal(
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

  callset <- pcgrr::load_dna_variants(
    fname = fname_cpsr_tsv,
    cols = cols,
    ref_data = ref_data,
    vartype = 'snv_indel',
    variant_origin = 'Germline')

  callset[['variant_display']] <- callset[['variant']] |>
    pcgrr::append_cancer_gene_evidence(
      ref_data = ref_data) |>
    dplyr::mutate(
      CLINVAR_TRAITS_ALL = paste(
        stringr::str_to_title(.data$CLINVAR_VARIANT_ORIGIN),
        .data$CLINVAR_PHENOTYPE,
        sep = " - ")) |>
    pcgrr::append_annotation_links() |>
    dplyr::select(
      -dplyr::contains("_RAW")
    ) |>
    dplyr::mutate(
      CONSEQUENCE = stringr::str_replace_all(
        .data$CONSEQUENCE,"&",", ")) |>
    dplyr::rename(
      SOURCE = .data$CPSR_CLASSIFICATION_SOURCE,
      CLINICAL_SIGNIFICANCE = .data$FINAL_CLASSIFICATION
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
         "VAR_ID",
         "ENTREZGENE"))) |>
    dplyr::select(
      dplyr::any_of(
        c("SYMBOL","ALTERATION","GENOTYPE","CONSEQUENCE",
        "CLINICAL_SIGNIFICANCE","SOURCE","PROTEIN_DOMAIN",
        "HGVSc", "HGVSc_RefSeq", "HGVSp", "CDS_CHANGE",
        "CODING_STATUS",
        "LOSS_OF_FUNCTION", "DP_CONTROL",
        "VARIANT_CLASS","GENENAME",
        "ONCOGENE","TUMOR_SUPPRESSOR","ENSEMBL_GENE_ID",
        "ENSEMBL_TRANSCRIPT_ID","REFSEQ_TRANSCRIPT_ID",
        "DBSNP_RSID",
        "CLINVAR","CLINVAR_CLASSIFICATION","CLINVAR_CONFLICTED",
        "CLINVAR_REVIEW_STATUS_STARS",
        "CPSR_PATHOGENICITY_SCORE",
        "CPSR_CLASSIFICATION",
        "CPSR_CLASSIFICATION_CODE")),
      dplyr::everything()
    )

  if(NROW(callset[['variant_display']]) > 0){
    callset[['variant_display']] <- callset[['variant_display']] |>
      dplyr::filter((.data$CLINICAL_SIGNIFICANCE == "Pathogenic" |
                      .data$CLINICAL_SIGNIFICANCE == "Likely_Pathogenic" |
                      .data$CLINICAL_SIGNIFICANCE == "VUS") &
                      .data$CODING_STATUS == "coding")|>
      dplyr::distinct()
    if(ignore_vus == TRUE){
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
    vartype = 'snv_indel',
    primary_site = "Any",
    retained_info_tags = "None",
    variant_origin = "Somatic") {

  pcgrr::check_file_exists(fname)
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
  if (retained_info_tags != "None") {
    retained_cols <- stringr::str_split(
      retained_info_tags, pattern = ",")[[1]]
    for(c in retained_cols) {
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
    for(c in retained_cols) {
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

  results <- list()
  results[['variant']] <- calls
  results[['biomarker_evidence']] <- list()
  results[['biomarker_evidence']][['items']] <-
    data.frame()

  ## Rename annotations for more clarity
  if ("TSG" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        TUMOR_SUPPRESSOR = "TSG"
      )
  }

  ## Rename annotations for more clarity
  if ("TSG_SUPPORT" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        TUMOR_SUPPRESSOR_SUPPORT = "TSG_SUPPORT"
      )
  }

  if ("VEP_ALL_CSQ" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(
        VEP_ALL_CSQ = stringr::str_replace_all(
          .data$VEP_ALL_CSQ, ",",", "
        )
      )
  }

  if ("HGVSp_short" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        HGVSP = "HGVSp_short"
      )
  }

  if ("TSG_RANK" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::rename(
        TUMOR_SUPPRESSOR_RANK = "TSG_RANK"
      )
  }

  for(col in c('VAF_TUMOR','VAF_CONTROL','TPM',
               'TPM_GENE','TPM_MIN','consTPM')) {
    if (col %in% colnames(results[['variant']])) {
      significant_digits = 3
      if(col == 'TPM' | col == 'TPM_GENE' |
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
  if("TPM_GENE" %in% colnames(results[['variant']]) &
       "TPM" %in% colnames(results[['variant']])){
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

  if ("TPM" %in% colnames(results[['variant']])){
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
      "TPM" %in% colnames(results[['variant']])){
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

  if ("BIOMARKER_MATCH" %in% colnames(results[['variant']]) &
     "VAR_ID" %in% colnames(results[['variant']])) {

    ## find all sample variants that have been matched against
    ## biomarkers - (genomic/hgvsp/codon/exon/gene)
    biomarker_set <-
      results[['variant']] |>
      dplyr::filter(!is.na(.data$BIOMARKER_MATCH))

    if (NROW(biomarker_set) > 0) {

      ## get citations of all biomarkers
      citations <- as.data.frame(
        ref_data[['biomarker']][['literature']] |>
          dplyr::select(
            c("EVIDENCE_ID",
              "SOURCE_ID",
              "LINK",
              "NAME")
          ) |>
          dplyr::mutate(
            CITATION = paste0(
              .data$SOURCE_ID,"|",
              .data$NAME
            )
          ) |>
          tidyr::separate_rows(
            c("EVIDENCE_ID"),
            sep=";"
          ) |>
          dplyr::group_by(
            .data$EVIDENCE_ID
          ) |>
          dplyr::summarise(
            CITATION = paste(
              unique(.data$CITATION), collapse = ", "
            ),
            CITATION_HTML = paste(
              unique(.data$LINK), collapse = ", "
            ),
            .groups = "drop"
          )
      )

      ## append more details of matched biomarkers from
      ## reference data on biomarkers (ref_dat$biomarker$variant)
      results[['biomarker_evidence']][['items']] <-
        as.data.frame(
          biomarker_set |>
            dplyr::select(
              c("VAR_ID",
                "VARIANT_CLASS",
                "BIOMARKER_MATCH"),
              ) |>
            dplyr::distinct() |>
            tidyr::separate_rows(
              "BIOMARKER_MATCH", sep=",") |>
            tidyr::separate(
              "BIOMARKER_MATCH",
              into = c("BIOMARKER_SOURCE",
                       "VARIANT_ID",
                       "EVIDENCE_ITEMS",
                       "BIOMARKER_MATCHTYPE"),
              sep = "\\|"
            ) |>
            dplyr::mutate(
              VARIANT_ID = as.character(.data$VARIANT_ID)) |>
            dplyr::left_join(
              dplyr::select(
                ref_data[['biomarker']][['variant']],
                c("VARIANT_ID", "ENTREZGENE","BIOMARKER_SOURCE")),
              by = c("VARIANT_ID","BIOMARKER_SOURCE"),
              relationship = "many-to-many") |>
            dplyr::rename(BIOMARKER_MATCH = BIOMARKER_MATCHTYPE) |>
            dplyr::mutate(BIOMARKER_RESOLUTION = dplyr::case_when(
              stringr::str_detect(.data$BIOMARKER_MATCH,"by_cna_segment") ~ "gene",
              stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") ~ "genomic",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_principal") ~ "hgvsp",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_principal") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_principal") ~ "codon",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_principal") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_principal") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_nonprincipal")~ "hgvsp_nonprincipal",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_principal") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_principal") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_exon_(mut|insertion|deletion)_nonprincipal") ~ "exon_nonprincipal",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_principal") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_principal") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_exon_(mut|insertion|deletion)_principal") ~ "exon",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_exon_") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_aa_region") ~ "gene_region_mut",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_exon_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_aa_") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_gene_mut_lof") ~ "gene_lof",
              !stringr::str_detect(.data$BIOMARKER_MATCH,"by_genomic_coord") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_hgvsp_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_codon_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_exon_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_aa_") &
                !stringr::str_detect(.data$BIOMARKER_MATCH,"by_gene_mut_lof") &
                stringr::str_detect(.data$BIOMARKER_MATCH,"by_gene_mut") ~ "gene_mut",
              TRUE ~ as.character('other')
            )) |>

            #dplyr::select(-)
            tidyr::separate_rows(
              "EVIDENCE_ITEMS", sep="&") |>
            tidyr::separate(
              "EVIDENCE_ITEMS",
              into = c("EVIDENCE_ID",
                       "PRIMARY_SITE",
                       "CLINICAL_SIGNIFICANCE",
                       "EVIDENCE_LEVEL",
                       "EVIDENCE_TYPE",
                       "VARIANT_ORIGIN"),
              sep = ":"
            ) |>
            dplyr::select(
              -c("CLINICAL_SIGNIFICANCE",
                 "EVIDENCE_LEVEL",
                 "EVIDENCE_TYPE")) |>
            dplyr::mutate(PRIMARY_SITE = stringr::str_replace_all(
              .data$PRIMARY_SITE, "_"," "
            )) |>
            dplyr::left_join(
              citations, by = "EVIDENCE_ID",
              relationship = "many-to-many"
            ) |>
            dplyr::left_join(
              dplyr::select(
                ref_data[['biomarker']][['clinical']],
                c("EVIDENCE_ID",
                "CANCER_TYPE",
                "EVIDENCE_TYPE",
                "THERAPEUTIC_CONTEXT",
                "CLINICAL_SIGNIFICANCE",
                "DISEASE_ONTOLOGY_ID",
                "EVIDENCE_LEVEL",
                "EVIDENCE_DESCRIPTION",
                "MOLECULAR_PROFILE_NAME",
                "MOLECULAR_PROFILE_TYPE",
                "RATING",
                "EVIDENCE_DIRECTION")
              ), by = c("EVIDENCE_ID"),
              relationship = "many-to-many"
            ) |>
            dplyr::rename(
              BM_VARIANT_ID = VARIANT_ID,
              BM_EVIDENCE_ID = EVIDENCE_ID,
              BM_SOURCE_DB = BIOMARKER_SOURCE,
              BM_RESOLUTION = BIOMARKER_RESOLUTION,
              BM_MATCH = BIOMARKER_MATCH,
              BM_PRIMARY_SITE = PRIMARY_SITE,
              BM_EVIDENCE_TYPE = EVIDENCE_TYPE,
              BM_CANCER_TYPE = CANCER_TYPE,
              BM_DISEASE_ONTOLOGY_ID = DISEASE_ONTOLOGY_ID,
              BM_VARIANT_ORIGIN = VARIANT_ORIGIN,
              BM_EVIDENCE_LEVEL = EVIDENCE_LEVEL,
              BM_EVIDENCE_DESCRIPTION = EVIDENCE_DESCRIPTION,
              BM_THERAPEUTIC_CONTEXT = THERAPEUTIC_CONTEXT,
              BM_CLINICAL_SIGNIFICANCE = CLINICAL_SIGNIFICANCE,
              BM_CITATION = CITATION,
              BM_REFERENCE = CITATION_HTML,
              BM_RATING = RATING,
              BM_EVIDENCE_DIRECTION = EVIDENCE_DIRECTION,
              BM_MOLECULAR_PROFILE = MOLECULAR_PROFILE_NAME,
              BM_MOLECULAR_PROFILE_TYPE = MOLECULAR_PROFILE_TYPE
            ) |>
            dplyr::mutate(
              BM_RATING = dplyr::if_else(
                is.na(.data$BM_RATING),
                as.integer(2),
                as.integer(.data$BM_RATING)
              )
            ) |>
            dplyr::select(
              c("VAR_ID",
                "VARIANT_CLASS",
                "ENTREZGENE",
                "BM_SOURCE_DB",
                "BM_VARIANT_ID",
                "BM_EVIDENCE_ID",
                "BM_EVIDENCE_TYPE",
                "BM_EVIDENCE_LEVEL",
                "BM_EVIDENCE_DESCRIPTION",
                "BM_EVIDENCE_DIRECTION",
                "BM_CLINICAL_SIGNIFICANCE",
                "BM_VARIANT_ORIGIN",
                "BM_REFERENCE",
                "BM_CANCER_TYPE",
                "BM_DISEASE_ONTOLOGY_ID",
                "BM_PRIMARY_SITE",
                "BM_MATCH",
                "BM_RESOLUTION"),
              dplyr::everything()
            ) |>

            ## Make sure variant origin of patient/tumor
            ## matches reported variant origin of biomarker
            ## - Do not consider complex molecular profile types
            dplyr::filter(
                .data$BM_VARIANT_ORIGIN == variant_origin &
                  .data$BM_MOLECULAR_PROFILE_TYPE == "Single") |>
            dplyr::distinct() |>
            dplyr::mutate(
              BM_MOLECULAR_PROFILE = dplyr::if_else(
                .data$BM_SOURCE_DB == "cgi",
                paste0(
                  "<a href='https://www.cancergenomeinterpreter.org/biomarkers'",
                  " target='_blank'>",
                  stringr::str_replace_all(
                    .data$BM_MOLECULAR_PROFILE,
                    ",",", "),
                  "</a>"
                ),
                paste0(
                  "<a href='https://civicdb.org/evidence/",
                  stringr::str_replace(.data$BM_EVIDENCE_ID,"EID",""),
                  "' target='_blank'>",
                  stringr::str_replace_all(
                    .data$BM_MOLECULAR_PROFILE,
                    ",",", "),
                  "</a>"
                )
              )
            )
        )


      if(variant_origin == "Somatic" &
         vartype == "snv_indel"){

        ## Only consider oncogenic variants
        ## or variants in hotspots for variants
        ## matching biomarkers at the "gene" level
        results$biomarker_evidence$items <-
          results$biomarker_evidence$items |>
          dplyr::inner_join(
              dplyr::select(
                results$variant,
                c("VAR_ID",
                  "ONCOGENICITY",
                  "MUTATION_HOTSPOT",
                  "VARIANT_CLASS",
                  "ENTREZGENE")),
              by = c("VAR_ID","VARIANT_CLASS",
                     "ENTREZGENE")) |>
          dplyr::filter(
            .data$BM_MATCH != "by_gene_mut" |
            (.data$BM_MATCH == "by_gene_mut" &
            stringr::str_detect(
              .data$ONCOGENICITY,"Oncogenic") |
              !is.na(.data$MUTATION_HOTSPOT))) |>
          dplyr::select(
            -c("ONCOGENICITY","MUTATION_HOTSPOT")
          )

      }

    }



    ## Assign each variant a tier according to
    ## AMP/ASCO/CAP guidelines for clinical actionability
    ## of somatic variants
    if (variant_origin == "Somatic") {

      pcgrr::log4r_info(
        paste0("Assigning variants to tiers of clinical significance",
               " - somatic actionability guidelines (AMP/ASCO/CAP)"))
      results <- pcgrr::assign_amp_asco_tiers(
        vartype = vartype,
        variants_df = results$variant,
        primary_site = primary_site,
        biomarker_items =
          results$biomarker_evidence$items
      )
    }

  }
  #else{
    #log4r_fatal("Input data does not contain 'BIOMARKER_MATCH' column - fatal")
  #}

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
load_expression_similarity <- function(settings = NULL){

  ## Load expression similarity results for input sample
  ## against reference collections
  expression_sim <- list()

  if(settings$conf$expression$similarity_analysis == 0){
    return(expression_sim)
  }

  if(file.exists(
    settings$molecular_data$fname_expression_similarity_tsv)){
    pcgrr::log4r_info(
      paste0("Loading expression similarity results for sample ",
             settings$sample_id))

    expression_similarity <- suppressWarnings(readr::read_tsv(
      settings$molecular_data$fname_expression_similarity_tsv,
      show_col_types = F,
      na = ".", guess_max = 100000))

    for(source in unique(expression_similarity$EXT_DB)){
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
    z_score_cutoff = 1.5){

  ## Load expression outlier results for input sample
  ## against reference collections

  tumor_site <- settings$conf$sample_properties$site
  expression_outliers <- data.frame()
  if(file.exists(
    settings$molecular_data$fname_expression_outliers_tsv)){
    pcgrr::log4r_info(
      paste0("Loading expression outlier results for sample ",
             settings$sample_id))

    ## Read raw expression outlier data from Python step of PCGR
    outlier_data <- suppressWarnings(readr::read_tsv(
      settings$molecular_data$fname_expression_outliers_tsv,
      show_col_types = F,
      na = ".", guess_max = 100000))

    if(!is.null(ref_data) &
       !is.null(ref_data$gene) &
       !is.null(ref_data$gene$gene_xref) &
       NROW(outlier_data) > 0 &
       "ENSEMBL_GENE_ID" %in% colnames(outlier_data)){
      outlier_data <-
        outlier_data |>
        dplyr::filter(!is.na(.data$ENSEMBL_GENE_ID))

      ## Append gene annotations from reference data
      ## (SYMBOL, ENTREZGENE, GENENAME, TSG, GENE_BIOTYPE, ONCOGENE)
      if(NROW(outlier_data) > 0){
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


        if(NROW(outlier_data) == 0){
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
          pcgrr::append_cancer_gene_evidence(
            ref_data = ref_data) |>
          pcgrr::append_cancer_association_ranks(
            ref_data = ref_data,
            primary_site = tumor_site) |>
          dplyr::mutate(VAR_ID = dplyr::row_number()) |> ## add unique ID
          pcgrr::append_drug_var_link(
            primary_site = tumor_site,
            ref_data = ref_data) |>
          pcgrr::append_annotation_links(
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

        if(NROW(expression_outliers) == 0){
          return(expression_outliers)
        }

        expression_outliers <- expression_outliers |>
          dplyr::filter(!is.na(.data$EXPR_OUTLIER))

        if(NROW(expression_outliers) == 0){
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

        if(NROW(expression_outliers) == 0){
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
load_expression_csq <- function(settings = NULL){

  ## Load expression consequence results for input sample
  ## against reference collections
  expression_csq <- data.frame()

  if(!is.null(settings$conf$expression$consequence_db)){

    if(!is.null(settings$molecular_data[['fname_csq_expression_tsv']]) &
       file.exists(settings$molecular_data[['fname_csq_expression_tsv']])){
      expression_csq <- readr::read_tsv(
        settings$molecular_data[['fname_csq_expression_tsv']],
        show_col_types = F, na = ".")
    }
  }

  return(expression_csq)

}
