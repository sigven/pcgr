#' Function that reads and validates a fully annotated CNA file from PCGR
#' pre-report pipeline
#'
#' @param fname Path to file name
#' @param ref_data Object with reference data
#' @param settings Object with PCGR report configuration
#'
#' @export
load_somatic_cna <- function(
    fname,
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

  callset_cna <- pcgrr::load_dna_variants(
    fname = fname,
    cols = pcgrr::data_coltype_defs$cna_somatic_raw,
    ref_data = ref_data,
    vartype = 'cna',
    primary_site =
      tumor_site,
    retained_info_tags = "None",
    variant_origin = "Somatic")


  if (NROW(callset_cna$variant) > 0) {
    callset_cna[['variant']] <- callset_cna[['variant']] |>
      dplyr::mutate(CN_TOTAL =
                      as.integer(.data$N_MAJOR + .data$N_MINOR)) |>
      dplyr::rename(CN_MINOR = "N_MINOR",
                    CN_MAJOR = "N_MAJOR") |>
      dplyr::mutate(MOLECULAR_ALTERATION = paste(
        .data$SYMBOL, .data$VARIANT_CLASS)) |>
      pcgrr::append_cancer_association_ranks(
        ref_data = ref_data,
        primary_site = tumor_site,
        pos_var = 'SEGMENT_START') |>
      pcgrr::append_targeted_drug_annotations(
        ref_data = ref_data,
        primary_site = tumor_site) |>
      dplyr::arrange(
        .data$TIER,
        dplyr::desc(.data$TISSUE_ASSOC_RANK),
        dplyr::desc(.data$GLOBAL_ASSOC_RANK))

    pcgrr::log4r_info("Generating data frame with hyperlinked variant/gene annotations")

    callset_cna[['variant_display']] <- callset_cna[['variant']] |>

      dplyr::mutate(
        SEGMENT = glue::glue(
          "<a href='http://genome.ucsc.edu/cgi-bin/hgTracks?db={hgname}&position=",
          "chr{CHROM}:{SEGMENT_START}-{SEGMENT_END}' target='_blank'>",
          "chr{CHROM}:{SEGMENT_START}-{SEGMENT_END}</a>"
        )
      ) |>
      pcgrr::append_cancer_gene_evidence(
        ref_data = ref_data,
        pos_var = 'SEGMENT_START'
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
        TRANSCRIPT_OVERLAP,
        c("ENSEMBL_TRANSCRIPT_ID",
          "REFSEQ_TRANSCRIPT_ID",
          "TRANSCRIPT_OVERLAP_PERCENT"),
        sep="|", remove = T
      ) |>
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
      dplyr::rename(
        TARGETED_INHIBITORS = "TARGETED_INHIBITORS2",
        TARGETED_INHIBITORS_ALL = "TARGETED_INHIBITORS_ALL2"
      ) |>

      dplyr::arrange(
        .data$TIER,
        dplyr::desc(.data$GLOBAL_ASSOC_RANK),
        dplyr::desc(.data$TISSUE_ASSOC_RANK))
  }

  return(callset_cna)

}

#' Function that reads and validates an annotated somatic SNV/InDel
#' file from PCGR pre-reporting pipeline
#'
#' @param fname Path to file name
#' @param ref_data Object with reference data
#' @param settings Object with PCGR report configuration
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
      primary_site = tumor_site,
      pos_var = 'POS') |>
    pcgrr::append_targeted_drug_annotations(
         ref_data = ref_data,
         primary_site = tumor_site) |>
    tidyr::separate(.data$HGVSc, c("ENST", "tmp_HGVSc"),
                    sep = ":", remove = F) |>
    dplyr::mutate(
      MOLECULAR_ALTERATION = dplyr::case_when(
        is.na(.data$SYMBOL) &
          !is.na(.data$CONSEQUENCE) ~
          as.character(.data$CONSEQUENCE),
        !is.na(.data$CONSEQUENCE) &
        stringr::str_detect(
          .data$CONSEQUENCE, "^(splice_acceptor|splice_donor)") &
          !is.na(.data$SYMBOL) &
          !is.na(.data$tmp_HGVSc) ~
        paste0(.data$SYMBOL," ",
              .data$CONSEQUENCE, " - ",
              .data$tmp_HGVSc),
        .data$EXONIC_STATUS == "exonic" &
          !is.na(.data$CONSEQUENCE) &
          !is.na(.data$HGVSP) ~
          paste0(.data$SYMBOL," ",
                .data$CONSEQUENCE, " - ",
                .data$HGVSP),
        TRUE ~ as.character(paste0(
          .data$SYMBOL," ",.data$CONSEQUENCE))
        )) |>
    dplyr::select(-c("tmp_HGVSc","ENST"))


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
          !is.na(.data$TIER) &
            .data$TIER <= 2 &
            .data$SOMATIC_CLASSIFICATION != "SOMATIC") |>
        NROW()

      if(n_actionable_filtered > 0){
        pcgrr::log4r_warn(
          paste0(
            "A total of n = ", n_actionable_filtered,
            " clinically actionable ",
            "variants were filtered as likely germline events"))
      }

      pcgrr::log4r_info(
        paste0(
          "Excluded n = ",
          NROW(callset$variant_unfiltered) - NROW(callset$variant),
          " variants aftering filtering "))

      if(as.logical(
        settings$conf$somatic_snv$tumor_only[["exclude_nonexonic"]]) == TRUE &
        NROW(callset[['variant']]) > 0 &
        "EXONIC_STATUS" %in% colnames(callset[['variant']])){

        callset[['variant']] <- callset[['variant']] |>
          dplyr::filter(
            .data$EXONIC_STATUS == "exonic")
      }

    }else{
      pcgrr::log4r_fatal(
        "Variant data.frame is lacking a 'SOMATIC_CLASSIFICATION' column")
    }
  }

  if(NROW(callset[['variant']]) > 0){
    callset[['variant']] <- callset[['variant']] |>
      dplyr::arrange(
        .data$TIER,
        dplyr::desc(.data$GLOBAL_ASSOC_RANK),
        dplyr::desc(.data$TISSUE_ASSOC_RANK),
        dplyr::desc(.data$ONCOGENICITY_SCORE))

    ## Make data frame with columns for display
    ## in HTML output

    pcgrr::log4r_info("Generating data frame with hyperlinked variant/gene annotations")


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
          .data$MUTATION_HOTSPOT_CANCERTYPE,",",", "))
  }

  # callset <-
  #   pcgrr::expand_biomarker_items(
  #     callset = callset,
  #     variant_origin = "somatic")

  return(callset)


}

#' Function that reads and validates CNA or SNV/InDel TSV files
#' file from PCGR/CPSR pre-report (Python) pipeline
#'
#' @param fname Path to raw input file with DNA aberrations (PCGR/CPSR)
#' @param cols column type definitions of raw input file
#' @param ref_data reference data object
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
    tibble::has_name(calls_raw, compulsary_cols)
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

  if ("VAF_CONTROL" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(
        VAF_CONTROL = round(as.numeric(.data$VAF_CONTROL), 3)
      )
  }

  if ("VAF_TUMOR" %in% colnames(results[['variant']])) {
    results[['variant']] <-
      results[['variant']] |>
      dplyr::mutate(
        VAF_TUMOR = round(as.numeric(.data$VAF_TUMOR), 3)
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
        -c("ENSEMBL_TRANSCRIPT_ID",
           "REFSEQ_TRANSCRIPT_ID",
           "TRANSCRIPT_START",
           "TRANSCRIPT_END")) |>
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
      )

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
              by = c("VARIANT_ID","BIOMARKER_SOURCE")) |>
            dplyr::rename(BIOMARKER_MATCH = .data$BIOMARKER_MATCHTYPE) |>
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
              BM_VARIANT_ID = .data$VARIANT_ID,
              BM_EVIDENCE_ID = .data$EVIDENCE_ID,
              BM_SOURCE_DB = .data$BIOMARKER_SOURCE,
              BM_RESOLUTION = .data$BIOMARKER_RESOLUTION,
              BM_MATCH = .data$BIOMARKER_MATCH,
              BM_PRIMARY_SITE = .data$PRIMARY_SITE,
              BM_EVIDENCE_TYPE = .data$EVIDENCE_TYPE,
              BM_CANCER_TYPE = .data$CANCER_TYPE,
              BM_DISEASE_ONTOLOGY_ID = .data$DISEASE_ONTOLOGY_ID,
              BM_VARIANT_ORIGIN = .data$VARIANT_ORIGIN,
              BM_EVIDENCE_LEVEL = .data$EVIDENCE_LEVEL,
              BM_EVIDENCE_DESCRIPTION = .data$EVIDENCE_DESCRIPTION,
              BM_THERAPEUTIC_CONTEXT = .data$THERAPEUTIC_CONTEXT,
              BM_CLINICAL_SIGNIFICANCE = .data$CLINICAL_SIGNIFICANCE,
              BM_CITATION = .data$CITATION,
              BM_REFERENCE = .data$CITATION_HTML,
              BM_RATING = .data$RATING,
              BM_EVIDENCE_DIRECTION = .data$EVIDENCE_DIRECTION,
              BM_MOLECULAR_PROFILE_NAME = .data$MOLECULAR_PROFILE_NAME,
              BM_MOLECULAR_PROFILE_TYPE = .data$MOLECULAR_PROFILE_TYPE
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
                  .data$BM_MOLECULAR_PROFILE_TYPE == "Any") |>
            dplyr::distinct() |>
            dplyr::mutate(
              BM_MOLECULAR_PROFILE_NAME = dplyr::if_else(
                BM_SOURCE_DB == "cgi",
                paste0(
                  "<a href='https://www.cancergenomeinterpreter.org/biomarkers'",
                  " target='_blank'>",
                  stringr::str_replace_all(
                    .data$BM_MOLECULAR_PROFILE_NAME,
                    ",",", "),
                  "</a>"
                ),
                paste0(
                  "<a href='https://civicdb.org/evidence/",
                  .data$BM_EVIDENCE_ID,
                  "' target='_blank'>",
                  stringr::str_replace_all(
                    .data$BM_MOLECULAR_PROFILE_NAME,
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
    ## ACMG/AMP guidelines for clinical actionable
    ## of somatic variants
    if (variant_origin == "Somatic") {

      pcgrr::log4r_info(
        paste0("Assigning variants to tiers of clinical significance",
               " - ACMG/AMP guidelines"))
      results <- pcgrr::assign_acmg_tiers(
        vartype = vartype,
        variants_df = results$variant,
        primary_site = primary_site,
        biomarker_items =
          results$biomarker_evidence$items
      )
    }

  }else{
    log4r_fatal("Input data does not contain 'BIOMARKER_MATCH' column - fatal")
  }

  results[['retained_info_tags']] <- paste(
    retained_cols_renamed, collapse=","
  )

  return(results)

}

load_expression_sim <- function(ref_data = NULL,
                                settings = NULL){

  ## Load expression similarity results for input sample
  ## against reference collections
  expression_sim <- list()

  if(!is.null(settings$conf$gene_expression$similarity_db)){
    for(db in names(settings$conf$gene_expression$similarity_db)){
      fname_db <- paste0("fname_expression_sim_",db)
      if(!is.null(settings$molecular_data[[fname_db]]) &
        file.exists(settings$molecular_data[[fname_db]])){
          expression_sim[[db]] <- readr::read_tsv(
            settings$molecular_data[[fname_db]],
            show_col_types = F, na = ".")
      }
    }
  }
  return(expression_sim)
}
